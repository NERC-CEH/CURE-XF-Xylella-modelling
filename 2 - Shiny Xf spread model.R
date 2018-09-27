rm(list=ls())

#setwd("C:/Users/dcha/OneDrive - NERC/CURE-XF/summer school 2018")

# install packages if needed
if(!"shiny" %in% installed.packages()){install.packages("shiny", dep=TRUE)}
if(!"shinyjs" %in% installed.packages()){install.packages("shinyjs", dep=TRUE)}
if(!"data.table" %in% installed.packages()){install.packages("data.table", dep=TRUE)}
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2", dep=TRUE)}
if(!"reshape" %in% installed.packages()){install.packages("reshape", dep=TRUE)}
if(!"viridis" %in% installed.packages()){install.packages("viridis", dep=TRUE)}

# load libraries
library(shiny)
library(shinyjs)
library(data.table)
library(ggplot2)
library(reshape)
library(viridis)

# read in grid cell data.frame...
X = fread("Xf_spread_grid_data.csv", data.table=FALSE)

# re-code zones
X$zone = sapply(gsub("-", "_", as.character(X$zone)), switch,
                containment="buffer",
                disease_free="beyond",
                infected="infected",
                surveillance="buffer")
X$zone = factor(X$zone, levels=c("infected", "buffer", "beyond"))


#######################################################
### function to simulate an Xf epidemic in Puglia, with management in the buffer zone
sim_spread = function(df, r=1, beta=0.4, number_of_years=15, 
                      d50_short = 0.5, d50_long = 15,
                      pLong = 0.0025, 
                      initial_location="fixed",
                      samplesPerCell=0,
                      testFNR=0){
  
  n0 = 1/9106 # initial infection (1 tree per grid cell)
  
  gridRes = 1 # grid cell resolution in km
  
  # add necessary columns to df
  df$initial_olive = df$olive # olive cover at start of infection (to record loss of trees)
  df$Xf = 0 # level of Xylella (initially 0 everywhere)
  df$year_infected = NA # year grid cell is infected (NA if not infected)
  df$year_removed = NA # year trees are removed from grid cell after infection is detected (NA if not removed) 
  
  # initialise the infection at the specified location...
  if(initial_location == "fixed") 
    initXY = c(760743, 4428222)
  if(initial_location == "high") 
    initXY = df[sample(which(df$olive > 0.8 & df$y < 4440000 & df$x < 770000), 1), 1:2]
  if(initial_location == "mid")
    initXY = df[sample(which(df$olive > 0.45 & df$olive < 0.55 & df$y < 4440000 & df$x < 770000), 1), 1:2]
  if(initial_location == "low")
    initXY = df[sample(which(df$olive < 0.2 & df$y < 4440000 & df$x < 770000), 1), 1:2]
  initXY = as.numeric(initXY)
  initPt = which.min(sqrt((df$x-initXY[1])^2 + (df$y-initXY[2])^2))
  df$Xf[initPt] = n0
  df$year_infected[initPt] = 0
  
  df$Xf_yr0 = df$Xf
  
  # run the spread simulation over the specified number of years...
  for(year in 1:number_of_years){
    
    # 1 - Do surveillance and removal...
    # make an index of all infected cells in the buffer zone
    bufferIdx = which(df$Xf>0 & df$zone %in% c("containment", "surveillance")) 
    # simulate sampling and diagnostic testing...
    positiveTests = sapply(bufferIdx, function(i){
      prPos = df$Xf[i]*(1-testFNR) # probability of taking a positive sample
      rbinom(1, size=samplesPerCell, prob=prPos) # Bernoulli trial with that probability to decide if Xf is detected
    })
    removeIdx = bufferIdx[positiveTests > 0] # which cells to remove trees
    df$Xf[removeIdx] = 0 # eradicate the infection
    df$olive[removeIdx] = 0 # clear all the trees
    df$year_removed[removeIdx] = year # record the removal event
    
    # 2 - Grow local infection growth in each grid cell using Gompertz model...
    df$Xf = df$olive * df$Xf^exp(-r)
    
    # 3 - Infection of new grid cells by dispersal...
    infected = df$Xf>0 # Boolean for which grid cells are infected
    notInfected = which(!infected) # integer index for grid cells not infected
    infected = which(infected) # # integer index for grid cells infected
    sourceXf = df$Xf[infected]^beta # infection level in source grid cells
    # calculate probability of new infection...
    prInfection = sapply(notInfected, function(j){
      dij = sqrt((df$x[j] - df$x[infected])^2 + (df$y[j] - df$y[infected])^2)/1000 # distance to all sources (km)
      kShort = exp((dij-gridRes)*log(0.5)/d50_short) # short distance kernel
      #kLong = exp((dij-gridRes)*log(0.5)/d50_long) # long distance kernel
      kLong = (exp(dij*log(0.5)/d50_long)/(2*pi*dij)) / (exp(log(0.5)/d50_long)/(2*pi)) # 2d long distance kernel
      kBoth = ((1-pLong)*kShort + pLong*kLong)  # combined kernel
      pij = kBoth * sourceXf # probability of infection from each individual source
      (1 - exp(sum(log(1 - pij)))) # probability of infection from any source
    })
    # Bernoulli trials to make new infections...
    newInfections = rbinom(n=length(prInfection), size=1, prob=prInfection)
    newInfections = notInfected[newInfections==1] # index of grid cells with new infections
    df$Xf[newInfections] = n0 # set the initial infection level
    df$Xf[df$Xf > df$olive] = df$olive[df$Xf > df$olive] # ensure Xf is never above carrying capacity
    df$year_infected[newInfections] = year # record the year of infection
    
    # save Xf infection in each year
    df = cbind(df, df$Xf)
    names(df)[ncol(df)] = paste0("Xf_yr", year)
    
    # output to screen...
    message("Done year ", year)
  }
  
  # return the data.frame...
  return(df)
}

ui <- fluidPage(
  # *Input() functions
  
  shinyjs::useShinyjs(),
  
  title = "Spatial simulation of Xf epidemic in Puglia",
  
  titlePanel("Spatial simulation of Xf epidemic in Puglia"),
  
  sidebarLayout(
    sidebarPanel(
      
      actionButton("run_button", "Run simulation"),
      
      textOutput("text"),
      
      h4("Settings for Xylella epidemic"),
      
      radioButtons("initial_location", "Location of first infection:",
                   c("Centre of first known cluster" = "fixed",
                     "Random (high olive density)" = "high",
                     "Random (low olive density)" = "low")),
      
      sliderInput(inputId="number_of_years",
                  label="Number of years to simulate",
                  value=15, min=1, max=40, step=1),
      
      sliderInput(inputId="r",
                  label="Infection growth rate",
                  value=1, min=0, max=2, step=0.05),
      
      sliderInput(inputId="d50_short",
                  label="Median local dispersal distance (km)",
                  value=0.5, min=0.01, max=5, step=0.1),
      
      sliderInput(inputId="d50_long",
                  label="Median long range jump distance (km)",
                  value=15, min=5, max=40, step=1),
      
      sliderInput(inputId="pLong",
                  label="Proportion of long range range jumps",
                  value=0.0025, min=0, max=0.01, step=0.0005),
      
      h4("Settings for management in buffer zone"),
      
      sliderInput(inputId="samplesPerCell",
                  label="Trees tested per grid cell",
                  value=0, min=0, max=100, step=5),
      
      sliderInput(inputId="testFNR",
                  label="Diagnostic testing false negative rate",
                  value=0, min=0, max=1, step=0.01)
      
    ),
    
    # *Output() functions
    mainPanel(tabsetPanel(
      tabPanel("Time trends", plotOutput(outputId = "time_plot", height = "550px", width = "700px")),
      tabPanel("Final infection", plotOutput(outputId = "infection_plot", height = "550px", width = "700px")),
      tabPanel("Year of infection", plotOutput(outputId = "infection_year_plot", height = "550px", width = "700px")),
      tabPanel("Year of removal", plotOutput(outputId = "removal_year_plot", height = "550px", width = "700px")),
      tabPanel("Olive cover", plotOutput(outputId = "olive_plot", height = "550px", width = "700px")),
      tabPanel("Zones", plotOutput(outputId = "zone_plot", height = "550px", width = "700px")),
      tabPanel("Model output", 
               downloadButton("downloadData", "Download"),
               tableOutput(outputId = "data_table") )
    ))
  )
)

server <- function(input, output){
  
  output$olive_plot = renderPlot(expr={ 
    print(ggplot(X, aes(x=x, y=y)) + 
            geom_tile(aes(fill=olive)) + coord_fixed() +
            scale_fill_viridis(option="viridis", direction=1, limits=c(0,1)) +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  legend.position = "right",
                  legend.title = element_blank(),
                  legend.text = element_text(size=14))
    )
  })
  
  output$zone_plot = renderPlot(expr={ 
    print(ggplot(X, aes(x=x, y=y)) + 
            geom_tile(aes(fill=zone)) + coord_fixed() +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  legend.position = "right",
                  legend.title = element_blank(),
                  legend.text = element_text(size=16))
    )
  })
  
  observeEvent(input$run_button, {
    
    simEpidemic <- withCallingHandlers(
      {
      shinyjs::html(id = "text", html = "")
      sim_spread(df=X, r=input$r, number_of_years=input$number_of_years, 
                 d50_short = input$d50_short, d50_long = input$d50_long, 
                 pLong = input$pLong, initial_location = input$initial_location,
                 samplesPerCell=input$samplesPerCell, testFNR=input$testFNR)
      },
      message = function(m) {
        shinyjs::html(id = "text", html = m$message, add = FALSE)
      },
      warning = function(m) {
        shinyjs::html(id = "text", html = m$message, add = FALSE)
      }
    )
    
    output$time_plot = renderPlot(expr={
        nYrs = input$number_of_years
        simEpidemic$year_infected = factor(simEpidemic$year_infected, levels=0:nYrs)
        simEpidemic$year_removed = factor(simEpidemic$year_removed, levels=0:nYrs)
        plotdf = data.frame(year=0:nYrs,
                   total_infection=colSums(simEpidemic[,paste0("Xf_yr", 0:nYrs)]),
                   t(apply(simEpidemic[,paste0("Xf_yr", 0:nYrs)], 2, function(x){
                     tapply(x, simEpidemic$zone, sum)})),
                   total_removed = as.numeric(cumsum(table(simEpidemic$year_removed))),
                   new_removals = as.numeric(table(simEpidemic$year_removed)) )
        plotdf = melt(plotdf, id.vars="year")
        print(ggplot(plotdf, aes(x=year, y=value)) + geom_line(col="tomato", size=1) + 
              facet_wrap(~variable, scales="free_y", ncol=2) +
              xlab("Year") +
              theme(axis.title.y = element_blank(),
                    axis.title.x = element_text(size=16),
                    axis.text = element_text(size=12),
                    strip.text = element_text(size=16)))
      })
    
    output$data_table <- renderTable({ 
      nYrs = input$number_of_years
      simEpidemic$year_infected = factor(simEpidemic$year_infected, levels=0:nYrs)
      simEpidemic$year_removed = factor(simEpidemic$year_removed, levels=0:nYrs)
      df = data.frame(year=0:nYrs,
                          total_infection=colSums(simEpidemic[,paste0("Xf_yr", 0:nYrs)]),
                          t(apply(simEpidemic[,paste0("Xf_yr", 0:nYrs)], 2, function(x){
                            tapply(x, simEpidemic$zone, sum)})),
                          total_removed = as.numeric(cumsum(table(simEpidemic$year_removed))),
                          new_removals = as.numeric(table(simEpidemic$year_removed)) )
      }, 
      digits=3 
    )
    
    output$downloadData <- downloadHandler(
      filename = paste0("Xf_spread-", Sys.time(), ".csv"),
      content = function(con) {
        nYrs = input$number_of_years
        simEpidemic$year_infected = factor(simEpidemic$year_infected, levels=0:nYrs)
        simEpidemic$year_removed = factor(simEpidemic$year_removed, levels=0:nYrs)
        df = data.frame(year=0:nYrs,
                        total_infection=colSums(simEpidemic[,paste0("Xf_yr", 0:nYrs)]),
                        t(apply(simEpidemic[,paste0("Xf_yr", 0:nYrs)], 2, function(x){
                          tapply(x, simEpidemic$zone, sum)})),
                        total_removed = as.numeric(cumsum(table(simEpidemic$year_removed))),
                        new_removals = as.numeric(table(simEpidemic$year_removed)) )
                        
        write.csv(df, con, row.names = FALSE)
      },
      contentType="text/csv"
    )
    
    output$infection_plot = renderPlot(expr={ 
        #simEpidemic$Xf[simEpidemic$Xf==0] = NA
        print(ggplot(simEpidemic, aes(x=x, y=y)) + 
              geom_tile(aes(fill=Xf)) + coord_fixed() +
              #scale_fill_gradientn(colours=c("yellow", "orange", "tomato", "red", "red4"), na.value = "grey90", limit=c(0,1)) +
              #scale_fill_viridis(option="inferno", direction=1, limits=c(0,1)) +
              scale_fill_gradientn(colours=rev(c("yellow", "orange", "tomato", "red", "red4", "black", "white")), na.value = "grey90", limit=c(0,1)) +
                geom_point(data=simEpidemic[which(simEpidemic$year_infected==0),], 
                         aes(x=x,y=y,Xf=NULL), shape=1, size=4, col="#009E73", stroke=1) +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "right",
                    legend.title = element_blank(),
                    legend.text = element_text(size=14),
                    panel.background = element_rect(fill="lightblue"))
        ) 
      })
    
    output$infection_year_plot = renderPlot(expr={ 
        print(ggplot(simEpidemic, aes(x=x, y=y)) + 
              geom_tile(aes(fill=year_infected)) + coord_fixed() +
              #scale_fill_gradientn(colours=c("yellow", "orange", "tomato", "red", "red4"), na.value = "grey90") +
              scale_fill_viridis(option="viridis", direction=1) +
              geom_point(data=simEpidemic[which(simEpidemic$year_infected==0),], 
                         aes(x=x,y=y,Xf=NULL), shape=1, size=4, col="#009E73", stroke=1) +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "right",
                    legend.title = element_blank(),
                    legend.text = element_text(size=14),
                    panel.background = element_rect(fill="lightblue"))
        )
      })
    
    output$removal_year_plot = renderPlot(expr={ 
        print(ggplot(simEpidemic, aes(x=x, y=y)) + 
              geom_tile(aes(fill=year_removed)) + coord_fixed() +
              #scale_fill_gradientn(colours=c("yellow", "orange", "tomato", "red", "red4"), na.value = "grey90") +
              scale_fill_viridis(option="viridis", direction=1) +
              geom_point(data=simEpidemic[which(simEpidemic$year_infected==0),], 
                         aes(x=x,y=y,Xf=NULL), shape=1, size=4, col="#009E73", stroke=1) +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "right",
                    legend.title = element_blank(),
                    legend.text = element_text(size=14),
                    panel.background = element_rect(fill="lightblue"))
        )
      })
      
  })   
}

#shinyApp(ui=ui, server=server)
runApp(shinyApp(ui, server), launch.browser = TRUE)
