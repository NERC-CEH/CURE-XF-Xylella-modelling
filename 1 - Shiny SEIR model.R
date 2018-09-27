rm(list=ls())

# install packages if needed
if(!"shiny" %in% installed.packages()){install.packages("shiny", dep=TRUE)}
if(!"deSolve" %in% installed.packages()){install.packages("deSolve", dep=TRUE)}
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2", dep=TRUE)}
if(!"reshape" %in% installed.packages()){install.packages("reshape", dep=TRUE)}

# load libraries
library(shiny)
library(deSolve)
library(ggplot2)
library(reshape)


#####################################################
### continuous SEIR with latent infections
SEIR_model = function(initial_infection=0.001, 
                      transmission_rate_E = 1,
                      transmission_rate_I = 2,
                      latent_period_yrs = 0,
                      yrs_to_death_E = Inf,
                      yrs_to_death_I = 5,
                      vaccination_rate = 0,
                      number_of_years = 20,
                      plot=TRUE){
  # Solve continuous SEIR with latent infections
  # Adapted from Keeling and Rohani, Modeling Infectious Diseases in Humans and Animals: 
  # http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/index.html  
  
  # SEIR function 
  seir <- function(time, state, parameters) {
    require(deSolve)
    with(as.list(c(state, parameters)), {
      dS <- -beta*S*I - betaL*S*E - S*v
      dE <- beta*S*I + betaL*S*E - sigma*E - gammaL*E
      dI <- sigma*E - gamma*I
      dR <-  gamma*I + gammaL*E
      dV <- S*v
      return(list(c(dS, dE, dI, dR, dV)))
    })
  }
  
  ## Initial proportion in each compartment
  init <- c(S = 1-initial_infection, E = initial_infection, I=0.0, R = 0.0, V=0.0)
  ## parameter list (beta: infection parameter; gamma: recovery parameter)
  parameters <- c(beta = transmission_rate_I, betaL=transmission_rate_E, 
                  gamma=1/yrs_to_death_I, gammaL=1/yrs_to_death_E,
                  sigma=1/latent_period_yrs,
                  v = vaccination_rate)
  ## Time frame
  times <- seq(0, number_of_years, by = 0.1)
  
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  out <- as.data.frame(ode(y = init, times = times, func = seir, parms = parameters))
  out$uninfected = out$S + out$V
  out$diseased = 1 - out$S - out$V
  
  if(plot){
    df = melt(out, id.vars="time")
    print(ggplot(df, aes(x=time, y=value, colour=variable, shape=variable)) +
            geom_line(lwd=1) +
            geom_point(data=df[df$time%%2==0,], inherit.aes=TRUE, size=3) +
            ylab("Proportion of population") +
            xlab("Years since infection") +
            theme(legend.title = element_blank()))
  }
  
  return(out)
}


ui <- fluidPage(
  # *Input() functions
  
  title = "Non-spatial SEIR model for disease in a population",
  titlePanel("Non-spatial SEIR model for disease in a population"),
  
  sidebarLayout(
    sidebarPanel(
      
      h3("Settings for asymptomatic infecteds"),
      
      sliderInput(inputId="transmission_rate_E",
                  label="Transmission rate (beta_E)",
                  value=0.5, min=0, max=10.0, step=0.05),
      
      sliderInput(inputId="mortality_E",
                  label="Mortality per year (v_E)",
                  value=0, min=0, max=1, step=0.01),
      
      sliderInput(inputId="latent_period_yrs",
                  label="Average years asymptomatic (1/sigma)",
                  value=1, min=0.1, max=3.0, step=0.1),
      
      h3("Settings for symptomatic infecteds"),
      
      sliderInput(inputId="transmission_rate_I",
                  label="Transmission rate (beta_I)",
                  value=3, min=0, max=10, step=0.1),
      
      sliderInput(inputId="mortality_I",
                  label="Mortality per year (v_I)",
                  value=0.2, min=0, max=1, step=0.01),
      
      h3("Other settings"),
      
      sliderInput(inputId="number_of_years",
                  label="Number of years",
                  value=20, min=1, max=100, step=1),
      
      sliderInput(inputId="vaccination_rate",
                  label="Vaccination rate",
                  value=0, min=0, max=1, step=0.01)
      
      
    ),
    
    # *Output() functions
    mainPanel(tabsetPanel(
      tabPanel("Simple plot", plotOutput(outputId = "simple_plot")),
      tabPanel("Detailed plot", plotOutput(outputId = "infection_plot")),
      tabPanel("Model output", 
               downloadButton("downloadData", "Download"),
               tableOutput(outputId = "data_table") )
    ))
  )
)

server <- function(input, output){
  
  sim <- reactive({
    SEIR_model( transmission_rate_E = input$transmission_rate_E,
                transmission_rate_I = input$transmission_rate_I,
                latent_period_yrs = input$latent_period_yrs,
                yrs_to_death_E = 1/input$mortality_E,
                yrs_to_death_I = 1/input$mortality_I,
                vaccination_rate = input$vaccination_rate,
                number_of_years = input$number_of_years,
                plot=FALSE)
  })
  
  output$infection_plot <- renderPlot(expr={ 
    out = sim()
    out = out[,!(names(out) %in% c("uninfected", "diseased"))]
    df = melt(out, id.vars="time")
    df$label = factor(c("Susceptible", "Infected (asymptomatic)", "Infected (symptomatic)", 
                        "Dead", "Vaccinated")[as.numeric(df$variable)],
                      levels=c("Susceptible", "Infected (asymptomatic)", "Infected (symptomatic)", 
                               "Dead", "Vaccinated"))
    myPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                  "#D55E00", "#CC79A7","#999999")[c(7,3,6,2,4)]
    print(ggplot(df, aes(x=time, y=value, colour=label, shape=label)) +
            geom_line(lwd=1) +
            geom_point(data=df[df$time%%2==0,], inherit.aes=TRUE, size=5) +
            scale_colour_manual(values=myPalette) +
            ylab("Proportion of olive trees") +
            xlab("Years since infection") +
            theme(legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.text = element_text(size=20),
                  axis.title = element_text(size=20),
                  axis.text = element_text(size=12),
                  legend.key.size = unit(1, "cm") )
          )
  }) 
  
  output$simple_plot <- renderPlot(expr={ 
    out = sim()
    out = out[,names(out) %in% c("time", "uninfected", "diseased")]
    df = melt(out, id.vars="time")
    df$label = factor(c("Uninfected", "Diseased")[as.numeric(df$variable)],
                      levels=c("Uninfected", "Diseased"))
    myPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                  "#D55E00", "#CC79A7","#999999")[c(3,7)]
    print(ggplot(df, aes(x=time, y=value, colour=label, shape=label)) +
            geom_line(lwd=1) +
            scale_colour_manual(values=myPalette) +
            ylab("Proportion of olive trees") +
            xlab("Years since infection") +
            theme(legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.text = element_text(size=20),
                  axis.title = element_text(size=20),
                  axis.text = element_text(size=12),
                  legend.key.size = unit(1, "cm") )
    )
  }) 
  
  output$data_table <- renderTable({ 
    df = sim()
    names(df)[1] = "Year"
    df = df[df$Year %% 1 == 0,]}, 
    digits=3 
  )
  
  output$downloadData <- downloadHandler(
    filename = paste0("SEIR_dynamics-", Sys.time(), ".csv"),
    content = function(con) {
      df = sim()
      names(df)[1] = "Year"
      df = df[df$Year %% 1 == 0,]
      write.csv(df, con, row.names = FALSE)
    },
    contentType="text/csv"
  )
}

#shinyApp(ui=ui, server=server)
runApp(shinyApp(ui, server), launch.browser = TRUE)
