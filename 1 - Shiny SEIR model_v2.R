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
                      latent_period_yrs = 1,
                      yrs_to_death_E = 200,
                      yrs_to_death_I = 20,
                      number_of_years = 20,
                      plot=TRUE){
  # Solve continuous SEIR with latent infections
  # Adapted from Keeling and Rohani, Modeling Infectious Diseases in Humans and Animals: 
  # http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/index.html  
  
  # SEIR function 
  seir <- function(time, state, parameters) {
    require(deSolve)
    with(as.list(c(state, parameters)), {
      dS <- -beta*S*I - betaL*S*E
      dE <- beta*S*I + betaL*S*E - sigma*E
      dI <- sigma*E - gamma*I
      dR <-  gamma*I
      return(list(c(dS, dE, dI, dR)))
    })
  }
  
  ## Initial proportion in each compartment
  init <- c(S = 1-initial_infection, E = initial_infection, I=0.0, R = 0.0)
  ## parameter list (beta: infection parameter; gamma: recovery parameter)
  parameters <- c(beta = transmission_rate_I, betaL=transmission_rate_E, 
                  gamma=1/yrs_to_death_I, gammaL=1/yrs_to_death_E,
                  sigma=1/latent_period_yrs)
  ## Time frame
  times <- seq(0, number_of_years, by = 0.1)
  
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  out <- as.data.frame(ode(y = init, times = times, func = seir, parms = parameters))
  out$diseased = 1-out$S
  
  if(plot){
    df = melt(out, id.vars="time")
    print(ggplot(df, aes(x=time, y=value, colour=variable)) +
            geom_line(lwd=1) +
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
    
      h4("Settings for asymptomatic infecteds"),
               
      sliderInput(inputId="transmission_rate_E",
                           label="Transmission rate",
                           value=0.5, min=0, max=6.0, step=0.1),
               
      sliderInput(inputId="yrs_to_death_E",
                           label="Average time to death (years)",
                           value=200, min=0.1, max=200.0, step=0.1),
               
      sliderInput(inputId="latent_period_yrs",
                           label="Average length of asymptomatic period (years)",
                           value=1, min=0.1, max=3.0, step=0.1),
    
      h4("Settings for symptomatic infecteds"),
    
      sliderInput(inputId="transmission_rate_I",
              label="Transmission rate",
              value=2, min=0, max=6.0, step=0.1),
    
      sliderInput(inputId="yrs_to_death_I",
                label="Average time to death (years)",
                value=5, min=0.1, max=10.0, step=0.1),
    
      h4("General settings"),
    
      sliderInput(inputId="number_of_years",
                label="Number of years",
                value=25, min=1, max=100, step=1)
    ),
  
    # *Output() functions
    mainPanel(tabsetPanel(
      tabPanel("Plot", plotOutput(outputId = "infection_plot", height = "700px", width = "900px")),
      tabPanel("Data", 
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
                yrs_to_death_E = input$yrs_to_death_E,
                yrs_to_death_I = input$yrs_to_death_I,
                number_of_years = input$number_of_years,
                plot=FALSE)
  })
  
  output$infection_plot <- renderPlot(expr={ 
    df = melt(sim(), id.vars="time")
    print(ggplot(df, aes(x=time, y=value, colour=variable)) +
            geom_line(lwd=1) +
            ylab("Proportion of population") +
            xlab("Years since infection") +
            theme(legend.position = "top",
                  legend.title = element_blank(),
                  legend.text = element_text(size=20),
                  axis.title = element_text(size=20),
                  axis.text = element_text(size=12)))
    }) 
  
  output$data_table <- renderTable({ 
    df = sim()
    names(df)[1] = "Year"
    df = df[df$Year %% 1 == 0,]}, digits=3 
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