library(shiny)
library(shinyjs)
# Define UI for SlowGoodness application
shinyUI(fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("Chi-Square Goodness-of-Fit Resampling"),
  # Sidebar
  sidebarPanel(
#   conditionalPanel(
#      condition="input.resample == 0 || output.totalPrev == output.total",
    inputPanel(id="setup",
      helpText("Enter the probabilities as decimal numbers.",
                  "If they do not sum to 1, then the",
                  "application will re-scale them for you."),
      textInput("nulls","Null Probabilities (separated by commas)",
                ".17,.17,.17,.17,.17,.17"),
      textInput("obs","Enter Observed Counts (separated by commas)",
                "8,18,11,7,9,7"),
      textInput("names","Enter Level Names (separated by commas)",
                "One,Two,Three,Four,Five,Six")
    ),
    helpText("One simulation means the machine will produce one table of",
             "counts, using the Null probabilities. How many simulations do",
             "you want the machine to perform at once? (Limit is 10000.)"),
    numericInput("sims","Number of Simulations at Once",1,min=0,step=1),
    actionButton("resample","Simulate Now"),
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      actionButton("reset","Start Over")
    )
  ),
  # Here comes the main panel
  # Here comes the main panel
  mainPanel(
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      plotOutput("barGraphInitial"),
      p(textOutput("remarksInitial")),
      tableOutput("obsTable")
    ),
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      tabsetPanel(selected="Latest Simulation",
                  tabPanel("Latest Simulation",
                           plotOutput("barGraphLatest"),
                           p(textOutput("remarksLatest1")),
                           tableOutput("summary1"),
                           p(textOutput("remarksProbBar"))),
                  tabPanel("Density Plot of Simulations",
                           plotOutput("densityplot"),
                           p(textOutput("remarksLatest2")),
                           tableOutput("summary2"),
                           p(textOutput("remarksProbDensity"))),
                  tabPanel("Probability Distribution",
                           plotOutput("chisqCurve"),
                           br(),
                           checkboxInput("compareDen",
                            HTML("Compare with simulated <br>chi-square distribution")),
                           p(textOutput("remarksProb"))
                  ),
                  id="myPanel"
      )
    )
  )
))