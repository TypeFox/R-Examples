library(shiny)

# Define UI for CoinFlip application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Simulating With Two Possible Outcomes"),
  
  # Sidebar
  sidebarPanel(
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      helpText("You can use your own names for each of the two outcomes, or just",
             "stick with the default names below"),
      textInput("success","Name of a Success","Success"),
      textInput("failure","Name of a Failure","Failure"),
      helpText("Enter the chance of success on each trial, as a decimal number",
             "between 0 and 1, not as a percentage."),
      numericInput("p","Chance of Success",0.50,min=0,max=1),
      numericInput("n","Number of Trials",16,min=1,step=1),
      helpText("Enter below the number of successes you actually observed in your",
             "trials"),
      numericInput("xObs","Observed Successes",value=15,min=0,step=1)
      
    ),
    helpText("One simulation means the machine will go though all the trials",
             "and count up the number of successes.  How many simulations do",
             "you want the machine to perform at once?  (Limit is 10000.)"),
    numericInput("sims","Number of Simulations at Once",1,min=0,step=1),
    actionButton("resample","Simulate Now"),
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      actionButton("reset","Start Over")
    )
    
    ),

  
  # Here comes the main panel
  
     mainPanel(
    
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      plotOutput("barGraphInitial")
#      p(textOutput("remarksInitial")) weird bug here I think
      ),
    
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      tabsetPanel(selected="Latest Simulation",
        tabPanel("Latest Simulation",
               conditionalPanel(
                  condition="input.n <= 20",
                  verbatimTextOutput("fullSim")
                  ),
               plotOutput("barGraphLatest"),
               p(textOutput("remarksLatest1")),
               tableOutput("summary1")),
        tabPanel("Histogram of Simulations",
               plotOutput("histogram"),
               p(textOutput("remarksLatest2")),
               tableOutput("summary2")),
        tabPanel("Probability Distribution",
                 plotOutput("probHistogram"),
                 p(textOutput("remarkProb"))
                 ),
        id="MyPanel"
    )
    )
    
    
  )
  
))
