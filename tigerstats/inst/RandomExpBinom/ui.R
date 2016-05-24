library(shiny)

# Define UI for RandomExpBinom application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("A Randomized Experiment With a Two-Value Response"),
  
  # Sidebar
  sidebarPanel(
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      helpText("You can use your own names for each of the two treatment, or",
             "stick with the default names below"),
      textInput("groupNames","Enter group names (separated by a comma)",
              "GreeterYawns,Control"),
      textInput("groupSizes","Enter group sizes (separated by a comma)",
              "34,16"),
      helpText("You can use your own names for the values of the response,",
             "variable or stick with the default below"),
      textInput("success","Name of a Success","yawns"),
      textInput("failure","Name of a Failure","none"),
      helpText("Enter the number of successes in each group",
             "separated by a comma."),
      textInput("successCounts","Enter success counts",
              "10,4")
    ),
    helpText("One simulation means the machine will randomly assign subjects",
             "to the two groups, with sizes as specified.  How many do",
             "you want to perform at once?  (Limit is 10000.)"),
    numericInput("sims","Number of Simulations at Once",1,min=0,step=1),
    br(),
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
      HTML("<p>Here is a two-way table of the observed results, along with the proportion
           of successes in each group.</p>"),
      tableOutput("initialTwoWay"),
      p(textOutput("remarksInitialMore")),
      plotOutput("barGraphInitial"),
      p(textOutput("remarksInitial"))
      ),
    
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      tabsetPanel(selected="Latest Simulation",
        tabPanel("Latest Simulation",
               HTML("<p>Here is a two-way table of the latest simulation, along with 
                  the proportion of successes in each group.</p>"),
               tableOutput("latestTwoWayBar"),
               p(textOutput("remarksLatest1")),
               plotOutput("barGraphLatest"),
               tableOutput("summary1")),
        tabPanel("Density Curve of Differences",
               HTML("<p>Here is a two-way table of the latest simulation, along with 
                  the proportion of successes in each group.</p>"),
               tableOutput("latestTwoWayDen"),
               p(textOutput("remarksLatest2")),
               plotOutput("density"),
               HTML("
                  <ul>
                    <li>The vertical line (if present) shows the difference in the actual experiment.</li>
                    <li>The red dot shows the difference from the last simulation..</li>
                  </ul>"
                    ),
               tableOutput("summary2")),
        tabPanel("Probability Distribution",
                 plotOutput("normalCurve"),
                 p(textOutput("remarksProb")),
                 HTML("<p><strong>Warning!</strong>  The approximation can be rather poor when
                      there is a small number of successes or a small number of failures in 
                      each group.</p>")
                 ),
        id="MyPanel"
    )
    )
    
    
  )
  
))
