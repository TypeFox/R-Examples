library(shiny)

# Define UI for FindRegLine application
navbarPage(
  title = "Can You Find the Regression Line?",
  tabPanel(
    title = "The Game",
    sidebarPanel(
      conditionalPanel(
        condition = "input.submit == 0 || output.beginning == true",
        textInput("player","Enter your name"),
        helpText("How good are you at finding where the regression line lies?",
             "Adjust the sliders below to change the y-intercept and the",
             "slope of the solid line.")
        ),
      conditionalPanel(
        condition = "input.submit == 0 || output.beginning == true || output.playing == true",
        uiOutput("aslider"),
        uiOutput("bslider"),
        br(),
        helpText("When you are satisfied with the line, press the button below."),
        actionButton("submit","Submit New Guess"),
        br(),
        helpText("Each press of the button counts as one turn.  Your score is the sum",
             "of the number of turns you have taken so far and a 'closeness' measurement.",
             "The closeness measurement is always 100 at the start of the game, and can ",
             "get as low as zero if your line matches the regression line exactly.",
              "Therefore, lower scores are better!")
        ),
      conditionalPanel(
        condition = "output.playing == true",
        actionButton("enditall",label="I quit:  show the regression line!")
      ),
      conditionalPanel(
        condition = "output.reporting == true",
        actionButton("reset", "Play Again")
      )
    ),
  mainPanel(
      conditionalPanel(
        condition = "input.player != '' && output.reporting == true",
        htmlOutput("rank")
    ),
      conditionalPanel(
        condition="input.submit == 0 || (output.beginning == true || output.playing == true)",
        plotOutput("gamecloud"),
        HTML("<h5>(Note: x and y axes for the above plot are in yellow.)</h5>"),
        br()
        ),
      conditionalPanel(
        condition="output.reporting == true",
        plotOutput("finalcloud")
        ),
      fluidRow(
        column(width = 6,
               tableOutput("score")),
        column(width = 3,
               conditionalPanel(
                 condition="output.reporting == true",
                 tableOutput("revelation")
                 )
        )
      )
    )
  ),
  tabPanel(
    title = "Leader Board",
    fluidRow(
      column(width = 3,
             actionButton("updateBoard", "Update (as others play)")),
      column(width = 3,
             downloadButton('downloadData', 'Download'))
    ),
    HTML("<p> </p>"),
    DT::dataTableOutput("leaders")
  ),
  tabPanel(
    title = "About",
    includeHTML("infoFiles/about.html")
  )
)
