shinyUI(pageWithSidebar(
  headerPanel('Item model explorer'),
  sidebarPanel(
    selectInput('model', 'Item Model:', c("dichotomous", "graded", "nominal")),
    conditionalPanel(
      condition = "input.model != 'dichotomous'",
      sliderInput("outcomes", "Outcomes:", min=3, max=20, value=5)),
    conditionalPanel(
      condition = "input.model == 'nominal'",
      selectInput("nominalTa", "T.a matrix:", c("trend", "id")),
      selectInput("nominalTc", "T.c matrix:", c("trend", "id", "partial credit"))),
    actionButton("drawNewParametersAction", label = "Draw new parameters"),
    checkboxInput("showParameters", label = "Show/Edit parameters", value = TRUE),
    conditionalPanel(
      condition = "input.showParameters",
      hr(),
      selectInput('editPar', 'Edit Parameter:', 'a'),
      sliderInput('editParValue', "Parameter Value:", min=-5, max=5, value=0, step=.01, ticks=FALSE),
      fluidRow(actionButton("setParValue0", label = "Set 0"),
               actionButton("setParValue1", label = "Set 1"),
               actionButton("setAllValue0", label = "Set all 0")),
      hr(),
      tableOutput("parView")
    )
  ),
  mainPanel(
    plotOutput('plot1'),
    plotOutput('info')
  )
))
