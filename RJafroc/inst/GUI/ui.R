shinyUI(fluidPage(
  titlePanel(sprintf("Welcome to RJafroc v%s!", packageVersion("RJafroc"))),
    mainPanel(
      fileInput("dataFile", "Select Data File:"),
      br(),
      uiOutput("ui")
    )
))