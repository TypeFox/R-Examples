library(shiny)

shinyUI(
    fluidPage(
        titlePanel("RGA package: realtime demo"),
        fluidRow(
            uiOutput("auth_button")
        ),
        fluidRow(
            textOutput("active_users")
        )
    )
)
