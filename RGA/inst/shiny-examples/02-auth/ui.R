library(shiny)

shinyUI(
    fluidPage(
        titlePanel("RGA package: authorization demo"),
        fluidRow(
            uiOutput("auth_links")
        ),
        fluidRow(
            DT::dataTableOutput("accounts")
        )
    )
)
