library(shiny)

shinyUI(
    fluidPage(
        titlePanel("Google Analytics: Dimensions & Metrics"),
        absolutePanel(style = "z-index: 1000", draggable = TRUE, fixed = TRUE, width = 280, height = 300,
                      top = "auto", left = "auto", right = 20, bottom = 20,
                      wellPanel(
                          checkboxGroupInput(inputId = "columns",
                                             label = "Columns to show:",
                                             choices = cn[!cn %in% selected])
                      )
        ),
        fluidRow(
            column(4,
                   checkboxInput(inputId = "segments", label = "Show allowed in segments only"),
                   checkboxInput(inputId = "status", label = "Show deprecated"),
                   checkboxInput(inputId = "calc", label = "Show culculated only")
            ),
            column(4,
                   selectInput(inputId = "group", label = "Group:",
                               choices = c("All", unique(ga$group)))
            ),
            column(4,
                   selectInput(inputId = "type", label = "Type:",
                               choices = c("All", unique(ga$type)))
            )
        ),
        fluidRow(
            DT::dataTableOutput(outputId = "table")
        )
    )
)
