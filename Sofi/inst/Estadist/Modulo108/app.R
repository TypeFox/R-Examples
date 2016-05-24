library(shiny)

source("linked_scatter.R")

ui <- fixedPage(
  h2("Module example"),
  linkedScatterUI("scatters"),
  textOutput("summary"),
  actionButton("sal", "Salir")
)

server <- function(input, output, session) {
  df <- callModule(linkedScatter, "scatters", reactive(mpg),
    left = reactive(c("cty", "hwy")),
    right = reactive(c("drv", "hwy"))
  )

  output$summary <- renderText({
    sprintf("%d observation(s) selected", nrow(dplyr::filter(df(), selected_)))
  })
  
  observeEvent(input$sal, {
    stopApp()
  })
}

shinyApp(ui, server)
