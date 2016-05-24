

library(shiny)
library(rangeMapper)
con = rangeMap.open ( options('rangeMapper.path')$rangeMapper.path )

cols =  c( list(rangeMapper = palette_rangemap('divergent')), rangeMapper:::brewer.pal.get () )

shinyServer(function(input, output) {

  output$plot <- renderPlot({

    if(!is.null(input$maps) )  {
      plot( rangeMap.fetch(con, input$maps, spatial = FALSE), colours = cols[input$colorRamp][[1]]  )
      }

  })

  output$summary <- renderTable({
    rangeMapProjInfo(con)
  })



})

