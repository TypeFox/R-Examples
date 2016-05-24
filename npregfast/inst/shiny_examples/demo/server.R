detach("package:npregfast")
library(shiny)
#library(shinyjs)
library(miniUI)
library(wesanderson)
library(npregfast)


shinyServer(function(input, output) {
  
  data  <- barnacle
  # For storing which rows have been excluded
  vals <- reactiveValues(
    keeprows = rep(TRUE, nrow(data))
  )
  
  
  output$distPlot <- renderPlot({
    
    # Plot the kept and excluded points as two separate data sets
    keep    <- data[ vals$keeprows, , drop = FALSE]
    exclude <- data[!vals$keeprows, , drop = FALSE]
    
    
    
    if(input$type == "with"){
      form <- DW ~ RC : F
    }else{
      form <- DW ~ RC
    }
    
    if(input$selband == "cv"){
      h0 = -1
      h = -1
    }else{
      h0 = input$band
      h = input$band
    }
    
    if(input$poly == 1) {der <- as.numeric(input$der1)}
    if(input$poly == 2) {der <- as.numeric(input$der2)}
    if(input$poly == 3) {der <- as.numeric(input$der3)}
    
    
    
  
    fit <- frfast(form, data = keep, nboot = 100, kernel = input$kernel,
                  h0 = h0, h = h, p = input$poly)
    plot(fit, der = der, points = input$show_points, 
         CIcol = input$colci, col = input$colmu, CIlwd = 2, 
         ablinecol = "#24281A", pch = 16, pcol = input$pcol)
    
  })
  
  # Toggle points that are clicked
  observeEvent(input$plot1_click, {
    res <- nearPoints(data, input$plot1_click, allRows = TRUE,
                      xvar = "RC", yvar = "DW")
    
    vals$keeprows <- xor(vals$keeprows, res$selected_)
  })
  
  # Toggle points that are brushed, when button is clicked
  observeEvent(input$exclude_toggle, {
    res <- brushedPoints(data, input$plot1_brush, allRows = TRUE,
                         xvar = "RC", yvar = "DW")
    
    vals$keeprows <- xor(vals$keeprows, res$selected_)
  })
  
  # Reset all points
  observeEvent(input$exclude_reset, {
    vals$keeprows <- rep(TRUE, nrow(data))
  })
  
  observeEvent(input$info_btn, {
    shinyjs::info("This plot supports mouse based-interaction, via clicking and brushing. The points selected or included in the selected area will be deleted and will not be considered in the analysis. In order to use this option correctly, the selection of the points must be carried out with only one graphical output marked and without interaction. Once the points have been deleted, the other graphical output and estimation options can be marked.")
      })
  
  
  
  # hide the loading message
  shinyjs::hide("loading-content", TRUE, "fade")
})


