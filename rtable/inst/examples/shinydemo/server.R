library(shiny)

shinyServer(function(input, output) {
  
  ft = reactive({
    data = lsmtab
    data$"least-squares means" = sprintf("%.3f (%.2f)", data$lsmean, data$SE)
    data$"cont. int." = sprintf("%.2f - %.2f", data$lower.CL, data$upper.CL)
    data$"df" = sprintf("%.1f", data$df )
    data$"colorindex" = cut( lsmtab$lsmean, 
                             breaks = quantile( lsmtab$lsmean, probs = c(0:4)/4 ), 
                             include.lowest = T, labels = F )    
    data$"color" = c(input$q1, input$q2, input$q3, input$q4 )[data$"colorindex"]
    data
  })
  
  output$flextable <- renderFlexTable({
    
    data = ft()
    
    options( "ReporteRs-fontsize" = input$fontsize )
    ft = FlexPivot(data, 
      id = "nitro", transpose = "Variety", space.table = input$tablestyle, 
      columns = c("least-squares means", "df", "cont. int."), 
      color = c("least-squares means"="color") )
    
    ft[,] = parCenter(padding.left=input$margin, padding.right = input$margin)
    ft
  })

  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste("FlexTable", input$filetype, sep = "")
    },
    
    content = function(file) {
      data = ft()
      
      options( "ReporteRs-fontsize" = input$fontsize )
      ft = FlexPivot(data, 
                     id = "nitro", transpose = "Variety", space.table = input$tablestyle, 
                     columns = c("least-squares means", "df", "cont. int."), 
                     color = c("least-squares means"="color") )
      
      ft[,] = parCenter(padding.left=input$margin, padding.right = input$margin)
      
      if( input$filetype == ".docx" ){
        doc = docx()
        doc = addSection(doc, landscape = TRUE)
      } else {
        doc = pptx()
        doc = addSlide(doc, "Title and Content")
      }
      doc = addFlexTable(doc, ft, par.properties = parCenter() )
      if( input$filetype == ".docx" ){
        doc = addSection(doc, landscape = FALSE)
      }
      writeDoc( doc, file = file )
    }
  )
})
