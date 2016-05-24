shinyServer(function(input, output) {

  vals <- reactive({
    n <- as.numeric(input$n)
    x <- rnorm(input$n)
    
    if(input$algorithm == "DIY")
      breaks <- as.numeric(input$n_bins)
    else
      breaks <- input$algorithm

    list(n=n, x=x, breaks=breaks)

  })
 
  
  output$main_plot <- renderPlot({
    l <- vals()
    with(vals, 
         hist(x,
              breaks = breaks,
              xlab = "normal data",
              main = sprintf("n: %s, (n^(1/3): %03f", n, n^(1/3)))
         )
  })

  ## "DRY" this up
  output$summary <- renderTable({
    l <- vals()

    with(vals, {

      out <- hist(x, breaks = breaks, plot=FALSE)
    
      x <- c("Algorithm"=input$algorithm,
             "n" = n,
             "No bins" = length(out$breaks) - 1
             )
      
      data.frame(Value=x, stringsAsFactors=FALSE)
    })
    
  })
})
