shinyServer(function(input, output) {

  vals <- reactive(function() {
    n <- as.numeric(input$n)
    x <- switch(input$family,
                "Normal" = rnorm(n),
                "Exponential" = rexp(n),
                "Symmetric, long-tailed" = rt(n, df=3),
                "Skew, long-tailed" = rlnorm(n)
                )
    if(input$algorithm == "DIY")
      bw <- as.numeric(input$bw)
    else
      bw <- input$algorithm
    list(n=n, x=x, kernel=input$kernel, bw=bw)
  })
  
  output$main_plot <- reactivePlot(function() {
    l <- vals()
    with(l, {
         plot(density(x, kernel=kernel, bw=bw), main=sprintf("Kernel: %s, bw: %s", kernel, bw))
         points(x, abs(jitter(rep(0, length(x)))), pch=16, col="gray80")
         
         })
  })


  output$summary <- reactiveTable(function() {
    l <- vals()
    
    x <- c("algorithm"=input$algorithm,
           "n" = l$n,
           "kernel" = l$kernel
           )

    data.frame(Value=x, stringsAsFactors=FALSE)
        
  })
})
