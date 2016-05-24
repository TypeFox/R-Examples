shinyServer(function(input, output) {

  output$main_plot <- renderPlot(function() {
    set.seed(100)

    n <- as.numeric(input$n_points)
    pts <- rnorm(n-1, mean=10, sd=4)
    
    lst <- as.numeric(input$nth)

    d <- data.frame(points=c(pts, lst), color=c(rep("blue", n-1),  "red"))


    plot.new()
    plot.window(xlim=c(0,50), ylim=c(0,1))
    axis(1)
    points(pts, .5 + 0*pts, col="blue")
    points(lst, .5, col="red", cex=2)

    abline(v=median(d$points), col="gray50")
    abline(v=mean(d$points), col="black", lty=2)
    abline(v=mean(d$points, trim=.1), col="black", lty=3)
    
    
  })

})
