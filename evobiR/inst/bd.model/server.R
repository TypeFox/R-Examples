library(shiny)
library(phytools)
shinyServer(function(input, output) {
  tree <- reactive({
      set.seed(input$seed.val)
      pbtree(b = input$birth, 
             d = input$death,
             t = input$time, 
             scale = 1,
             nsim = 4,
             extant.only = input$extinct) 
  })
  output$treePlot <- renderPlot({
    par(mfcol=c(2,2))
    for(i in 1:4){
      plot.phylo(tree()[[i]], show.tip.label=F)
      mtext(paste("N =", length(tree()[[i]]$tip.label)), side = 1, line = 0)
    }
  })  
})


