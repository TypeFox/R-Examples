library(shiny)
library(chipPCR)

# server for the Shiny app
shinyServer(function(input, output) {
  res.AmpSim <- reactive({AmpSim(cyc = 1:input$cycles, 
                                 b.eff = input$b.eff, 
                                 bl = input$bl,
                                 ampl = input$ampl, 
                                 Cq = input$Cq,
                                 noise = input$noise, 
                                 nnl = input$nnl, 
                                 nnl.method = input$nnl.method)
  })
  # Use bg.max to calculate the SDM and alike
  res.bg <- reactive({bg.max(res.AmpSim())
  })
  
  res.th.cyc <- reactive({th.cyc(res.AmpSim()[, 1], res.AmpSim()[, 2], 
                                 r = input$th.r, auto = input$th.auto, linear = input$th.lin)})
  
  
  # Create a plot
  output$AmpSimPlot <- renderPlot({
    plot(res.AmpSim(), main = "Simulated curve", type = "b")
    lines(c(res.AmpSim()[1, 1], res.th.cyc()[1]), c(res.th.cyc()[2], res.th.cyc()[2]), 
          lty = "dashed", col = "orange")
    lines(c(res.th.cyc()[1], res.th.cyc()[1]), 
          c(min(res.AmpSim()[, 2]), res.th.cyc()[2]), lty = "dashed", col = "orange")
    points(res.th.cyc(), pch = 20, col = "orange")
  })
  
  output$inderPlot <- renderPlot({
    plot(res.bg(), main = "Calculation of curve parameters")
  })
  
  output$bgTable <- renderTable({
    res.bg()
  })
  
  output$bgSummary <- renderPrint({
    summary(res.bg())
    cat("\nThreshold cycle: ", format(res.th.cyc()[1], digits = 4))
  })
  
})

