shinyServer(function(input, output) {
  library(rivr)
  So = reactive({as.numeric(input$So)})
  n = reactive({input$n})
  Q = reactive({input$Q})
  Cm  = reactive({as.numeric(input$Cm)})
  B = reactive({input$B})
  SS = reactive({input$SS})
  stepdist = reactive({as.integer(input$stepdist)})
  totaldist = reactive({input$totaldist})
  y0 = reactive({input$y0}) 
  g = reactive({ifelse(Cm() > 1.0, 32.2, 9.81)})
  z0 = reactive({input$z0})
  x0 = reactive({input$x0})
  res = reactive({compute_profile(So(), n(), Q(), y0(), Cm(), g(), B(), SS(), 
    z0(), x0(), stepdist(), totaldist())})
  output$main_plot <- renderPlot({plot(res())})
  output$main_table <- renderTable({
    as.data.frame(res())
  }, include.rownames = FALSE)
  output$info_table = renderTable({  
    t(setNames(
      data.frame(
        rivr:::get_profile(So(), n(), Q(), g(), Cm(), B(), SS(), y0()),
        round(normal_depth(So(), n(), Q(), y0(), Cm(), B(), SS()), 2),
        round(critical_depth(Q(), y0(), g(), B(), SS()), 2)
      ),
      c("Profile type",
        "Normal depth", 
        "Critical depth"
      )
    ))
  }, include.colnames = FALSE)
  output$example_table = renderTable({setNames(
    data.frame(
      c("M1", "M2", "S2", "S3"),
      c("2.7", "0.64", "2.65", "0.5"),
      c("250", "250", "250", "250"),
      c("0.001", "0.001", "0.005", "0.005"),
      c("100", "100", "10", "10"),
      c("0.045", "0.045", "0.01", "0.01"),
      rep("US Customary", 4)
    ),
    c("Profile type", "Control section depth", "Channel flow", 
    "Channel slope", "channel width", "Manning's roughness", "Units"))
  }, include.rownames = FALSE)
})
