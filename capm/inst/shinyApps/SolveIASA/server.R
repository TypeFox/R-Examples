library(ggplot2)
shinyServer(function(input, output) {
  
  #### Functions ####
  
  # Parameters
  parms <- function() {
    return(c(
      b1 = input$b1, b2 = input$b2, df1 = input$df1,
      dm1 = input$dm1, df2 = input$df2, dm2 = input$dm2,
      sf1 = input$sf1, sf2 = input$sf2, sm1 = input$sm1,
      sm2 = input$sm2, k1 = input$k1, k2 = input$k2,
      h1 = input$h1, h2 = input$h2, ab = input$ab,
      ad = input$ad, v = input$v, z = input$z))
  } 
  
  # Initial conditions
  inits <- function() {
    return(c(
      f1 = input$f1, fs1 = input$fs1,
      m1 = input$m1, ms1 = input$ms1,
      f2 = input$f2, fs2 = input$fs2,
      m2 = input$m2, ms2 = input$ms2))
  }
  
  # Solution for point estimates
  SIASA <- function() {
    sol <- SolveIASA(pars = parms(),
                     init = inits(),
                     time = seq(0, input$time,
                                by = input$t_steps))
    return(sol)
  }
  
  # Scenarios
  SIASA2 <- function() {
    sol <- SolveIASA(
      pars = parms(),
      init = inits(),
      time = seq(0, input$time,
                 by = input$t_steps),
      s.range = seq(min(input$s.range), max(input$s.range),
                    l = input$s.intr),
      ab.range <- c(min(input$ab.range), max(input$ab.range)),
      ad.range <- c(min(input$ad.range), max(input$ad.range)),
      im.range <- c(input$im.1, input$im.2),
      s.fm <- !input$s.fm)
    return(sol)
  }
  
  outputVar <- reactive({
    switch(input$output_var,
           'Owned intact animals (n1)' = 'n1',
           'Owned sterilized animals (ns1)' = 'ns1',
           'Stray intact animals (n2)' = 'n2',
           'Stray sterilized animals (ns2)' = 'ns2',
           'Owned animals (N1)' = 'N1',
           'Stray animals (N2)' = 'N2',
           'Total population (N)' = 'N',
           'Owned intact females (f1)' = 'f1',
           'Owned sterilized females (fs1)' = 'fs1',
           'Owned intact males (m1)' = 'm1',
           'Owned sterilized males (ms1)' = 'ms1',
           'Stray intact females (f2)' = 'f2',
           'Stray sterilized females (fs2)' = 'fs2',
           'Stray intact males (m2)' = 'm2',
           'Stray sterilized males (ms2)' = 'ms2')
  })
  
  outputVar2 <- reactive({
    switch(input$output_var2,
           'Intact females (f)' = 'f',
           'Sterilized females (fs)' = 'fs',
           'Intact males (m)' = 'm',
           'Sterilized males (ms)' = 'ms',
           'Intact animals (n)' = 'n',
           'Sterilized animals (ns)' = 'ns',
           'Total population (N)' = 'N')
  })
  
  Local <- function() {
    return(CalculateLocalSens(SIASA(), sensv = outputVar()))
  }
  
  GlobalAll <- function() {
    ranges <- SetRanges(parms(), range = input$range)
    glob.all <- CalculateGlobalSens(SIASA(), ranges = ranges,
                                    sensv = outputVar(), all = T)
    return(glob.all)
  }
  
  Global <- function() {
    ranges <- SetRanges(parms(), range = input$range)
    glob <- CalculateGlobalSens(SIASA(), ranges = ranges,
                                sensv = outputVar())
    return(glob)
  }
  
  #### Outputs ####
  output$points_p <- renderPlot({
    plot(PlotModels(SIASA(), variable = outputVar()))
  })
  
  output$points_t <- renderDataTable({
    SIASA()$results
  }, options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 5))
  
  output$scenarios <- renderPlot({
    if (input$create.scenarios == 0) {
      return()
    }
   PlotModels(SIASA2(), variable = outputVar2())
  })
  
  output$local <- renderPlot({
    if (input$sensitivities == 0) {
      return()
    }
    
    isolate({
      PlotLocalSens(Local(), ax.size = 15) +
        theme(text=element_text(size=15),
              legend.title = element_blank())
    })
  })
  
  output$globalall <- renderPlot({
    if (input$sensitivities == 0) {
      return()
    }
    
    isolate({
      pl <- PlotGlobalSens(GlobalAll()) +
        theme(text = element_text(size = 15),
              legend.title = element_blank())
      plot(pl)
    })
  })
  
  output$global <- renderPlot({
    if (input$sensitivities == 0) {
      return()
    }
    
    isolate({
      pl <- PlotGlobalSens(Global()) +
        theme(text = element_text(size = 15),
              legend.title = element_blank())
      plot(pl)
    })
  })
})