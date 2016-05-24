library(shiny)
library(dpcR)

shinyServer(function(input, output) {
  full_model <- reactive({
    adpcr1 <- sim_adpcr(m = input[["m_1"]], n = 765, times = 1000, 
                        pos_sums = FALSE, n_panels = input[["nr_1"]])
    adpcr2 <- sim_adpcr(m = input[["m_2"]], n = 765, times = 1000, 
                        pos_sums = FALSE, n_panels = input[["nr_2"]])
    adpcr2 <- rename_dpcr(adpcr2, exper = "Experiment2")
    adpcr3 <- sim_adpcr(m = input[["m_3"]], n = 765, times = 1000, 
                        pos_sums = FALSE, n_panels = input[["nr_3"]])
    adpcr3 <- rename_dpcr(adpcr3, exper = "Experiment3")
    dat <- bind_dpcr(adpcr1, adpcr2, adpcr3)
    test_res <- test_counts(dat)
    
    list(dat = dat, test_res = test_res)
  })
  
  
  output[["input_data"]] <- renderTable({
    whole_data <- full_model()[["dat"]]
    storage.mode(whole_data) <- "integer"
    slot(whole_data, ".Data")
  })
  
  output[["input_data_summ"]] <- renderPrint({
    summary(full_model()[["dat"]])
  })
  
  output[["test_plot_agg"]] <- renderPlot({
    plot(full_model()[["test_res"]], aggregate = TRUE)
  })
  
  output[["test_plot"]] <- renderPlot({
    other_sums <- summary(full_model()[["dat"]], print = FALSE)[["summary"]]
    plot(full_model()[["test_res"]], aggregate = FALSE)
  })
  
  output[["test_summ"]] <- renderPrint({
    summary(full_model()[["test_res"]])
  })
  
  output[["test"]] <- renderPrint({
    summary(full_model()[["dat"]])
  })
})