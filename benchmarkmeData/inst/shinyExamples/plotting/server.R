library("benchmarkmeData")
library("ggplot2")

get_data = function(past_results, test) {
  past_results = past_results[past_results$test_group %in% test,]
  past_results = past_results[order(past_results$time),]
  past_results$rank = 1:nrow(past_results)
  past_results  
}

get_plot = function(past_results, results,
                    facet_x="None", facet_y = "None", 
                    test="prog") {
  past_results = get_data(past_results, test)
  g = ggplot(past_results) + 
    geom_point(aes(rank, time)) + 
    scale_y_log10() + 
    theme_bw() + 
    xlab("Rank") + 
    ylab("Time (secs)")
  
  g1 = g + get_facet(facet_x, facet_y)
  
  if(!is.null(results)) {
    results = results[results$test_group %in% test,]
    rank = which(past_results$time > results$time)[1] - 1/2
    if(is.na(rank)) rank = nrow(past_results) + 1
    results$rank = rank
    g1 = g1 + geom_point(data=results, aes(rank, time), 
                         colour="red", size=3)
  }
  g1
}

# Change name on ui.R to colname
lookup = function(facet_name) {
  to = c(".", "byte", "r_minor", "sysname", "blas")
  from = c("None", "Byte Optimised", "R Version", "OS", "BLAS Optimised")
  to[which(facet_name == from)]
}

## ui -> ggplot2. Ignore None
get_facet = function(facet_x, facet_y) {
  if(facet_x == "None" && facet_y == "None") return(NULL)
  if(facet_x == facet_y) facet_y = "None"
  x = paste(lookup(facet_y), lookup(facet_x), sep= " ~ ")
  facet_grid(x)
}

server = function(input, output){
  tmp_env = new.env()
  data(past_results, package="benchmarkmeData", envir=tmp_env)
  past_results = tmp_env$past_results
  test = past_results$test_group[1]
  
  past_results$byte = ifelse(past_results$byte_optimize > 0.5, "Byte Optimised", "Non-Byte Optimised")
  past_results$blas = ifelse(past_results$blas_optimize, "BLAS Optimised", "Standard")
  res = past_results
  
  results = get("results", envir=benchmarkmeData:::.bme_env)
  if(!is.null(results)) {
    results$byte = ifelse(results$byte_optimize > 0.5, "Byte Optimised", "Not Byte Optimised")
    results$blas = ifelse(results$blas_optimize, "BLAS Optimised", "Standard")
    test = results$test_group[1]
  }
  
  rv = reactiveValues(plot = get_plot(res, results, test=test), 
                      facet_x = "None", facet_y = "None")
  
  observeEvent(input$facet_x, 
               {rv$plot <- get_plot(res, results, input$facet_x, input$facet_y, input$test)})
  observeEvent(input$facet_y, 
               {rv$plot <- get_plot(res, results, input$facet_x, input$facet_y, input$test)})
  observeEvent(input$test, 
               {rv$plot <- get_plot(res, results, input$facet_x, input$facet_y, input$test)})
  
  output$plot = renderPlot(rv$plot)
}

