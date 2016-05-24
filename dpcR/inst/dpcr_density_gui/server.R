library(shiny)
library(binom)

shinyServer(function(input, output) {
  # Defines the input for the total number of positive molecules (k) and the total number of 
  # partitions (n).
  output$summary <- renderPrint({
    cat(paste0("Positive partitions: ", input$k))
    cat(paste0("\nTotal partitions: ", input$n))
  })
  
  # Calculates the concentration of the molecules and the total volume.
  output$volume <- renderPrint({
    cat(paste0("Concentration: ", round(- log(1 - input$k/input$n) / input$vol, 4), 
               " [molecules/nL]"))
    cat(paste0("\nTotal volume: ", input$vol * input$n, " [nL]"))
  })
  
  # Generates a plot of the number of positive partitions depending on the total number of positive 
  # molecules (k) and the total number of partitions (n).
  # NOTE: binom.confint uses different expression for n and k. 
  # k is assigned to x and n to n in here.
  # The confidence level is the only parameter that can be chosen.
  
  output$distr <- renderPlot({
    dat <- dpcR:::dpcr_calculator(input$k, input$n)
    conf <- binom.confint(input$k, input$n, methods = input$methods, 
                          conf.level = input$conf.level)[, 4:6] * input$n
    dpcR:::plot_distr(dat, ylab = "Density", xlab = "Positive partitions",
               main = "Number of positive partitions")
    dpcR:::plot_conf_int(conf, dat, "left", conf_int_col = adjustcolor("cyan4", alpha = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha = 0.15))
    dpcR:::plot_conf_int(conf, dat, "right", conf_int_col = adjustcolor("cyan4", alpha = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha = 0.15))
  })

  # Generates a plot of the number of molecules per partition depending on the total number of 
  # positive molecules (k) and the total number of partitions (n).
  # NOTE: binom.confint uses different expression for n and k. 
  # k is assigned to x and n to n in here.

  output$pois1 <- renderPlot({
    dat <- dpcR:::dpcr_calculator(input$k, input$n, average = TRUE)
    conf <- - log(1 - binom.confint(input$k, input$n, methods = input$methods, 
                          conf.level = input$conf.level)[, 4:6])
    #is that correct? || Methinks YES, (Stefan)
    dpcR:::plot_distr(dat, ylab = "Density", xlab = "Molecules/partition", 
               main = "Number of molecules per partition")
    dpcR:::plot_conf_int(conf, dat, "left", conf_int_col = adjustcolor("cyan4", alpha = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha = 0.15))
    dpcR:::plot_conf_int(conf, dat, "right", conf_int_col = adjustcolor("cyan4", alpha = 0.15), 
                  conf_int_border = adjustcolor("cyan4", alpha = 0.15))
  })

  # Generates an output of all Binomial confidence intervals depending on the positive 
  # molecules (k) and the total number of partitions (n).
  # NOTE: binom.confint uses different expression for n and k. 
  # k is assigned to x and n to n in here.

  output$table <- renderTable({
    conf_all <- binom.confint(input$k, input$n, methods = "all", 
                          conf.level = input$conf.level)
    
    conf_all[, 4:6] <- conf_all[, 4:6] * input$n
    levels(conf_all$method)[levels(conf_all$method) == 
                              input$methods] <- paste0("*", toupper(input$methods), "*")
    symmetry <- round(conf_all[, 6], 4) - 
      round(conf_all[, 4], 4) == round(conf_all[, 4], 4) - round(conf_all[, 5], 4)
    dat <- data.frame(conf_all[, -c(2:3)], - log(1 - conf_all[, 4:6]/input$n), symmetry)
    names(dat) <- c("Method", "Number of positive mol. (PM)", "Lower PM", "Upper PM", 
                    "Avg. mol./partition (AM/P)", "Lower AM/P", "Upper AM/P", "symmetric")
    dat
  })
})