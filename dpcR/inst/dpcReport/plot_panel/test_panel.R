ppp_data <- dpcR:::create_ppp(data_vector = single_array, nx_a = ncol(single_array), 
                              ny_a = nrow(single_array), marks = TRUE, plot = FALSE)

res <- spatstat::quadrat.test(ppp_data, 
                              nx = ifelse(is.null(input[["nx"]]), 5, input[["nx"]]), 
                              ny = ifelse(is.null(input[["ny"]]), 5, input[["ny"]]),
                              "two.sided", method = "Chisq", TRUE, nsim = 1999)
