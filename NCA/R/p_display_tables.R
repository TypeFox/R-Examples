p_display_tables <- 
function (results.nca, loop.data, use.title=TRUE, pdf=FALSE) {
  # Create the three tables
  data1  <- c(loop.data$x.low, loop.data$x.high, 
              loop.data$y.low, loop.data$y.high)
  remainder <- sprintf("%f", sum(data1) %% 1)
  digits <- nchar(sub("0+$", "", remainder)) - 2
  
  # Properties of the dataset
  num1 <- p_pretty_number(min(length(loop.data$x), length(loop.data$y)), "", 0, TRUE)
  num2 <- p_pretty_number(loop.data$scope, "", "auto", TRUE)
  num3 <- matrix(sapply(data1, p_pretty_number, "", digits, TRUE), nrow=4)
  table1 <- rbind(num1, num2, num3)
  colnames(table1) <- ""
  rownames(table1) <- p_TABLE1ROWS
  
  # Effect per alalysis type
  data <- c(sapply(results.nca, function(x) p_pretty_number(x$ceiling)),
            sapply(results.nca, function(x) p_pretty_number(x$effect)),
            sapply(results.nca, function(x) p_pretty_number(x$above, prec=0)))
  table2 <- matrix(data, ncol=3)
  
  # Add col/row names, also remove the OLS row
  rownames(table2) <- names(results.nca)
  m <- match("ols", names(results.nca), 0)
  if (m > 0) {
    table2 <- matrix(table2[-m,], ncol=3)
    rownames(table2) <- names(results.nca)[-m]
  }
  colnames(table2) <- p_TABLE2COLS
  
  # Efficiencies
  data  <- c(sapply(results.nca, function(x) p_pretty_number(x$slope)),
             sapply(results.nca, function(x) p_pretty_number(x$intercept)),
             sapply(results.nca, function(x) p_pretty_number(x$ineffs$abs)),
             sapply(results.nca, function(x) p_pretty_number(x$ineffs$rel)),
             sapply(results.nca, function(x) p_pretty_number(x$ineffs$x)),
             sapply(results.nca, function(x) p_pretty_number(x$ineffs$y)))
  table3 <- matrix(data, ncol=6)
  
  # Add col/row names, also remove the OLS row
  rownames(table3) <- names(results.nca)
  m <- match("ols", names(results.nca), 0)
  if (m > 0) {
    table3 <- matrix(table3[-m,], ncol=6)
    rownames(table3) <- names(results.nca)[-m]
  }
  colnames(table3) <- p_TABLE3COLS

  # If not PDF, create new graphics
  title = paste("NCA Parameters: ", p_generate_title(loop.data))
  if (!pdf) {
    p_new_window(title=title)
  }
  par(mfrow=c(3, 1))
  par(family="mono")
  
  # Plot the three tables
  textplot(table1, cex=1.5, halign="left", mar=c(0,0,4,0), rmar=1)
  if (use.title) {
    par(family="")
    title(title, cex.main=1.5)
    par(family="mono")
  }
  if (nrow(table2) > 0) {
    textplot(table2, cex=1.5, halign="left", mar=c(0,0,0,0), rmar=1)
  }
  if (nrow(table3) > 0) {
    textplot(table3, cex=1.5, halign="left", mar=c(0,0,0,0), rmar=1)
  }
  
  if (pdf) {
    dev.off()
  }
}
