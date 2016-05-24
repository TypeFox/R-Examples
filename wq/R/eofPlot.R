eofPlot <-
function(x, type = c("coef", "amp"), rev = FALSE, ord = FALSE) {

  # Validate args
  type <- match.arg(type)
  num <- ncol(x$REOF) - 1

  # Plot
  if (type == "coef") {
    d1 <- x$REOF
    if (ord) d1[, "id"] <- reorder(levels(d1[, "id"]), d1[, "EOF1"])
    if (rev) d1[, 2:(num + 1)] <- -d1[, 2:(num + 1)]
    m1 <- melt(d1, id = "id")
    ggplot(m1, aes_string(x = "value", y = "id")) +
      geom_vline(xintercept = 0, colour = "red", size = 0.2) +
      geom_point(colour = "blue") +
      facet_wrap(~ variable, ncol = num) +
      labs(y = "", x = "Coefficient")
  } else {
    d1 <- x$amplitude
    if (rev) d1[, 2:num] <- -d1[, 2:num]
    d1 <- within(d1, id <- as.numeric(as.character(id)))
    m1 <- melt(d1, id = "id")
    ggplot(m1, aes_string(x = "id", y = "value")) +
      geom_hline(aes(yintercept = 0), colour = "red", size = 0.2) +
      geom_line(colour = "blue") +
      geom_point(colour = "blue") +
      facet_wrap(~ variable, nrow = num) +
      labs(x = "", y = "Amplitude")
  }
}
