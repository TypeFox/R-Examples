sr <- function(object, ncomps = object$ncomp) {
  Coefs <- as.matrix(object$coefficients)[, ncomps]
  new.t <- as.matrix(object$Xdata) %*% ((Coefs) / sqrt(crossprod((Coefs))))
  new.p <- (t(object$Xdata) %*% new.t) / as.numeric((t(new.t) %*% new.t))
  X.hat <- new.t %*% t(new.p)
  X.error <- object$Xdata - X.hat
  n <- nrow(X.hat)
  SR.df <- data.frame(SRS = ((apply(X.hat, 2, function(x) sum((x)^2)))) / 
                        (((apply(X.error, 2, function(x) sum((x)^2))))))
  SR.df$p.value <- 1-pf(SR.df$SRS, n - 2, n - 3)
  SR.df$p.value <- 1-pf(SR.df$SRS, n - 2, n - 3)
  SR.df$f.value <- qf(.95, n - 2, n - 3)
  SR.df$Significant <- ifelse(SR.df$p.value < 0.05, "Yes", "No")
  Results <- list(sr = SR.df, modeled = X.hat, error = X.error)
  class(Results) <- "sr"
  Results
}

print.sr <- function(x, ...) {
  cat("\nsr summary:\n")
  print(x$sr[order(x$sr$SRS, decreasing = TRUE), ])
}

sr.modeled <- function(x, ...) {
  cat("\nsr model space:\n")
  print(x$modeled)
}

sr.error <- function(x, ...) {
  cat("\nsr model error:\n")
  print(x$error)
}

plot.sr <- function(x, variables = "all", ...) {
  if(variables == "all") {
    df <- x$sr[order(x$sr$SRS, decreasing = TRUE), ]
  } else {
    df <- x$sr[order(x$sr$SRS, decreasing = TRUE), ][1:variables, ]
  }
  df$Variables <- factor(row.names(df), levels = row.names(df)[order(df$SRS)])
print(with(df, ggplot(df, aes(SRS, Variables)) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          geom_point(aes(col = Significant), size = 5) + 
          ggtitle("Selectivity Ratio Results") + 
          xlab("F-statistics") + 
          # theme(legend.position = "none") + 
          theme(plot.title = element_text(size = 20)) + 
          theme(axis.title.x = element_text(size = 20)) +
          theme(axis.title.y = element_text(size = 20, angle = 90)) + 
          theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
          theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
}

