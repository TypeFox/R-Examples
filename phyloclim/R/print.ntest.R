print.ntest <- function(x, ...){
  
  ## format for pretty output
  ## ------------------------
  fillup <- function(x){
    x <- as.character(round(x, 3))
    n <- 5 - nchar(x)
    if (n > 0) x <- paste(x, rep(0, n), sep = "")
    x
  }
  
  if (x$method == "niche identity test"){
    cat(x$method, "\n")
    cat("species:", x$species, "\n\n")
    cat("D = ", fillup(x$statistic[1]), ", p = ", x$p.value[1], "\n", sep = "")
    cat("I = ", fillup(x$statistic[2]), ", p = ", x$p.value[2], "\n", sep = "")
    cat("(p-value based on", nrow(x$null.dist) + 1, "permutations)\n")
    cat("alternative hypothesis:", x$null)
  }
  if (x$method == "background similarity test"){
    
    conf.level <- gsub("%", "", rownames(x$ci.x.randomY))
    conf.level <- diff(as.numeric(conf.level))
    conf.level <- paste(conf.level, "%", sep = "")
    
    cat(x$method, "\n")
    cat("species:", paste(c("(X)", "  (Y)"), x$species), "\n\n")
    cat(paste(rep(" ", 6), sep = ""), paste("-----", conf.level, 
                                            "confidence intervals -----\n"))
    cat(paste(rep(" ", 6), sep = ""), "X vs. random Y      Y vs. random X\n")
    
    xy <- apply(x$ci.x.randomY, 1:2, fillup)
    yx <- apply(x$ci.y.randomX, 1:2, fillup)
    
    cat("D = ", fillup(x$statistic[1]), "   ( ", 
        paste(xy[, 1], collapse = ", "), " )    ( ", 
        paste(yx[, 1], collapse = ", "), " )\n", sep = "")
    cat("I = ", fillup(x$statistic[2]),  "   ( ", 
        paste(xy[, 2], collapse = ", "), " )    ( ", 
        paste(yx[, 2], collapse = ", "), " )\n", sep = "")
    cat("(confidence intervals based on", nrow(x$nd.x.randomY) + 1, "permutations)\n")
    cat("alternative hypothesis:", x$null)
  }
}