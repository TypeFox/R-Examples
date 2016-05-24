### show method for glimQL objects
setMethod("show", signature(object = "glimQL"),
  function(object){
    cat("Quasi-likelihood generalized linear model\n")
    cat("-----------------------------------------\n")
    print(object@CALL)
    cat("\nFixed-effect coefficients:\n")
    b  <- coef(object)
    se <- sqrt(diag(vcov(object)))
# function to insert NAs in se when there are NAs in b
    insNA <- function(b, se){
    if(any(is.na(b))){
      nb <- length(b)
      SE <- rep(NA, nb)
      j <- 1
      for(i in seq(nb))
        if(!is.na(b[i])){
          SE[i] <- se[j]
          j <- j + 1
          }
      se <- SE
      }
    se
    }
# model coefficients
    se <- insNA(b, se)
    summry <- data.frame(b = b, se = se, z = b / se, P = 2 * (1 - pnorm(abs(b / se))))
# model display
    nam <- rownames(summry)
    List <- vector(mode = "list", length = 4)
    for(i in 1:4){
      x <- summry[,i]
      List[[i]] <- if(i < 4)
                     format(round(x, digits = 4), nsmall = 3)
                   else
                     ifelse(x < 1e-4, "< 1e-4", format(round(x, digits = 4), nsmall = 3))
      }
    summ <- as.data.frame(t(do.call("rbind", List)))
    rownames(summ) <- nam
    colnames(summ) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    print(summ)
    cat("\nOverdispersion parameter:\n")
    print(c(phi = round(object@phi, digits = 4)))
    cat("\nPearson's chi-squared goodness-of-fit statistic =",
        round(sum(residuals(object, type = "pearson")^2), digits = 4), "\n\n")
    invisible(summry)
    })
