`fitdAICrc` <-
function(x, modelset = c("pureBirth", "bd", "DDL", "DDX", "yule2rate"), ints = NULL)
{
  checkbasal(x)
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  res <- list()
  
  temp <- list()
  
  for (i in 1:length(modelset))
  {
    temp <- 0
    if (modelset[i] != "yule2rate" && modelset[i] != "rbvd" && modelset[i] != "yule3rate")
        method <- parse(text = paste(modelset[i], "(x)", sep = ""))
    else
        method <- parse(text = paste(modelset[i], "(x, ints)", sep = ""))
    temp <- suppressWarnings(eval(method))
    #temp <- eval(method)

    #cat(method, temp$LH, "\n")
    res$model[i] <- modelset[i]
    if (modelset[i] == "pureBirth") {res$params[i] = "r1"
      res$np[i] <- 1
      res$mtype[i] <- "RC"}
    if (modelset[i] == "bd") {res$params[i] = "r1, a"
      res$np[i] <- 2
      res$mtype[i] <- "RC"}
    if (modelset[i] == "DDX") {res$params[i] = "r1, X"
      res$np[i] <- 2
      res$mtype[i] <- "RV"}
    if (modelset[i] == "DDL") {res$params[i] = "r1, k"
      res$np[i] <- 2
      res$mtype[i] <- "RV"}
    if (modelset[i] == "yule2rate") {res$params[i] = "r1, r2, ts"
      res$np[i] <- 3
      res$mtype[i] <- "RV"
      temp=as.list(temp)
      }
    if (modelset[i] == "rvbd") {res$params[i] = "r1, r2, a, ts"
      res$np[i] <- 4
      res$mtype[i] <- "RV"}
    if (modelset[i] == "yule3rate"){res$params[i] = "r1, r2, r3, st1, st2"
      res$np[i] <- 5
      res$mtype[i] <- "RV"
      temp=as.list(temp)
     }
    #fill res with estimated parameter values

    res$LH[i] <- temp$LH
    res$r1[i] <- temp$r1
    if (!is.null(temp$r2))
      res$r2[i] <- temp$r2
    else
      res$r2[i] <- NA
    if (!is.null(temp$a))
      res$a[i] <- temp$a
    else
      res$a[i] <- NA
    if (!is.null(temp$xp))
      res$xp[i] <- temp$xp
    else
      res$xp[i] <- NA
    if (!is.null(temp$k))
      res$k[i] <- temp$k
    else
      res$k[i] <- NA
    if (!is.null(temp$st))
      res$st[i] <- temp$st
    else
      res$st[i] <- NA
    if (!is.null(res$LH[i]))
      res$AIC[i] <- (-2*res$LH[i]) + 2*res$np[i]
    else
      res$AIC[i] <- NA
    if (!is.null(temp$st2))
      res$st2[i] <- temp$st2
    else
      res$st2[i] <- NA
    if (!is.null(temp$r3))
      res$r3[i] <- temp$r3
    else
      res$r3[i] <- NA
  }
  aicmin <- min(res$AIC)
  res$dAIC <- res$AIC - aicmin
  rcbest <- min(res$AIC[res$mtype == "RC"])
  rvbest <- min(res$AIC[res$mtype == "RV"])
  daicstat <- rcbest - rvbest


  cat("\n--------------Model Summary----------------\n\n")
  pvec <- c("LH", "AIC", "r1", "r2", "r3", "a", "x", "k", "st", "st2")
  for (i in 1:length(modelset))
  {
    cat("MODEL", res$model[i], "\n\nParameters: ", res$params[i], "\n\n")
    for (j in 1:length(pvec))
    {
     if (!is.na(eval(parse(text = paste("res$", pvec[j], "[i]", sep = "")))))
        cat(pvec[j], eval(parse(text = paste("res$", pvec[j], "[i]", sep = ""))), "\n\n")

    }
    cat("\n--------------------------\n")

  }

  cat("\nBest Constant Rate Model =", res$model[res$AIC == min(res$AIC[res$mtype == "RC"])], " AIC ", rcbest, "\n")
  cat("\nBest Rate Variable Model =", res$model[res$AIC == min(res$AIC[res$mtype == "RV"])], " AIC ", rvbest, "\n")
  cat("\ndelta AICrc = ", daicstat, "\n")
  res <- as.data.frame(res)
  return(res)

}

