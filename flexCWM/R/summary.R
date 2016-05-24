summary.cwm <-function(object,criterion="BIC", concomitant=FALSE,digits = getOption("digits")-2, ...)
{
  criterion <- match.arg(criterion,.ICnames())
  best <- getBestModel(object,criterion=criterion,...)
  obj <- best$models[[1]]
  title1 <- paste0("Best fitted model according to ",criterion)
  
  nch <- nchar(title1)
  cat(rep("-",nch ),"\n",sep="")
  cat(title1, "\n")
  cat(rep("-", nch),"\n",sep="")
  #
  tab <- data.frame("loglikelihood" = obj$logLik, "n" = length(obj$cluster), 
                    "df" = obj$df,row.names = "")
  tab[[criterion]] <- obj$IC[[criterion]]
                    
  print(tab, digits = digits)
  #
  cat("\nClustering table:")
  print(table(obj$cluster), digits = digits)
  #
  cat("\nPrior: ")
  cat(paste(names(obj$prior), format(obj$prior,digits=digits), sep = " = ", 
            collapse = ", "), sep = "")
  cat("\n")
  #
  if (length(obj$GLModel)>0){
    cat("\n")
    cat(paste0("Distribution used for GLM: ",.getFamily(best,1),". Parameters:"))
    #
    cat("\n")
    for(i in seq_len(obj$k)){
      par <- obj$GLModel[[i]]
      cat("\n")
      cat(paste("Component",i))
      cat("\n")
      printCoefmat(coef(summary(par$model)))
      for(j in seq_len(length(par)-1)){
        .mycat(par[j+1],digits=digits)
        cat("\n")
      }
    }
    cat("\n")
  }

  # 
  if(!is.null(obj$concomitant$normal.model)){
    cat("Model for normal concomitant variables: ", as.character(obj$concomitant$normal.model), " (", .ModelNames(obj$concomitant$normal.model)$type, 
        ") with ", obj$k, ifelse(obj$k > 1, " components\n", " component\n"),
        sep="")
    cat("\n")
  }
  #

  if(concomitant & !is.null(obj$concomitant$normal.mu)){
    cat("Normal concomitant variables")
    cat("\n  Means:\n  ")
    print(obj$concomitant$normal.mu, digits = digits)
    cat("\n  Variance-covariance matrices:\n") 
    for(i in seq_len(obj$k)){
      cat(paste0("  Component ",i,"\n"))
      print(obj$concomitant$normal.Sigma[,,i], digits = digits) 
    }
    cat("\n")
  }

  #
  if(concomitant & !is.null(obj$concomitant$multinomial.prob)){
    cat("Categorical concomitant variables: multinomial probabilities")
    print((obj$concomitant$multinomial.prob), digits = digits)
    cat("\n")
  }

  if(concomitant & !is.null(obj$concomitant$poisson.lambda)){
    cat("Poisson concomitant variables: lambda parameter \n")
    print(obj$concomitant$poisson.lambda, digits = digits)
    cat("\n")
  }
if(concomitant & !is.null(obj$concomitant$binomial.p)){
  cat("Binomial concomitant variables: p parameter \n")
  print(obj$concomitant$binomial.p, digits = digits)
  cat("\n")
}
}
.mycat <- function(x,digits) cat(paste(names(x), format(x,digits=digits), sep = " = ", collapse = ", "), sep = "")

print.cwm <- function(x,...){
  if (length(x$models) >0) {
    best <- whichBest(x,...)
    best.unique <- unique(best)  
    for (i in seq_len(length(best.unique))){
      if (length(x$models)>1){
        b <- best==best.unique[i]        
        m <- paste(names(best)[b],collapse=", ")
        m <- paste("\nBest model according to", m, "is obtained with")
      } else m <- "\nEstimated model with"
      m <- paste(m, "k =", x$models[[best.unique[i]]]$k,"group(s)")
      if (!is.null(x$models[[best.unique[[i]]]]$concomitant$normal.model)) 
        m <- paste(m, "and parsimonious model",  x$models[[best.unique[[i]]]]$concomitant$normal.model)
      fam <- .getFamily(x,best.unique[[i]])
      if (!is.null(fam)) 
        m <- paste(m, "and family",.getFamily(x,best.unique[[i]]))
      
      cat(m,"\n")
    }
  }
  else cat("No models have been estimated.")
  invisible(x)
}
