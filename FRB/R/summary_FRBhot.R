
summary.FRBhot <- function(object, digits=5, ...) {
  res <- list(object=object, digits=digits)
  class(res) <- "summary.FRBhot"
  res
}

print.summary.FRBhot <- function(x, ...) {

  digits <- x$digits
  x <- x$object

  if (is.null(x$Mu1))
    {cat(x$method, "\n\n")
     #cat("data: ", deparse(x$data),"\n")
     cat("data: ", x$data.name,"\n")
     cat("T^2_R = ",round(x$statistic,2), "\n")
     cat("p-value = ", round(x$p.value,4), "\n")
     cat("Alternative hypothesis :",x$alt,"\n")
     cat("\n", x$conf*100, "% simultaneous confidence intervals for components of mean :\n")
     print(x$CI,digits=digits)
     cat("\nSample estimates :\n")
     cat("   location:\n")
     print(x$estimate, digits=digits)
     cat("\n   covariance:\n")
     print(x$Sigma, digits=3)
  }
  else { cat(x$meth, "\n\n")
#       cat("data: ", deparse(x$data[[1]]), " and ", deparse(x$data[[2]]),"\n")
       cat("data: ", x$data.name,"\n")
       cat("T^2_R = ",round(x$statistic,2),  "\n")
       cat("p-value = ", round(x$p.value,4), "\n")
       cat("Alternative hypothesis :",x$alt, "\n")
       cat("\n", x$conf*100, "% simultaneous confidence intervals for components of difference in means :\n")
       print(x$CI,digits=3)
       cat("\nSample estimates :\n")
       cat("   location:\n")
       printmat <- rbind(x$Mu1,x$Mu2,x$Mu1-x$Mu2)
       rownames(printmat) <- c("   Sample 1","   Sample 2", " difference")
       print(printmat, digits=digits)
       cat("\n   covariance:\n")
       print(x$Sigma, digits=3)
    }
}

