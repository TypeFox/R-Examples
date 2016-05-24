##' @export
print.summary.biprobit <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  printCoefmat(x$coef,digits=digits,...)
  S <- x$score;  names(S) <- rep("",length(S))
  cat("\nlogLik:", sum(x$logLik));
  cat("  mean(score^2):", formatC(mean(S^2),...), "\n");
  print(x$N,quote=FALSE)
  if (!is.null(x$msg)) {
      cat(x$msg,"\n")
  }

  for (i in seq(x$ncontrasts)) {
      corcontr <- x$model$zlen>1
      mcontr2 <- !x$model$eqmarg && x$model$blen>2
      mcontr1 <- x$model$eqmarg & x$model$blen>1
      if (corcontr || mcontr2 || mcontr1 || x$contrast) cat("\nContrast:\n")
      if (corcontr || x$contrast) {
          cat("\tDependence   ", x$par[[i]]$corref, "\n")
      }
      if (mcontr2 || (x$contrast & !x$model$eqmarg)) {
          cat("\tMean 1       ", x$par[[i]]$mref1, "\n")
          cat("\tMean 2       ", x$par[[i]]$mref2, "\n")
      } 
      if (mcontr1 || (x$contrast & x$model$eqmarg))
          cat("\tMean         ", x$par[[i]]$mref1, "\n")

      if (!is.null(x$varcomp)) {
          cat("\n")
          ##res <- x$varcomp
          res <- c()
          P <- x$prob[seq(x$nstat)+x$nstat*(i-1),,drop=FALSE]
          if (!is.null(P)) {
              res <- rbind(res,P)
          }
          idx <- unlist(sapply(c("Concordance","Marginal","P\\(Y"),function(x) grep(x,rownames(res))))
          idx2 <- setdiff(seq(nrow(res)),idx)
          res2 <- rbind(res[idx2,],rep(NA,ncol(res)),res[idx,])
          nn <- rownames(res2)
          rownames(res2) <- unlist(lapply(nn,function(x) gsub(paste("c",i,":",sep=""),"",x)))
          print(RoundMat(res2,digits=digits,na=FALSE),quote=FALSE)
      }

  }
  if (!is.null(x$time)) {
      cat("\n")
      cat("Event of interest before time ", x$time, "\n", sep="")
  }
}
