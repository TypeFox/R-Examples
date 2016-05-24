summary.ivbma <- function(object,nms.U=NULL,nms.V=NULL,...)
  {
    x <- object
    p.U <- dim(x$lambda.bar)[1]
    r <- dim(x$lambda.bar)[2]
    p.V <- length(x$rho.bar)
    tbl.s1 <- list()
    tbl.s2 <- matrix(0,p.V,5)

    
    for(j in 1:r)
      {
        tbl.s1[[j]] <- matrix(0,p.U,5)
        for(i in 1:p.U)
          {
            tbl.s1[[j]][i,] <- c(x$M.bar[i,j],x$lambda.bar[i,j],quantile(x$lambda[i,j,],c(.025,.5,.975)))
          }
        if ( !is.null (nms.U))
          {
            rownames(tbl.s1[[j]]) <- nms.U
          }
        colnames(tbl.s1[[j]]) <- c("Prob","Mean", "Lower","Med","Upper")
      }


    for(i in 1:p.V)
      {
        tbl.s2[i,] <- c(x$L.bar[i], x$rho.bar[i],quantile(x$rho[,i],c(.025,.5,.975)))
      }
    if ( !is.null (nms.V))
      {
        rownames(tbl.s2) <- nms.V
      }
    colnames(tbl.s2) <- c("Prob","Mean","Lower","Med","Upper")

    l <- NULL
    l$tbl.s1 <- tbl.s1
    l$tbl.s2 <- tbl.s2
    l$tbl.cov <- x$Sigma.bar
    if(x$run.diagnostics)
      {
        l$Bayes.Sargan <- x$Bayesian.Sargan
        l$Sargan <- x$Sargan
      }
    return(l)
  }


