fam_env <-
function(formula, BbBw = NULL, regression = 'linear', cluster = 'default', adjust='robcov', robcov_method='huber', bootcov_B = 200, data, ...){
  formula0 <- formula
  data0 <- data[,all.vars(formula0)]
  convert <- all.vars(formula0)[1 + which(all.vars(formula0)[2:length(all.vars(formula0))] %in% BbBw)]
  if(cluster=='default'){cluster0 <- rep(1:(nrow(data0)/2), each=2)}
    else {cluster0 <- cluster}
  if( length(convert) >= 1 ){
    newnamesBb <- rep("Bb.", length(convert))
    newnamesBw <- rep("Bw.", length(convert))
      data0[,(1+ncol(data0)):(ncol(data0)+2*length(convert))] <- NA
        for(j in 1:length(convert)){
        newnamesBb[j] <- paste(newnamesBb[j], convert[j], sep="")
        newnamesBw[j] <- paste(newnamesBw[j], convert[j], sep="")
          for(k in seq(1, nrow(data0), by=2)){
          data0[k, j + ncol(data0)-2*length(convert)] <- mean(c(data0[k,convert[j]], data0[k+1,convert[j]]))
          data0[k+1, j + ncol(data0)-2*length(convert)] <- data0[k, j + ncol(data0)-2*length(convert)]
          data0[k, j + ncol(data0)-length(convert)] <- data0[k,convert[j]] - data0[k, j + ncol(data0)-2*length(convert)]
          data0[k+1, j + ncol(data0)-length(convert)] <- data0[k+1,convert[j]] - data0[k+1, j + ncol(data0)-2*length(convert)]
          }
        }
      colnames(data0)[(1+ncol(data0)-2*length(convert)):(ncol(data0))] <- c(newnamesBb, newnamesBw)
    formula1 <- paste(all.vars(formula0)[1], "~", paste(setdiff(all.vars(formula0)[2:length(all.vars(formula0))], BbBw), collapse=" + "), "+", paste(c(newnamesBb, newnamesBw), collapse=" + "))
      if(regression=='linear'){unadjFIT <- ols(formula = as.formula(formula1), data = data0, x=T, y=T, ...)
      } else if(regression=='logistic'){unadjFIT <- lrm(formula = as.formula(formula1), data = data0, x=T, y=T, ...)}  
        if(adjust=="robcov"){adjFIT<-robcov(unadjFIT, cluster0, method=robcov_method)
        } else if(adjust=="bootcov"){adjFIT<-bootcov(unadjFIT, cluster0, B=bootcov_B)}
  } else{
    if(regression=='linear'){unadjFIT <- ols(formula = formula0, data = data0, x=T, y=T, ...)
    } else if(regression=='logistic'){unadjFIT <- lrm(formula = formula0, data = data0, x=T, y=T, ...)}
      if(adjust=="robcov"){adjFIT<-robcov(unadjFIT, cluster0, method=robcov_method)
      } else if(adjust=="bootcov"){adjFIT<-bootcov(unadjFIT, cluster0, B=bootcov_B)}
    }
}
