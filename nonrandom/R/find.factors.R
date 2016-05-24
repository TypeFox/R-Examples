find.factors <- function(regr,
                         name.regr,
                         data,
                         lr.form,
                         name.sep.var,
                         name.resp.var)
{  
  ## ######################
  ## check regr for factors
  regr.fac <- vector(length=dim(regr)[2])
  
  for (i in 1: dim(regr)[2]){
    regr.fac[i] <- ifelse(class(regr[,i])=="factor", 1, 0)
  }  
  non.fac.matrix <- as.data.frame(regr[,which(regr.fac==0)])
  names(non.fac.matrix) <- names(regr)[which(regr.fac==0)]

  
  ## ##################################
  ## create dummy variables for factors
  fac         <- as.data.frame(regr[,regr.fac == 1])
  name.fac    <- name.regr[regr.fac == 1]
  nlevels.fac <- as.numeric(apply(as.data.frame(fac), 2,
                                  function(x) nlevels(as.factor(x))))
  levels.fac  <- list(apply(as.data.frame(fac), 2,
                            function(x) as.numeric(levels(as.factor(x)))))
  
  dummy.list <- list(); fac.matrix <- c()
  all.fac.names <- vector(length=sum(nlevels.fac)-length(nlevels.fac))

  for (i in 1:dim(fac)[2]){
    dummy.list[[i]] <- matrix(NA, nrow=dim(data)[1], ncol=(nlevels.fac[i]-1))
    
    for (j in 2:nlevels.fac[i]){ ## Kat.1 ist Referenzkat.
      dummy.list[[i]][apply(as.data.frame(fac[,i]), 1,
                            function(x) x==levels.fac[[1]][j]) == TRUE, j-1] <- 1
      dummy.list[[i]][,j-1] <- ifelse(is.na(dummy.list[[i]][,j-1]),0,1)
      dummy.list[[i]][,j-1] <- ifelse(is.na(fac[,i]),fac[,i],dummy.list[[i]][,j-1])

      if (i==1)
        all.fac.names[j-1] <- paste(name.fac[i],j-1, sep="")
      else
        all.fac.names[(nlevels.fac[i-1]-1)+(j-1)] <- paste(name.fac[i],j-1, sep="")
    }
    fac.matrix <- cbind(fac.matrix, dummy.list[[i]])
  }

  colnames(fac.matrix) <- all.fac.names
  regr1 <- cbind(fac.matrix, non.fac.matrix)

  lr.form1 <- formula(paste(name.resp.var,
                            paste(c(name.sep.var,names(regr1)),
                                  collapse="+"),
                            sep="~"))
  data1 <- cbind(data, fac.matrix)
 
  out <- list(regr1,
              data1,
              lr.form1)

  return(out)
}
