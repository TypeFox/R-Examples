balance.stand.diff <- function(sel,
                               treat,
                               index,    
                               method,
                               cat.levels,
                               match.T,
                               alpha,
                               equal)
{
  data <- sel
  
  ## Vectors for significance results before/after stratification
  table.before <- table.after <- vector(length = dim(data)[2])
  names(table.before) <- names(table.after) <- names(data)

  ## Vectors for methods used for each covariate
  meth <- vector(length=dim(data)[2])
  names(meth) <- names(data)
  
  if (!match.T){## stratification
    means0 <- means1 <- sd0 <- sd1 <- stdf <- matrix(NA, ncol=dim(data)[2],
                                                     nrow=nlevels(as.factor(index))+1)
  }else{ ## matching
    means0 <- means1 <- sd0 <- sd1 <- stdf <- matrix(NA, ncol=dim(data)[2],
                                                     nrow=nlevels(as.factor(index)))
  }
  colnames(means1) <- colnames(means0) <- colnames(sd0) <-
    colnames(sd1) <- colnames(stdf) <- names(data)


  
  ## #############
  ## Loop over sel
  for (i in 1:dim(data)[2]){
    
    cov <- as.numeric(data[,i])
    
    if (nlevels(as.factor(cov))   == 1 |
        nlevels(as.factor(treat)) == 1){ ## balance can not be computed

      meth[i] <- "none"; table.before[i] <- table.after[i] <- NA
      next
   
    }else{ ## balance can be computed

      ## (1) binary/categorical covariates
      if (nlevels(as.factor(cov)) <= cat.levels){
        
        meth[i] <- ifelse(cat.levels==2, "bin", "cat")
        
        
        if (!match.T){ ## if stratification
          
          table.cov <- table(cov, treat)
          
          if (dim(table.cov)[2] != nlevels(as.factor(treat)) |
              dim(table.cov)[1] != nlevels(as.factor(cov))){  ## if only (un)treated obs or
                                                              ## only one cov category, no tests
            means0[1,i] <- means1[1,i] <- sd0[1,i] <-
              sd1[1,i] <- stdf[1,i] <- table.before[i] <- NA

          }else{ ## tests are possible

            if (dim(table.cov)[1]>2){
                p.table.cov <- as.data.frame(prop.table(table.cov, 2)[-1,])
              }else{
                p.table.cov <- as.data.frame(t(prop.table(table.cov, 2)[-1,]))
              }
            p0 <- p.table.cov[,names(p.table.cov)==min(treat, na.rm=TRUE)]
            p1 <- p.table.cov[,names(p.table.cov)==max(treat, na.rm=TRUE)]
            
            ## covariance matrix of p0 and p1
            S.mat <- matrix(NA, nrow=length(p0), ncol=length(p0))
            
            if (equal){
              for( k in 1:dim(S.mat)[1] ){
                for (l in 1:dim(S.mat)[2] ){    
                  if (k==l){ ## diagonal elements
                    S.mat[k,l] <- ( p0[k]*(1-p0[k]) + p1[k]*(1-p1[k]) ) / 2
                  }else{
                    S.mat[k,l] <- ( p0[k]*p0[l] + p1[k]*p1[l] ) / 2 
                  }}}
              stdf[1,i] <- 100*sqrt(mahalanobis(t(p0),p1,S.mat))
            }else{ ## not equal    
              n0 <- length(cov[treat == min(treat, na.rm=TRUE)])
              n1 <- length(cov[treat == max(treat, na.rm=TRUE)])
              n  <- length(cov)
              
              for( k in 1:dim(S.mat)[1] ){
                for (l in 1:dim(S.mat)[2] ){    
                  if (k==l){
                    S.mat[k,l] <- (n0/n)*p0[k]*(1-p0[k]) + (n1/n)*p1[k]*(1-p1[k])
                  }else{
                    S.mat[k,l] <- (n0/n)*p0[k]*p0[l] + (n1/n)*p1[k]*p1[l]
                  }}}            
              stdf[1,i] <- 100*sqrt(mahalanobis(t(p0),p1,S.mat))
            }
            table.before[i]   <- ifelse(stdf[1,i] > alpha, 0, 1) ## 0: sign, 1: not sign
            
          } ## end (tests are possible)

          ## tests per stratum
          table.cov.s <- table(cov, treat, index) ## 3-dim table per stratum
          
          help.NA <- rep(0,length = nlevels(as.factor(index))) 
          
          ## loop over strata(=index)
          for (j in 1:nlevels(as.factor(index))){
           
            table.cov.s.j <- table.cov.s[,,j]
            
            if (any(rowSums(table.cov.s.j) == 0) |
                any(colSums(table.cov.s.j) == 0)){ ## no tests if only (un)treated obs 
                                                   ## or only one cov category  
              stdf[j+1,i] <- NA
              ## j+1 = 1st stratum
              
            }else{ ## tests are possible

               if (dim(table.cov.s.j)[1]>2){
                p.table.cov.s.j <- as.data.frame(prop.table(table.cov.s.j, 2)[-1,])
              }else{
                p.table.cov.s.j <- as.data.frame(t(prop.table(table.cov.s.j, 2)[-1,]))
              }
              p0.j <- p.table.cov.s.j[,names(p.table.cov.s.j)==min(treat, na.rm=TRUE)]
              p1.j <- p.table.cov.s.j[,names(p.table.cov.s.j)==max(treat, na.rm=TRUE)]
              
              ## covariance matrix of p0.j and p1.j
              S.mat.j <- matrix(NA, nrow=length(p0.j), ncol=length(p0.j))
              
              if (equal){
                for( k in 1:dim(S.mat.j)[1] ){
                  for (l in 1:dim(S.mat.j)[2] ){    
                    if (k==l){
                      S.mat.j[k,l] <-
                        ( p0.j[k]*(1-p0.j[k]) + p1.j[k]*(1-p1.j[k]) ) / 2
                    }else{
                      S.mat.j[k,l] <-
                        ( p0.j[k]*p0.j[l] + p1.j[k]*p1.j[l] ) / 2 
                    }}}
                stdf[j+1,i] <- 100*sqrt(mahalanobis(t(p0.j),p1.j,S.mat.j))
                
              }else{ ## not equal
                
                n0.j <- length(cov[treat == min(treat, na.rm=TRUE) &  index == j])
                n1.j <- length(cov[treat == max(treat, na.rm=TRUE) &  index == j])
                n.j  <- length(cov[index == j])
                
                for( k in 1:dim(S.mat.j)[1] ){
                  for (l in 1:dim(S.mat.j)[2] ){    
                    if (k==l){
                      S.mat.j[k,l] <-
                        (n0.j/n.j)*p0.j[k]*(1-p0.j[k]) + (n1.j/n.j)*p1.j[k]*(1-p1.j[k])
                    }else{
                      S.mat.j[k,l] <-
                        (n0.j/n.j)*p0.j[k]*p0.j[l] + (n1.j/n.j)*p1.j[k]*p1.j[l]
                    }}}         
                stdf[j+1,i] <- 100*sqrt(mahalanobis(t(p0.j),p1.j,S.mat.j))
              } ## end if(equal)
            } ## end (tests possible)
            
            if (sum(table.cov.s.j) == 0){help.NA[j] <- j}
            
          } ## end loop over index
          
          help.NA <- which(help.NA != 0)+1
          
          table.after[i] <- ifelse(sum(stdf[-c(1, help.NA),i] > alpha) == 0, 1, 0)
          ## 0: at least one sign diff in strata, 1: not sign
          ## NA: in at least one stratum no test was available 
          
          
        }else{ ## if matching (handled as 'per stratum')
          
          table.cov.s <- table(cov, treat, index) ## 3-dim table
          
          help.NA <- rep(0,length = nlevels(as.factor(index)))
          
          for (j in 1:nlevels(as.factor(index))){ ## j=1: before, j=2: after

            table.cov.s.j <- table.cov.s[,,j]
            
            if (any(rowSums(table.cov.s.j) == 0) |
                any(colSums(table.cov.s.j) == 0)){ ## if only (un)treated obs or
                                                   ## only one cov category is
                                                   ## present, no tests            
              stdf[j,i] <- table.before[i] <- NA
              
            }else{ ## tests are possible
              
              if (dim(table.cov.s.j)[1]>2){
                p.table.cov.s.j <- as.data.frame(prop.table(table.cov.s.j, 2)[-1,])
              }else{
                p.table.cov.s.j <- as.data.frame(t(prop.table(table.cov.s.j, 2)[-1,]))
              }
              p0.j <- p.table.cov.s.j[,names(p.table.cov.s.j)==min(treat, na.rm=TRUE)]
              p1.j <- p.table.cov.s.j[,names(p.table.cov.s.j)==max(treat, na.rm=TRUE)]
              
              ## covariance matrix of p0 and p1
              S.mat.j <- matrix(NA, nrow=length(p0.j), ncol=length(p0.j))
              
              if (equal){
                for( k in 1:dim(S.mat.j)[1] ){
                  for (l in 1:dim(S.mat.j)[2] ){    
                    if (k==l){
                      S.mat.j[k,l] <- ( p0.j[k]*(1-p0.j[k]) + p1.j[k]*(1-p1.j[k]) ) / 2
                    }else{
                      S.mat.j[k,l] <- ( p0.j[k]*p0.j[l] + p1.j[k]*p1.j[l] ) / 2 
                    }}}
                stdf[j,i] <- 100*sqrt(mahalanobis(t(p0.j),p1.j,S.mat.j))
                
              }else{ ## not equal
                
                n0.j <- length(cov[treat == min(treat, na.rm=TRUE) &  index == j])
                n1.j <- length(cov[treat == max(treat, na.rm=TRUE) &  index == j])
                n.j  <- length(cov[index == j])
                
                for( k in 1:dim(S.mat.j)[1] ){
                  for (l in 1:dim(S.mat.j)[2] ){    
                    if (k==l){
                      S.mat.j[k,l] <-
                        (n0.j/n.j)*p0.j[k]*(1-p0.j[k]) + (n1.j/n.j)*p1.j[k]*(1-p1.j[k])
                    }else{
                      S.mat.j[k,l] <-
                        (n0.j/n.j)*p0.j[k]*p0.j[l] + (n1.j/n.j)*p1.j[k]*p1.j[l]
                    }}}         
                stdf[j,i] <- 100*sqrt(mahalanobis(t(p0.j),p1.j,S.mat.j))
                
              } ## end if(equal)
            } ## end (tests possible?)
            
            if (sum(table.cov.s.j) == 0){help.NA[j] <- j}      
            
          } ## end loop over index
          
          help.NA <- which(help.NA != 0)+1
          table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)
          table.after[i]  <- ifelse(stdf[2,i] > alpha, 0, 1)
          ## 0: sign, 1: not sign
          
        } ## end (categorical variables)

        
      }else{ ## continuous covariates
          
        meth[i] <- "num"
        
        if (!match.T){ ## stratification
          
          ## before stratification
          if (any(length(na.omit(cov[treat==min(treat, na.rm=TRUE)])) == c(0,1)) |
              any(length(na.omit(cov[treat==max(treat, na.rm=TRUE)])) == c(0,1))){

            means0[j,i] <- means1[j,i] <- sd0[j,i] <-
              sd1[j,i] <- stdf[j,i] <- table.before[i] <- NA
            
          }else{

            means0[1,i] <- mean(cov[treat == min(treat, na.rm=TRUE)], na.rm=TRUE)
            means1[1,i] <- mean(cov[treat == max(treat, na.rm=TRUE)], na.rm=TRUE)

            sd0[1,i] <- sd(cov[treat == min(treat, na.rm=TRUE)], na.rm=TRUE)
            sd1[1,i] <- sd(cov[treat == max(treat, na.rm=TRUE)], na.rm=TRUE)

            if (equal){
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((sd1[1,i]^2 + sd0[1,i]^2) / 2))
            }else{
              n0 <- length(cov[treat == min(treat, na.rm=TRUE)])
              n1 <- length(cov[treat == max(treat, na.rm=TRUE)])
              
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((n1/length(cov))*sd1[1,i]^2 +
                                     (n0/length(cov))*sd0[1,i]^2))
            }
            
            table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)

          }
          
          ## per stratum
          help.NA <- rep(0, nlevels(as.factor(index)))
          
          for (j in 1:nlevels(as.factor(index))){
            if (length(cov[index==j]) == 0){help.NA[j] <- j} 
            
            if (length(cov[treat == min(treat, na.rm=TRUE) & index == j]) == 0 |
                length(cov[treat == max(treat, na.rm=TRUE) & index == j]) == 0){

              means0[j+1,i] <- means1[j+1,i] <- sd0[j+1,i] <-
                sd1[j+1,i] <- stdf[j+1,1] <- NA
 
            }else{

              means0[j+1,i] <-
                mean(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              means1[j+1,i] <-
                mean(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              sd0[j+1,i] <-
                sd(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              sd1[j+1,i] <-
                sd(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              if (equal){
                stdf[j+1,i] <- 100*(abs(means1[j+1,i] - means0[j+1,i]) /
                                    sqrt((sd1[j+1,i]^2 + sd0[j+1,i]^2) / 2))
              }else{
                n0.j <- length(cov[treat == min(treat, na.rm=TRUE) & index == j])
                n1.j <- length(cov[treat == max(treat, na.rm=TRUE) & index == j])
                
                stdf[j+1,i] <- 100*(abs(means1[j+1,i] - means0[j+1,i]) /
                                    sqrt((n1.j/length(cov[index == j]))*sd1[j+1,i]^2 +
                                         (n0.j/length(cov[index == j]))*sd0[j+1,i]^2))
              }
            }
          }
          help.NA <- which(help.NA != 0)+1
          
          table.after[i] <- ifelse(sum(stdf[-c(1, help.NA),i] > alpha) == 0, 1, 0)
          ## 0: at least one significant difference in strata, 1: no significance
          
        }else{ ## IF MATCHING
         
          ## matching internal
          for (j in 1:nlevels(as.factor(index))){
          
            if (length(cov[treat == min(treat, na.rm=TRUE) & index == j]) == 0 |
                length(cov[treat == max(treat, na.rm=TRUE) & index == j]) == 0){

              means0[j,i] <- means1[j,i] <- sd0[j,i] <-
                sd1[j,i] <- stdf[j,1] <- NA

            }else{
              
              means0[j,i] <-
                mean(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              means1[j,i] <-
                mean(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              sd0[j,i] <-
                sd(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              sd1[j,i] <-
                sd(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              if (equal){
                stdf[j,i] <- 100*(abs(means1[j,i] - means0[j,i]) /
                                  sqrt((sd1[j,i]^2 + sd0[j,i]^2) / 2))
              }else{
                n0.j <- length(cov[treat == min(treat, na.rm=TRUE) & index == j])
                n1.j <- length(cov[treat == max(treat, na.rm=TRUE) & index == j])
                
                stdf[j,i] <- 100*(abs(means1[j,i] - means0[j,i]) /
                                  sqrt((n1.j/length(cov[index == j]))*sd1[j,i]^2 +
                                       (n0.j/length(cov[index == j]))*sd0[j,i]^2))
              }
            }
            table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)
            table.after[i]  <- ifelse(stdf[2,i] > alpha, 0, 1)
            ## 0: at least one significant difference in strata, 1: no significance 

          } ## end loop over index 
        } ## end !match.T
      } ## end (categorical/continuous)
    } ## test whether balance can be computed
  } ## end loop over covariates


  
  tab <- rbind(table.before, table.after)
  colnames(tab) <- names(data)

  bal.tab <- matrix(NA,2,2)
  ## Prepare output as table           
  ##            | s before | ns before 
  ## -----------|----------|---------- 
  ##   s after  |          |           
  ##  ns after  |          |
  
  bal.tab[1,1] <- length(which(table.before == 0 & table.after ==0))
  bal.tab[2,2] <- length(which(table.before == 1 & table.after ==1))
  bal.tab[2,1] <- length(which(table.before == 0 & table.after ==1))
  bal.tab[1,2] <- length(which(table.before == 1 & table.after ==0))
  
  colnames(bal.tab) <- c("before: no balance (0)", "before: balance (1)") 
  rownames(bal.tab) <- c("after: no balance (0)" , "after: balance (1)") 
                                                                       
  cov.NA <- colnames(tab)[is.na(tab[1,]) | is.na(tab[2,])]
  cov.bal.before <- colnames(tab)[tab[1,]==1 & !is.na(tab[1,])]
  cov.bal.after  <- colnames(tab)[tab[2,]==1 & !is.na(tab[2,])]
       
  bal.list <- list(balance.table         = tab,
                   balance.table.summary = bal.tab,
                   covariates.NA         = cov.NA,
                   covariates.bal.before = cov.bal.before,
                   covariates.bal.after  = cov.bal.after,
                   ##means0                = means0,
                   ##means1                = means1,
                   ##sd0                   = sd0,
                   ##sd1                   = sd1,
                   stdf                  = stdf,
                   method                = meth,
                   alpha                 = alpha)

  names(bal.list)[6] <-paste("Stand.diff",sep="")
  
  return(bal.list)
}
