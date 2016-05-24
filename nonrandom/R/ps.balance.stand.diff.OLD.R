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
  
  ## test statistics, p-values
  ##if (nlevels(as.factor(index)) != 2){
  if (!match.T){
    ## stratification
    means0 <- means1 <-
      sd0 <- sd1 <- stdf <- matrix(NA,
                                   ncol=dim(data)[2],
                                   nrow=nlevels(as.factor(index))+1)
  }else{
    ## matching
    means0 <- means1 <-
      sd0 <- sd1 <- stdf <- matrix(NA,
                                   ncol=dim(data)[2],
                                   nrow=nlevels(as.factor(index)))
  }
  colnames(means1) <- colnames(means0) <- names(data)
  colnames(sd0) <- colnames(sd1) <- colnames(stdf) <- names(data)


  
  ## #############
  ## Loop over sel
  for (i in 1:dim(data)[2]){
    
    cov <- as.numeric(data[,i])
 
    if (nlevels(as.factor(cov))   == 1 |
        nlevels(as.factor(treat)) == 1){

      meth[i] <- "none"
      table.before[i] <- table.after[i] <- NA
      next
   
    }else{

      
      ## binary/categorical covariates
      if (nlevels(as.factor(cov)) == 2){

        meth[i]  <- "bin"
        
        if (!match.T){ ## if stratification

          if (max(cov, na.rm=TRUE) > 1){ ## na.rm=TRUE is necessary if
                                         ## cov contains NAs.
            help <- cov
            help <- ifelse(cov == max(cov, na.rm=TRUE), 1 , 0)
            cov  <- help
          }

          table.cov <- table(cov, treat)
        
          ## if only treated/untreated observations or only one category
          ## of cov is assigned, then no test is available        
          if (dim(table.cov)[2] != nlevels(as.factor(treat)) |
              dim(table.cov)[1] != nlevels(as.factor(cov))){           

            means0[1,i] <- means1[1,i]     <- NA
            sd0[1,i]    <- sd1[1,i]        <- NA
            stdf[1,i]   <- table.before[i] <- NA

          }else{

            means0[1,i] <- mean(cov[treat == min(treat, na.rm=TRUE)], na.rm=TRUE)
            means1[1,i] <- mean(cov[treat == max(treat, na.rm=TRUE)], na.rm=TRUE)

            sd0[1,i] <- sqrt(means0[1,i]*(1-means0[1,i]))
            sd1[1,i] <- sqrt(means1[1,i]*(1-means1[1,i]))

            if (equal){
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((sd1[1,i]^2 + sd0[1,i]^2) / 2))
            }else{

              ## Change 02/10/2012: switch n0 and n1
              n0 <- length(cov[treat == min(treat, na.rm=TRUE)])
              n1 <- length(cov[treat == max(treat, na.rm=TRUE)])
              
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((n1/sum(table.cov))*sd1[1,i]^2 +
                                     (n0/sum(table.cov))*sd0[1,i]^2))
            }
            
            table.before[i]   <- ifelse(stdf[1,i] > alpha, 0, 1)
            ## 0: sign, 1: not sign

          }
          
          ## strata internal  
          table.cov.s <- table(cov, treat, index)
          help.NA <- rep(0,length = nlevels(as.factor(index)))
          
          for (j in 1:nlevels(as.factor(index))){
            
            ## if only treated/untreated observations or not all
            ## categories of cov in strata j are assigned, then no test
            ## is performed
            if (any(rowSums(table.cov.s[,,j]) == 0) |
                any(colSums(table.cov.s[,,j]) == 0)){
              
              means0[j+1,i] <- means1[j+1,i] <- NA
              sd0[j+1,i] <- sd1[j+1,i] <- stdf[j+1,i] <- NA
              
            }else{

              means0[j+1,i] <-
                mean(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              means1[j+1,i] <-
                mean(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              sd0[j+1,i] <- sqrt(means0[j+1,i]*(1-means0[j+1,i]))
              sd1[j+1,i] <- sqrt(means1[j+1,i]*(1-means1[j+1,i]))

              if (equal){
                stdf[j+1,i] <- 100*(abs(means1[j+1,i] - means0[j+1,i]) /
                                    sqrt((sd1[j+1,i]^2 + sd0[j+1,i]^2) / 2))
              }else{
                
                ## Change 02/10/2012: swith n0 and n1
                n0 <- length(cov[treat == min(treat, na.rm=TRUE) &  index == j])
                n1 <- length(cov[treat == max(treat, na.rm=TRUE) &  index == j])
                
                stdf[j+1,i] <- 100*(abs(means1[j+1,i] - means0[j+1,i]) /
                                    sqrt((n1/sum(table.cov.s[,,j]))*sd1[j+1,i]^2 +
                                         (n0/sum(table.cov.s[,,j]))*sd0[j+1,i]^2))
              }
            }

            if (sum(table.cov.s[,,j]) == 0){
              help.NA[j] <- j
            } 
            
          }

          help.NA <- which(help.NA != 0)+1
          
          table.after[i] <- ifelse(sum(stdf[-c(1, help.NA),i] > alpha) == 0, 1, 0)
          ## 0: at least one significant difference in strata, 1: no significance
          ## NA: in at least one stratum no test was available 

          
        }else{ ## if matching
 
          if (max(cov, na.rm=TRUE) > 1){ ## re-coded from min/max to 0/1
            help <- cov
            help <- ifelse(cov == max(cov, na.rm=TRUE), 1 , 0)
            cov  <- help
          }
                   
          ## matching internal   
          table.cov.s <- table(cov, treat, index)
          
          for (j in 1:nlevels(as.factor(index))){
            
            ## if only treated/untreated observations or not all
            ## categories of cov in strata j iare assigned, then no test is
            ## performed
            if (any(rowSums(table.cov.s[,,j]) == 0) |
                any(colSums(table.cov.s[,,j]) == 0)){

              means0[j,i] <- means1[j,i] <- NA
              sd0[j,i] <- sd1[j,i] <- NA
              stdf[j,i] <- table.before[i] <- NA
              
            }else{

              means0[j,i] <-
                mean(cov[treat == min(treat, na.rm=TRUE) & index == j], na.rm=TRUE)
              means1[j,i] <-
                mean(cov[treat == max(treat, na.rm=TRUE) & index == j], na.rm=TRUE)

              sd0[j,i] <- sqrt(means0[j,i]*(1-means0[j,i]))
              sd1[j,i] <- sqrt(means1[j,i]*(1-means1[j,i]))

              if (equal){
                stdf[j,i] <- 100*(abs(means1[j,i] - means0[j,i]) /
                                  sqrt((sd1[j,i]^2 + sd0[j,i]^2) / 2))
              }else{
                ## Change 02/10/2012: swith n0 and n1
                n0 <- length(cov[treat == min(treat, na.rm=TRUE) &  index == j])
                n1 <- length(cov[treat == max(treat, na.rm=TRUE) &  index == j])
                
                stdf[j,i] <- 100*(abs(means1[j,i] - means0[j,i]) /
                                  sqrt((n1/sum(table.cov.s[,,j]))*sd1[j,i]^2 +
                                       (n0/sum(table.cov.s[,,j]))*sd0[j,i]^2))
              }
            }
          }
          table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)
          table.after[i]  <- ifelse(stdf[2,i] > alpha, 0, 1)
          ## 0: significance, 1: no significance
        }
        
      }else{ ## continuous covariates

        ##meth[i] <- "non-bin"
        meth[i] <- "num"
        
        if (!match.T){ ## stratification
          
          ## before stratification
          if (any(length(na.omit(cov[treat==min(treat, na.rm=TRUE)])) == c(0,1)) |
              any(length(na.omit(cov[treat==max(treat, na.rm=TRUE)])) == c(0,1))){


            means0[j,i] <- means1[j,i] <- NA
            sd0[j,i] <- sd1[j,i] <- NA
            stdf[j,i] <- table.before[i] <- NA
            
          }else{

            means0[1,i] <- mean(cov[treat == min(treat, na.rm=TRUE)], na.rm=TRUE)
            means1[1,i] <- mean(cov[treat == max(treat, na.rm=TRUE)], na.rm=TRUE)

            sd0[1,i] <- sd(cov[treat == min(treat, na.rm=TRUE)], na.rm=TRUE)
            sd1[1,i] <- sd(cov[treat == max(treat, na.rm=TRUE)], na.rm=TRUE)

            if (equal){
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((sd1[1,i]^2 + sd0[1,i]^2) / 2))
            }else{
              ## Change 02/10/2012: swith n0 and n1
              n0 <- length(cov[treat == min(treat, na.rm=TRUE)])
              n1 <- length(cov[treat == max(treat, na.rm=TRUE)])
              
              stdf[1,i] <- 100*(abs(means1[1,i] - means0[1,i]) /
                                sqrt((n1/length(cov))*sd1[1,i]^2 +
                                     (n0/length(cov))*sd0[1,i]^2))
            }
            
            table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)

          }
          
          ## strata internal
          help.NA <- rep(0, nlevels(as.factor(index)))
          
          for (j in 1:nlevels(as.factor(index))){

            if (length(cov[index==j]) == 0){
              help.NA[j] <- j
            } 
            
            if (length(cov[treat == min(treat, na.rm=TRUE) & index == j]) == 0 |
                length(cov[treat == max(treat, na.rm=TRUE) & index == j]) == 0){

              means0[j+1,i] <- means1[j+1,i] <- NA
              sd0[j+1,i] <- sd1[j+1,i] <- stdf[j+1,1] <- NA
 
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
                ## Change 02/10/2012: swith n0 and n1
                n0 <- length(cov[treat == min(treat, na.rm=TRUE) & index == j])
                n1 <- length(cov[treat == max(treat, na.rm=TRUE) & index == j])
                
                stdf[j+1,i] <- 100*(abs(means1[j+1,i] - means0[j+1,i]) /
                                    sqrt((n1/length(cov[index == j]))*sd1[j+1,i]^2 +
                                         (n0/length(cov[index == j]))*sd0[j+1,i]^2))
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

              means0[j,i] <- means1[j,i] <- NA
              sd0[j,i] <- sd1[j,i] <- stdf[j,1] <- NA

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
                ## Change 02/10/2012: swith n0 and n1
                n0 <- length(cov[treat == min(treat, na.rm=TRUE) & index == j])
                n1 <- length(cov[treat == max(treat, na.rm=TRUE) & index == j])
                
                stdf[j,i] <- 100*(abs(means1[j,i] - means0[j,i]) /
                                  sqrt((n1/length(cov[index == j]))*sd1[j,i]^2 +
                                       (n0/length(cov[index == j]))*sd0[j,i]^2))
              }
            }
            table.before[i] <- ifelse(stdf[1,i] > alpha, 0, 1)
            table.after[i]  <- ifelse(stdf[2,i] > alpha, 0, 1)
            ## 0: at least one significant difference in strata, 1: no significance 
          }
        } 
      } 
    } 
  }
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
                   means0                = means0,
                   means1                = means1,
                   sd0                   = sd0,
                   sd1                   = sd1,
                   stdf                  = stdf,
                   method                = meth,
                   alpha                 = alpha)

  names(bal.list)[6:10] <-c(paste("Means.treat.",min(treat, na.rm=TRUE), sep=""),
                            paste("Means.treat.",max(treat, na.rm=TRUE), sep=""),
                            paste("SDs.treat.",min(treat, na.rm=TRUE), sep=""),
                            paste("SDs.treat.",max(treat, na.rm=TRUE), sep=""),
                            paste("Standardized.differences",sep="")) 
  
  return(bal.list)
}
