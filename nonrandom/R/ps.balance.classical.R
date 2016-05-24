balance.classical <- function(sel,
                              treat,
                              index,    
                              method,
                              cat.levels,
                              match.T,
                              alpha,
                              equal)
{
  data  <- sel
  
  ## Vectors for significance results before/after stratification
  Before <- After <- vector(length=dim(data)[2])
  names(Before) <- names(After) <- names(data)
  
  ## Vectors for methods used for each covariate
  meth <- vector(length=dim(data)[2])
  names(meth) <- names(data)

  ## test statistics, p-values
  if (!match.T){
    ## stratification
    value.matrix <- p.matrix <- matrix(NA,
                                       ncol=dim(data)[2],
                                       nrow=nlevels(as.factor(index))+1)
  }else{
    ## matching
    value.matrix <- p.matrix <- matrix(NA,
                                       ncol=dim(data)[2],
                                       nrow=nlevels(as.factor(index)))
  }
  colnames(value.matrix) <- colnames(p.matrix) <- names(data)


  ## ###################
  ## Loop over sel
  for (i in 1:dim(data)[2]){
    
    cov <- data[,i]
 
    if ( nlevels(as.factor(cov)) == 1 ){

      ##meth[i] <- "none"
      meth[i] <- ""
      Before[i] <- After[i] <- NA
      next

    }else{
    
      ## binary/categorical covariates
      if ( nlevels(as.factor(cov)) <= cat.levels ){

        ## meth[i]  <- "cat"
        meth[i]  <- "chi^2"
        
        if (!match.T){ ## if stratification
          
          ## before stratification
          table.cov <- table(cov,treat)
          
          ## if only treated/untreated observations or only one category
          ## of cov is assigned, then no test is available        
          if (dim(table.cov)[2] != nlevels(as.factor(treat)) |
              dim(table.cov)[1] != nlevels(as.factor(cov))){
            
            test.cov <- Before[i] <- NA
            p.matrix[1,i] <- value.matrix[1,i] <- NA
            
          }else{
            
            test.cov          <- chisq.test(table.cov)
            Before[i]   <- ifelse(test.cov$p.value <= (alpha/100), 0, 1)
            ## 0: sign, 1: not sign
            p.matrix[1,i]     <- round(test.cov$p.value, 3)
            value.matrix[1,i] <- round(test.cov$stat, 3)
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
              
              test.cov.s <- p.matrix[j+1,i] <- value.matrix[j+1,i] <- NA
              
            }else{
              
              test.cov.s          <- chisq.test(table.cov.s[,,j])
              p.matrix[j+1,i]     <- round(test.cov.s$p.value, 3)
              value.matrix[j+1,i] <- round(test.cov.s$stat, 3)
            }

            if (sum(table.cov.s[,,j]) == 0){
              help.NA[j] <- j
            }            
          }

          help.NA <- which(help.NA != 0)+1
          
          After[i] <- ifelse(sum(p.matrix[-c(1, help.NA),i] <= (alpha/100)) == 0, 1, 0)
          ## 0: at least one significant difference in strata, 1: no significance
          ## NA: in at least one stratum no test was available 


          
        }else{ ## if matching
                  
          ## matching internal   
          table.cov.s <- table(cov, treat, index)
          
          for (j in 1:nlevels(as.factor(index))){
            
            ## if only treated/untreated observations or not all
            ## categories of cov in strata j iare assigned, then no test is
            ## performed
            if (any(rowSums(table.cov.s[,,j]) == 0) |
                any(colSums(table.cov.s[,,j]) == 0)){ 
              
              test.cov.s <- p.matrix[j,i] <- value.matrix[j,i] <- NA
              
            }else{
              
              test.cov.s        <- chisq.test(table.cov.s[,,j])
              p.matrix[j,i]     <- round(test.cov.s$p.value, 3)
              value.matrix[j,i] <- round(test.cov.s$stat, 3)
            }
          }
          Before[i] <- ifelse(p.matrix[1,i] <= (alpha/100), 0, 1)
          After[i]  <- ifelse(p.matrix[2,i] <= (alpha/100), 0, 1)
          ## 0: significance, 1: no significance
        }
      }else{ ## continuous covariates

        ##meth[i] <- "non-cat"
        meth[i] <- "t"
        
        if (!match.T){ ## stratification
          
          ## before stratification
          if ( any(length(na.omit(cov[treat==min(treat, na.rm=TRUE)])) == c(0,1)) |
               any(length(na.omit(cov[treat==max(treat, na.rm=TRUE)])) == c(0,1)) ){
            
            Before[i] <- p.matrix[1,i] <- value.matrix[1,i] <- NA
            
          }else{
            
            test.cov  <- t.test(cov[treat==min(treat, na.rm=TRUE)],
                                cov[treat==max(treat, na.rm=TRUE)])
            
            Before[i]   <- ifelse(test.cov$p.value <= (alpha/100), 0, 1)
            p.matrix[1,i]     <- round(test.cov$p.value, 3)
            value.matrix[1,i] <- round(test.cov$stat, 3)
            
          }
          
          ## strata internal
          help.NA <- rep(0, nlevels(as.factor(index)))
          
          for (j in 1:nlevels(as.factor(index))){

            if (length(cov[index==j]) == 0){
              help.NA[j] <- j
            } 
            
            if (any(length(na.omit(cov[treat==min(treat, na.rm=TRUE) & index==j])) == c(0,1)) |
                any(length(na.omit(cov[treat==max(treat, na.rm=TRUE) & index==j])) == c(0,1))){
              
              p.matrix[j+1,i] <- value.matrix[j+1,i] <- NA
              
            }else{
              
              test.cov.s  <-  t.test(cov[treat==min(treat, na.rm=TRUE) & index==j],
                                     cov[treat==max(treat, na.rm=TRUE) & index==j])
              
              p.matrix[j+1,i]     <- round(test.cov.s$p.value, 3)
              value.matrix[j+1,i] <- round(test.cov.s$stat, 3)
            }
          }
          help.NA <- which(help.NA != 0)+1
          
          After[i] <- ifelse(sum(p.matrix[-c(1, help.NA),i] <= (alpha/100)) == 0, 1, 0)
          ## 0: at least one significant difference in strata, 1: no significance
          
        }else{ ## IF MATCHING
                    
          ## matching internal
          for (j in 1:nlevels(as.factor(index))){
            
            if (any(length(na.omit(cov[treat==min(treat, na.rm=TRUE) & index==j])) == c(0,1)) |
                any(length(na.omit(cov[treat==max(treat, na.rm=TRUE) & index==j])) == c(0,1))){
              
              p.matrix[j,i] <- value.matrix[j,i] <- NA
              
            }else{
              
              test.cov.s  <-  t.test(cov[treat==min(treat, na.rm=TRUE) & index==j],
                                     cov[treat==max(treat, na.rm=TRUE) & index==j])
              
              p.matrix[j,i]     <- round(test.cov.s$p.value, 3)
              value.matrix[j,i] <- round(test.cov.s$stat, 3)
            }
            Before[i] <- ifelse(p.matrix[1,i] <= (alpha/100), 0, 1)
            After[i]  <- ifelse(p.matrix[2,i] <= (alpha/100), 0, 1)
            ## 0: at least one significant difference in strata, 1: no significance 
          }
        } 
      } 
    } 
  }
  tab <- rbind(Before, After)
  colnames(tab) <- names(data)
  

  bal.tab <- matrix(NA,2,2)
  ## Prepare output as table           
  ##            | s before | ns before 
  ## -----------|----------|---------- 
  ##   s after  |          |           
  ##  ns after  |          |
  
  bal.tab[1,1] <- length(which(Before == 0 & After ==0))
  bal.tab[2,2] <- length(which(Before == 1 & After ==1))
  bal.tab[2,1] <- length(which(Before == 0 & After ==1))
  bal.tab[1,2] <- length(which(Before == 1 & After ==0))
  
  colnames(bal.tab) <- c("Before: no bal (0)", "Before: bal (1)") 
  rownames(bal.tab) <- c("After: no bal (0)" , "After: bal (1)") 
                                                                       
  cov.NA <- colnames(tab)[is.na(tab[1,]) | is.na(tab[2,])]
  cov.bal.before <- colnames(tab)[tab[1,]==1 & !is.na(tab[1,])]
  cov.bal.after  <- colnames(tab)[tab[2,]==1 & !is.na(tab[2,])]

  bal.list <- list(balance.table         = tab,
                   balance.table.summary = bal.tab,
                   covariates.NA         = cov.NA,
                   covariates.bal.before = cov.bal.before,
                   covariates.bal.after  = cov.bal.after,
                   p.matrix              = p.matrix,
                   #value.matrix          = value.matrix,
                   method                = meth,
                   alpha                 = alpha)

  names(bal.list)[6] <-c(paste("p.value")) 
       

  return(bal.list)
  
}
