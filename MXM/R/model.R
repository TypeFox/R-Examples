# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.

model = function(target, dataset, sesObject, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesObject@signatures) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
  }
  
  if ( any(is.na(dataset) ) == TRUE )
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")
    {
      dataset = apply(dataset, 2, function(x) { x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
    }else{
      for(i in 1:ncol(dataset))
      {
        if( any( is.na(dataset[,i]) ) )
        {
          xi = dataset[, i]
          if(class(xi) == "numeric")
          {                    
            xi[ which( is.na(xi) ) ] = median(xi, na.rm = TRUE) 
          }else if ( class(xi) == "factor" ) {
            xi[ which( is.na(xi) ) ] = levels(xi)[ which.max(xi) ]
          }
          dataset[, i] = xi
        }
      }
    }
  }
  
  if ( is.null(test) ) {  
    ci_test = sesObject@test
  } else ci_test = test 

  rob = sesObject@rob
  
  if ( nsignat == 1 || ( nsignat > 1 & nrow(sesObject@signatures) == 1 ) ) {
    ypografi = sesObject@selectedVars  
    
    p <- length(ypografi)
    # mat1 <- mat2 <- numeric(p)
    
   if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
     if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
       target = log( target/(1 - target) ) 
     }  

     if (rob == TRUE) {
       # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
       # mod = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi ] ), control = cont )
       mod = MASS::rlm(target ~., data = as.data.frame(dataset[, ypografi ]), maxit = 2000 )
       bic = BIC( mod )

      # for (i in 1:p) {
      #   mi <- robust::lmRob( target ~ ., data = as.data.frame( dataset[, c(1:i) ] ), control = cont )
      #   es <- fitted(mi)
      #   mat1[i] <- cor(target, es)^2
      #   mi <- robust::lmRob( target ~ ., data = as.data.frame( dataset[, -i ] ), control = cont )
      #   es <- fitted(mi)
      #   mat2[i] <- cor(target, es)^2
      # }  
      
     } else {
      mod = lm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
      
      # for (i in 1:p) {
      #   mi <- lm( target ~ ., data = as.data.frame( dataset[, c(1:i) ] ) )
      #   es <- fitted(mi)
      #   mat1[i] <- cor(target, es)^2
      #   mi <- lm( target ~ ., data = as.data.frame( dataset[, -i ] ) )
      #   es <- fitted(mi)
      #   mat2[i] <- cor(target, es)^2
      # }  
      
     }  
     
   } else if (ci_test == "testIndRQ") {
       if ( all( target>0 & target<1 ) ) {  ## are they proportions?
         target = log( target/(1 - target) ) 
       }
      mod = quantreg::rq( target ~., data = as.data.frame(dataset[, ypografi ])  )
      bic =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
      
      # for (i in 1:p) {
      #   mi <- quantreg::rq( target ~ ., data = as.data.frame( dataset[, c(1:i) ] ) )
      #   es <- fitted(mi)
      #   mat1[i] <- cor(target, es)^2
      #   mi <- quantreg::rq( target ~ ., data = as.data.frame( dataset[, -i ] ) )
      #   es <- fitted(mi)
      #   mat2[i] <- cor(target, es)^2
      # }  
      
    } else if ( ci_test == "testIndBeta" ) {
      mod = betareg::betareg( target ~ ., data = as.data.frame(dataset[, ypografi ])  )
      bic = BIC(mod)
      
      # y = log(target / ( 1- target) )
      # for (i in 1:p) {
      #   mi <- betareg::betareg( target ~ ., data = as.data.frame( dataset[, c(1:i) ] ) )
      #   es <- fitted(mi)
      #   es1 <- log( es / (1 - es) )
      #   mat1[i] <- cor(y, es1)^2
      #   mi <- betareg::betareg( target ~ ., data = as.data.frame( dataset[, -i ] ) )
      #   es <- fitted(mi)
      #   es1 <- log( es / (1 - es) )
      #   mat2[i] <- cor(y, es)^2
      # }
      
    } else if ( ci_test == "testIndSpeedglm" ) {
      if ( length( unique(target) )  == 2 || ( length( unique(target) )  == 2  & sum( floor(target) - target) == 0 ) ) {
        mod = speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
        bic =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
      } else {
        mod = speedglm::speedlm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
        bic =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
      }  
      
    } else if ( ci_test == "testIndPois") {
      if ( rob == TRUE ) {
        mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , poisson, maxit = 100 )
        bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
      } else {
        mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]) , poisson )
        bic <- BIC( mod )
      }

    } else if ( ci_test == "testIndNB" ) {
      mod = MASS::glm.nb( target ~ ., data = as.data.frame(dataset[, ypografi ])  )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi] ) )
      bic = BIC(mod)
      
    } else if ( class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1]/(target[, 1]) ) 
       } 
      mod = lm( target ~., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndCR") {
      mod = survival::coxph( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length(unique(target)) == 2 ) {
        if ( rob == TRUE ) {
          mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , binomial, maxit = 100 )
          bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
        } else { 
          mod = glm( target ~., data = as.data.frame(dataset[, ypografi ]) , binomial ) 
          bic = BIC(mod)
        }

      } else if ( is.ordered(target) == FALSE ) { 
        target = as.factor( as.numeric( as.vector(target) ) )
        mod = nnet::multinom( target ~., data = as.data.frame(dataset[, ypografi ]) , trace = FALSE )
        bic = BIC(mod)
        
      } else if ( is.ordered(target) == TRUE ) {
        mod = ordinal::clm( target ~., data = as.data.frame(dataset[, ypografi ])  )
        bic = BIC(mod)
      }
      
      
    }
    
    # if ( is.null( colnames(dataset) ) ) {
    #   names(ypografi) = paste("X", ypografi, sep = "")
    #   names(mat1) = paste("+X", ypografi, sep = "")
    #   names(mat2) = paste("-X", ypografi, sep = "")
    # } else {
    #   nama = colnames(dataset)
    #   names(ypografi) = nama
    #   names(mat1) = paste("+", nama, sep = "")
    #   names(mat2) = paste("-", nama, sep = "")
    # } 
     
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
    
  } 

  #############  more than one signatures
 
   if ( nsignat > 1 & nrow(sesObject@signatures) > 1 ) {
      
      if ( nsignat > nrow(sesObject@signatures) ) {
        nsignat = nrow(sesObject@signatures)
      }
    
    bic = numeric(nsignat)
    ypografi = sesObject@signatures[1:nsignat, ] 
    ypografi = as.matrix(ypografi)
    mod = list()
    
    for ( i in 1:nsignat ) {
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
       target = log( target/(1 - target) ) 
      }  

      if (rob == TRUE) {
        # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
        # mod[[ i ]] = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), control = cont )
        mod[[ i ]] = MASS::rlm(target~., data=as.data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000)
        bic[i] = cor( target, fitted(mod[[ i ]]) )^2
      } else {
        mod[[ i ]] = lm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if ( ci_test == "testIndRQ" ) {
        if ( all( target>0 & target<1 ) ) {  ## are they proportions?
         target = log( target/(1 - target) ) 
        }
       mod[[ i ]] = quantreg::rq( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] =  -2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * log( length(target) )

     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] = betareg::betareg( target ~ ., as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndSpeedglm" ) {
       if ( length( unique(target) )  == 2 || ( length( unique(target) )  == 2  & sum( floor(target) - target) == 0 ) ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
       } 
       
     } else if ( ci_test == "testIndPois ") {
       if ( rob == TRUE ) {
         mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson, maxit = 100 )
         bic[i] = mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
       } else {  
         mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson )
         bic[i] = BIC( mod[[ i ]] )
       }

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic = BIC(mod[[ i ]])

     } else if ( class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1] / (target[, 1]) ) 
       } 
       mod = lm( target ~.,  data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length(unique(target)) == 2 ) {
        if ( rob == TRUE ) {
          mod[[ i ]] <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial, maxit = 100 )
          bic[[ i ]] <- mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
        } else { 
          mod[[ i ]] = glm( target ~., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial ) 
          bic[[ i ]] = BIC(mod[[ i ]])
        }

       } else if ( is.ordered(target) == FALSE ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE )
         bic[i] = BIC( mod[[ i ]] )

       } else if ( is.ordered(target) == TRUE ) {
         mod[[ i ]] = ordinal::clm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
      
    }
    
    # if ( is.null( colnames(dataset) ) ) {
    #   names(ypografi) = paste("X", ypografi, sep = "")
    #   colnames(mat1) = paste("+X", ypografi, sep = "")
    #   colnames(mat2) = paste("-X", ypografi, sep = "")
    # } else {
    #   nama = colnames(dataset)
    #   names(ypografi) = nama
    #   colnames(mat1) = paste("+", nama, sep = "")
    #   colnames(mat2) = paste("-", nama, sep = "")
    # } 
    
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"    
    
  }
   
  ####### all signatures

  if ( nsignat == "all" ) { 
    ypografi = sesObject@signatures
    bic = numeric( nrow(ypografi) )
    mod = list()

    for ( i in 1:nrow(ypografi) ) {
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
        target = log( target / (1 - target) ) 
      } 
 
      if (rob == TRUE) {
        # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
        # mod[[ i ]] = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), control = cont )
        mod[[ i ]] = MASS::rlm(target~., data=as.data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000)
                bic[[ i ]] = cor( target, fitted(mod[[ i ]]) )^2
      } else {
        mod[[ i ]] = lm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if (ci_test == "testIndRQ") {
        if ( all( target > 0 & target < 1 ) ) {  ## are they proportions?
           target = log( target/(1 - target) ) 
        }
       mod[[ i ]] = quantreg::rq( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] =  -2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * log( length(target) )

     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] = betareg::betareg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndSpeedglm" ) {
       if ( length( unique(target) )  == 2 || ( length( unique(target) )  == 2  & sum( floor(target) - target) == 0 ) ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
       } 

     } else if ( ci_test == "testIndPois ") {
       if ( rob == TRUE ) {
         mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson, maxit = 100 )
         bic[i] = mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
       } else {  
         mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson )
         bic[i] = BIC( mod[[ i ]] )
       }

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic = BIC(mod)

     } else if (class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1] / (target[, 1]) ) 
       } 
       mod = lm( target ~., data = dataset[, ypografi] )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length(unique(target)) == 2) {
        if ( rob == TRUE ) {
          mod[[ i ]] <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial, maxit = 100 )
          bic[[ i ]] <- mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
        } else { 
          mod[[ i ]] = glm( target ~., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial ) 
          bic[[ i ]] = BIC(mod[[ i ]])
        }

       } else if ( is.ordered(target) == FALSE ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE )
         bic[i] = BIC( mod[[ i ]] )

       } else if ( is.ordered(target) == TRUE ) {
         mod[[ i ]] = ordinal::clm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
      
    }
    
  }  
  
    # if ( is.null( colnames(dataset) ) ) {
    #   names(ypografi) = paste("X", ypografi, sep = "")
    #   colnames(mat1) = paste("+X", ypografi, sep = "")
    #   colnames(mat2) = paste("-X", ypografi, sep = "")
    # } else {
    #   nama = colnames(dataset)
    #   names(ypografi) = nama
    #   colnames(mat1) = paste("+", nama, sep = "")
    #   colnames(mat2) = paste("-", nama, sep = "")
    # } 
    
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
  
  list(mod = mod, ypografi = ypografi)  
}
  
 