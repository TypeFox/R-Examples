#' Train ANFIS network
#'
#' ANFIS on-line or off-line hybrid Jang dynamic learning training process. In 
#' addition for off-line learning there is also adaptive learning coefficient 
#' and momentum term.
#' 
#' @param object ANFIS' class object.
#' @param A internal matrix for Iterative Least Squares Estimation of AX=B.
#' @param B internal matrix for Iterative Least Squares Estimation of AX=B.
#' @param initialGamma numeric large number >> 0. Default 1000.
#' @param epochs the max number of training epochs. Default 5.
#' @param tolerance convergence error to stop training. Default 1e-5.
#' @param k numeric with the initial step size for learning rule. Default 0.01.
#' @param eta numeric learning rule coefficient. Default 0.05.
#' @param phi numeric momentum rule coefficient. Default 0.2.
#' @param a numeric step to increase eta if delta_e is < 0, i.e. descending. 
#'  Default value 0.01.
#' @param b numeric fraction to decrease eta if delta_e is > 0, i.e. ascending. 
#'  Default value is 0.1.
#' @param delta_alpha_t_1 list with numeric matrix with last time step. Default 
#'  list().
#' @param lamda 0 < numeric < 1 forgetting factor. Default 0.9.
#' @param S covariance matrix for on-line LSE. Default matrix(nrow=0, ncol=0).
#'
#' @return 
#'  \item{matrix}{with the system solution for LSE output.}
#'  \item{error}{numeric vector with training associated errors (pattern or 
#'    epoch) according to trainingType.}
#'  \item{convergence}{TRUE/FALSE if it reached convergence or not.}
#'  \item{updated}{trainingType, premises, consequents, error, residuals, 
#'    fitted.values and coefficient.}
#'
#' @include Anfis-predict.R
#' @exportMethod LSE
#' @docType methods
#' @name LSE
#' @rdname ANFIS-training
#' @aliases LSE-methods
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setGeneric(name="LSE", def=function(object, A, B, initialGamma=1000){
  standardGeneric("LSE")
})
#'
#' @name LSE
#' @rdname ANFIS-training
#' @aliases LSE,ANFIS-method
#' @inheritParams LSE
setMethod(f="LSE", signature="ANFIS", definition=function(object, A, B, 
  initialGamma=1000){
  ## Initial conditions
  S<-initialGamma*diag(ncol(A))
  x<-matrix(0, nrow=ncol(A), ncol=ncol(B))
  
  ## Recurrence formula
  invisible(lapply(1:nrow(A), function(pattern){
    a <- t(A[pattern,,drop=FALSE])
    b <- t(B[pattern,,drop=FALSE])
    S <<- S - (S%*%a%*%t(a)%*%S)/as.numeric(1+t(a)%*%S%*%a)
    x <<- x + S%*%a%*%(t(b)-t(a)%*%x)
    return(NULL)
  }))

  return(x=x)
})
#'
#' @exportMethod trainHybridJangOffLine
#' @docType methods
#' @name trainHybridJangOffLine
#' @rdname ANFIS-training
#' @aliases trainHybridJangOffLine-methods
setGeneric(name="trainHybridJangOffLine", def=function(object, epochs=5, 
  tolerance=1e-5, initialGamma=1000, k=0.01){
    standardGeneric("trainHybridJangOffLine")
})
#'
#' @name trainHybridJangOffLine
#' @rdname ANFIS-training
#' @aliases trainHybridJangOffLine,ANFIS-method
#' @inheritParams LSE
setMethod(f="trainHybridJangOffLine", signature="ANFIS", 
  definition=function(object, epochs=5, tolerance=1e-5, initialGamma=1000, 
  k=0.01){
  ##Initialization of variables
  nameObject <- deparse(substitute(object))
  object@call <- match.call()
  object@trainingType <- "trainHybridJangOffLine"
  convergence <- FALSE
  B <- object@Y ## for LSE consequents estimation
    
  ##For each epoch
  epoch <- 1
  while(epoch<epochs & !convergence){
  print(paste("epoch: ", epoch))  

  w <- NULL ## matrix R x nrow(X) with tita_2 output in each column
  wSum <- NULL ## vector with sum values
  wNormalized <- NULL ## For output generation
    
  A <- NULL ## for LSE consequents estimation
      
  ##Pattern Foward Half Pass  	 
  invisible(lapply(1:nrow(object@X),function(pattern){
    ##Obtain tita_1 output as matrix (mf x input)
    tita_1 <- sapply(seq(along=object@premises), function(input, pattern){
      unlist(lapply(object@premises[[input]], function(MF){
        evaluateMF(MF, object@X[pattern, input])}))
    }, pattern)
    ##If not homogenious MFs
    if(class(tita_1)=="list"){
      aux <- matrix(NA,nrow=max(sapply(tita_1, length)), ncol=ncol(object@X))
      invisible(sapply(seq(along=tita_1), function(input){
        aux[1:length(tita_1[[input]]), input] <<- tita_1[[input]]
        return(NULL)
      }))
      tita_1 <- aux
    }
    ##MF=1 for each input
    if(!is.matrix(tita_1)){
      tita_1 <- t(tita_1)
    } 
    ##Obtain tita_2 output a.k.a. w
    w <<- cbind(w, apply(object@rules, 1, function(rule, tita_1){
      prod(sapply(seq(along=rule), function(input, tita_1){
        tita_1[rule[input],input]
      },tita_1))
    },tita_1))
    
    ## Obtain tita_3 a.k.a. w normalization 
    wSum <<- c(wSum, sum(w[, pattern]))
    wNormalized <<- cbind(wNormalized, w[, pattern]/wSum[pattern])

    ## Obtain tita_4 
    ## Generete the LSE system for the consequents
    # AX = B with each of the patterns 
    # (w1 x)p11 + (w1 y)q11 + (w1)r11 + ..... + 
    # + (w16 x)p1_16 + (w16 y)q1_16 + (w16)r1_16 = J1 for output 1
    # A is (Pattern)x(Rules*(input+1)); 
    # X is (Rules*(input+1))x(outputs) and B is (Patterns)x(outputs)
    A <<- rbind(A, matrix(matrix(c(object@X[pattern, ], 1), ncol=1) %*% 
      wNormalized[, pattern], nrow=1, byrow=FALSE))

    return(NULL)
  }))

  ##Rest of the Forward Pass
  ##Solve de LSE system by aproximations
  object@consequents <- LSE(object, A, B, initialGamma) 
  
  ##Check if object is still valid
  validObject(object) 

  ## Obtain tita_5 output estimation
  tita_5 <- A%*%object@consequents 
  object@errors <- c(object@errors, sum((object@Y-tita_5)^2))

  if(length(object@errors)!=0){
    if(abs(object@errors[length(object@errors)]) < tolerance){
      convergence <<- TRUE
    }
  }

  ##Backward Propagation
  if(!convergence){
  ##For each input
  delta_eDelta_alpha <- lapply(1:ncol(object@X), function(input){
  ##For each MembershipFunction
  lapply(1:length(object@premises[[input]]), function(MF){
    ##For each alpha, obtain delta_Error/delta_alpha
    parameters <- unlist(lapply(seq(along=object@premises[[input]][[MF]][]),
      function(alpha){
      ##Once we have defined the alpha we can calculate the 
      ##derivative_E/derivate_alpha. For each pattern once the alpha is defined 
      ##for the partial derivative of the Error.
      sum(unlist(mclapply(1:nrow(object@X), function(pattern){
        ##For each output
        sum(unlist(lapply(1:ncol(object@Y), function(output){
    #Just to avoid over calculations, 
    #lets compute the derivate delta_w*/delta_alpha != 0 once for each consequent
    rulesWithAlpha <- which(object@rules[, input, drop=FALSE] == MF)
    delta_wDelta_aplha <- sapply(rulesWithAlpha, function(wj){
      if(ncol(object@rules)>1){
        #applies the prod(mf)
        # delta_wj/delta_alpha <- delta_MF/delta_alpha * MF!=f(alpha)  
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)* 
          prod(sapply((1:ncol(object@X))[-input],function(inPut){
            evaluateMF(object@premises[[inPut]][[object@rules[wj,inPut]]], 
            object@X[pattern,input])
          }))
      }
      else{
        #for the only mf in the rule
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)
      }
    })

    ##Only for the consequents that participate in the output. Notice that 
    ##consequent are stored as object@consequents[,ouput] and not as vector for 
    ##each node
    sum(unlist(lapply(1:nrow(object@rules), function(consequent,
      rulesWithAlpha,delta_wDelta_aplha){
      # f <- c(tita_0, 1) %*% consequent for each node
      fConsequent <- as.numeric(c(object@X[pattern,],1) %*% 
        object@consequents[1:(ncol(object@X)+1) + 
        (consequent-1)*(ncol(object@X)+1),output])

      #The chain rule for the derivates of wi
      #(delta_wj/delta_alpha * wSum[pattern] - 
      #   w[consequent,pattern] * sum(delta_wj/delta_alpha)) / wSum[pattern]^2

      acum <- 0
      #Check if wj=f(alpha)
      if(any(rulesWithAlpha %in% consequent)){
        # delta_wj/delta_alpha * wSum[pattern]
        acum<-delta_wDelta_aplha[rulesWithAlpha %in% consequent]*wSum[pattern]
      }
      #w[consequent,pattern] * sum(delta_wj/delta_alpha)
      acum <- acum - w[consequent,pattern] * sum(delta_wDelta_aplha)
      acum <- acum / wSum[pattern]^2
      
      return(fConsequent * acum)
    }, rulesWithAlpha, delta_wDelta_aplha))) * 
    #For the consequents that participate in the output
    # multiplied by the first part of delta_e/delta_alpha
    (object@Y[pattern,output]-tita_5[pattern,output])*(-2)
        })))#For each output, the sum of the outputs
      })))#For each pattern, the sum of the patterns
    }))#For each alpha, obtain the vector of parameters of the corresponding MF
    names(parameters) <- names(object@premises[[input]][[MF]][])
    return(parameters)
  })#For each MembershipFunction
      })#For each input
      } #If(!convergence)
      
      ##Once we have derivative_e/derivative_alpha we can calculate 
      ##  delta_alpha = - eta * derivative_e/derivate_alpha 
      ##eta <- k / sqrt((derivative_e/derivate_alpha)^2)
         
      ## Update k
      # Rule 1 increase 10% if error undergoes four consecutive times
      if(length(object@errors)>=4){
        comparison <- cbind(last=4:2,first=3:1)
        comparison <- length(object@errors) + comparison - 4
        if(all(apply(comparison, 1, function(index){
          object@errors[index[1]]< object@errors[index[2]]}))){
            k <- 1.1 *  k
        }
      } 
      # Rule 2 decrease 10% if error undergoes two consecutive combinations of 
      # one increase and one reduction
      if(length(object@errors)>=5){
        comparison <- cbind(last=c(5,3,3,1),first=c(4,4,2,2))
        comparison <- length(object@errors) + comparison - 5  
        if(all(apply(comparison,1,function(index){
          object@errors[index[1]]< object@errors[index[2]]}))){
            k <- 0.9 * k  
        }
      } 
      
      ##Obtain eta
      eta <- k / sqrt(sum(unlist(delta_eDelta_alpha))^2)
      if(eta==Inf){
        eta <- k
      }

      ##Finally update alpha, delta_alpha = - eta * delta_e/delta_alpha 
      delta_alpha <- lapply(delta_eDelta_alpha, function(input){
        lapply(input,function(MF,input){
          return(- eta * MF[])
        },input)#For each MF
      })#For each input

      ##At last update the premises
      invisible(lapply(seq(along=object@premises), function(input){
        invisible(lapply(seq(along=object@premises[[input]]), function(MF){
          object@premises[[input]][[MF]][] <<- object@premises[[input]][[MF]][] + 
            delta_alpha[[input]][[MF]]
          return(NULL)
        }))#For each MembershipFunction
      }))#For each input
      
      epoch <- epoch + 1
      }#For each epoch

      ##Update ANFIS object
      object@fitted.values <- predict(object, object@X)
      object@residuals <- object@Y - object@fitted.values
      assign(nameObject, object, envir=parent.frame())

      ##training results
      return(list(error=object@errors, convergence=convergence, k=k))
})
#'
#' @exportMethod trainHybridOffLine
#' @docType methods
#' @name trainHybridOffLine
#' @rdname ANFIS-training
#' @aliases trainHybridOffLine-methods
setGeneric(name="trainHybridOffLine", def=function(object, epochs=5, 
  tolerance=1e-5, initialGamma=1000, eta=0.05, phi=0.2, a=0.01, b=0.1, 
  delta_alpha_t_1=list()){
    standardGeneric("trainHybridOffLine")
})
#'
#' @name trainHybridOffLine
#' @rdname ANFIS-training
#' @aliases trainHybridOffLine,ANFIS-method
#' @inheritParams LSE
setMethod(f="trainHybridOffLine", signature="ANFIS", 
  definition=function(object, epochs=5, tolerance=1e-5, initialGamma=1000, 
  eta=0.05, phi=0.2, a=0.01, b=0.1, delta_alpha_t_1=list()){
  ##Initialization of variables
  nameObject <- deparse(substitute(object))
  object@call <- match.call()
  object@trainingType <- "trainHybridOffLine"
  convergence <- FALSE
  B <- object@Y ## for LSE consequents estimation
  # for momentum term
  if(length(delta_alpha_t_1)==0)
  delta_alpha_t_1 <- lapply(seq(along=object@premises), function(input){
    lapply(seq(along=object@premises[[input]]), function(MF){
      return(0*object@premises[[input]][[MF]][])
    })#For each MembershipFunction  
  })#For each input
      
  ##For each epoch
  epoch <- 1
  while(epoch<epochs & !convergence){
  print(paste("epoch: ", epoch))  
  
  w <- NULL ## matrix R x nrow(X) with tita_2 output in each column
  wSum <- NULL ## vector with sum values
  wNormalized <- NULL ## For output generation
    
  A <- NULL ## for LSE consequents estimation
      
  ##Pattern Foward Half Pass
  invisible(lapply(1:nrow(object@X), function(pattern){
    ##Obtain tita_1 output
    tita_1 <- sapply(seq(along=object@premises), function(input, pattern){
      unlist(lapply(object@premises[[input]], function(MF){
        evaluateMF(MF, object@X[pattern, input])}))
    },pattern)
    ##If not homogenious MFs
    if(class(tita_1)=="list"){
      aux <- matrix(NA, nrow=max(sapply(tita_1, length)), ncol=ncol(object@X))
      invisible(sapply(seq(along=tita_1), function(input){
        aux[1:length(tita_1[[input]]), input] <<- tita_1[[input]]
        return(NULL)
      }))
      tita_1 <- aux
    }
    ##MF=1 for each input
    if(!is.matrix(tita_1)){
      tita_1 <- t(tita_1)
    } 
    ##Obtain tita_2 output a.k.a. w
    w <<- cbind(w, apply(object@rules, 1, function(rule,tita_1){
      prod(sapply(seq(along=rule), function(input, tita_1){
        tita_1[rule[input],input]},tita_1))
    },tita_1))

    ## Obtain tita_3 a.k.a. w normalization 
    wSum <<- c(wSum, sum(w[, pattern]))
    wNormalized <<- cbind(wNormalized, w[, pattern]/wSum[pattern])

    ## Obtain tita_4 
    ## Generete the LSE system for the consequents
    # AX = B with each of the patterns 
    # (w1 x)p11 + (w1 y)q11 + (w1)r11 + ..... + 
    #  + (w16 x)p1_16 + (w16 y)q1_16 + (w16)r1_16 = J1 for output 1
    # A is (Pattern)x(Rules*(input+1)); 
    # X is (Rules*(input+1))x(outputs) and B is (Patterns)x(outputs)
    A <<- rbind(A, matrix(matrix(c(object@X[pattern, ], 1), ncol=1) %*% 
      wNormalized[, pattern], nrow=1, byrow=FALSE))

    return(NULL)
  }))

      ##Rest of the Forward Pass
  ##Solve de LSE system by aproximations
  object@consequents <- LSE(object, A, B, initialGamma) 
  ##Check if object is still valid
  validObject(object) 

  ## Obtain tita_5 output estimation
  tita_5 <- A%*%object@consequents 
  object@errors <- c(object@errors, sum((object@Y-tita_5)^2))

  if(length(object@errors)!=0){
    if(abs(object@errors[length(object@errors)]) < tolerance){
      convergence <<- TRUE
    }
  }

  ##Backward Propagation
  if(!convergence){
    ##For each input
    delta_eDelta_alpha <- lapply(1:ncol(object@X), function(input){
  ##For each MembershipFunction
  lapply(1:length(object@premises[[input]]), function(MF){
    ##For each alpha, obtain delta_Error/delta_alpha
    parameters <- unlist(lapply(seq(along=object@premises[[input]][[MF]][]),
      function(alpha){
      ##Once we have defined the alpha we can calculate the 
      ## derivative_E/derivate_alpha. For each pattern once the alpha is defined
      ##for the partial derivative of the Error
      sum(unlist(mclapply(1:nrow(object@X), function(pattern){
        ##For each output
        sum(unlist(lapply(1:ncol(object@Y), function(output){
    #Just to avoid over calculations, lets compute the 
    #derivate delta_w*/delta_alpha != 0 once for each consequent
    rulesWithAlpha <- which(object@rules[,input,drop=FALSE] == MF)
    delta_wDelta_aplha <- sapply(rulesWithAlpha, function(wj){
      if(ncol(object@rules)>1){
        #applies the prod(mf)
        # delta_wj/delta_alpha <- delta_MF/delta_alpha * MF!=f(alpha)  
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)* 
          prod(sapply((1:ncol(object@X))[-input], function(inPut){
            evaluateMF(object@premises[[inPut]][[object@rules[wj, inPut]]], 
            object@X[pattern,input])
          }))
      }else{
        #for the only mf in the rule
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)
      }
    })


    ##Only for the consequents that participate in the output. Notice that 
    ##consequent are stored as object@consequents[,ouput] and not as vector 
    ##for each node
    sum(unlist(lapply(1:nrow(object@rules), function(consequent, rulesWithAlpha,
      delta_wDelta_aplha){
      # f <- c(tita_0, 1) %*% consequent for each node
      fConsequent <- as.numeric(c(object@X[pattern,],1) %*% 
        object@consequents[1:(ncol(object@X)+1) + 
          (consequent-1)*(ncol(object@X)+1),output])

      #The chain rule for the derivates of wi
      #(delta_wj/delta_alpha * wSum[pattern] - 
      #     w[consequent,pattern] * sum(delta_wj/delta_alpha))/wSum[pattern]^2

      acum <- 0
      #Check if wj=f(alpha)
      if(any(rulesWithAlpha %in% consequent)){
        # delta_wj/delta_alpha * wSum[pattern]
        acum <- delta_wDelta_aplha[rulesWithAlpha %in% consequent]*wSum[pattern]
      }
      #w[consequent,pattern] * sum(delta_wj/delta_alpha)
      acum <- acum - w[consequent,pattern] * sum(delta_wDelta_aplha)
      acum <- acum / wSum[pattern]^2

      return(fConsequent * acum)
    }, rulesWithAlpha, delta_wDelta_aplha))) * 
    #For the consequents that participate in the output
    # multiplied by the first part of delta_e/delta_alpha
      (object@Y[pattern,output]-tita_5[pattern,output])*(-2)
        })))#For each output, the sum of the outputs
      })))#For each pattern, the sum of the patterns
    }))#For each alpha, obtain the vector of parameters of the corresponding MF
    names(parameters) <- names(object@premises[[input]][[MF]][])
    return(parameters)
  })#For each MembershipFunction
      })#For each input
      } #If(!convergence)
      
      ##Once we have derivative_e/derivative_alpha we can calculate delta_alpha 
      ##delta_alpha_t+1= - eta * derivative_e/derivate_alpha + phi delta_alpha_t
      
      delta_alpha <- lapply(seq(along=delta_eDelta_alpha), function(input){
        lapply(seq(along=delta_eDelta_alpha[[input]]),function(MF,input){
          return(- eta * delta_eDelta_alpha[[input]][[MF]][] + 
            phi * delta_alpha_t_1[[input]][[MF]][])
        },input)#For each MF
      })#For each input
  
      ##At last update the premises
      invisible(lapply(seq(along=object@premises), function(input){
        invisible(lapply(seq(along=object@premises[[input]]), function(MF){
          object@premises[[input]][[MF]][] <<- object@premises[[input]][[MF]][] + 
            delta_alpha[[input]][[MF]]
          return(NULL)
        }))#For each MembershipFunction
      }))#For each input
      
      ##Update eta 
      if(length(object@errors) >= 2){
        if(object@errors[length(object@errors)]-
          object@errors[length(object@errors)-1] < 0){
          eta <- eta + a
        }else{
          a <- eta * b
          eta <- eta * (1-b)
        }
      }

      ##Update delta_alpha_t_1
      delta_alpha_t_1 <- delta_alpha 

      epoch <- epoch + 1
      }#For each epoch

      ##Update ANFIS object
      object@fitted.values <- predict(object, object@X)
      object@residuals <- object@Y - object@fitted.values
      assign(nameObject,object,envir=parent.frame())

      ##training results
      return(list(error=object@errors, convergence=convergence, eta=eta, a=a, 
        delta_alpha_t_1=delta_alpha_t_1))
})
#'
#' @exportMethod trainHybridJangOnLine
#' @docType methods
#' @name trainHybridJangOnLine
#' @rdname ANFIS-training
#' @aliases trainHybridJangOnLine-methods
setGeneric(name="trainHybridJangOnLine", def=function(object, epochs=5, 
  tolerance=1e-15, initialGamma=1000, k=0.01, lamda=0.9, 
  S=matrix(nrow=0,ncol=0)){
    standardGeneric("trainHybridJangOnLine")
})
#'
#' @name trainHybridJangOnLine
#' @rdname ANFIS-training
#' @aliases trainHybridJangOnLine,ANFIS-method
#' @inheritParams LSE
setMethod(f="trainHybridJangOnLine", signature="ANFIS", 
  definition=function(object, epochs=5, tolerance=1e-15, initialGamma=1000, 
  k=0.01, lamda=0.9, S=matrix(nrow=0,ncol=0)){
  ##Initialization of variables
  nameObject <- deparse(substitute(object))
  object@call <- match.call()
  object@trainingType <- "trainHybridJangOnLine"
  convergence <- FALSE
  if(length(S)==0){
    ##Not Yet initialize for LSE
    S<-initialGamma*diag(nrow(object@rules)*(ncol(object@X)+1))
  }
  
  ##For each epoch and pattern
  epoch <- 1
  pattern <- 1
  print(paste("epoch: ", epoch))
  while((epoch < epochs) & (pattern <= nrow(object@X)) & !convergence){
  
  ##Pattern Foward Half Pass
    ##Obtain tita_1 output as matrix (mf x input)
    tita_1 <- sapply(seq(along=object@premises), function(input, pattern){
      unlist(lapply(object@premises[[input]], function(MF){
        evaluateMF(MF,object@X[pattern,input])
      }))
    },pattern)
    ##If not homogenious MFs
    if(class(tita_1)=="list"){
      aux <- matrix(NA,nrow=max(sapply(tita_1, length)), ncol=ncol(object@X))
      invisible(sapply(seq(along=tita_1), function(input){
        aux[1:length(tita_1[[input]]),input] <<- tita_1[[input]]
        return(NULL)
      }))
      tita_1 <- aux
    }
    ##MF=1 for each input
    if(!is.matrix(tita_1)){
      tita_1 <- t(tita_1)
    } 
    ##Obtain tita_2 output a.k.a. w
    w <- matrix(apply(object@rules,1,function(rule,tita_1){
      prod(sapply(seq(along=rule),function(input,tita_1){
        tita_1[rule[input],input]
      },tita_1))
    },tita_1),ncol=1)
    
    ## Obtain tita_3 a.k.a. w normalization 
    wSum <- sum(w)
    wNormalized <- w/wSum

    ## Obtain tita_4 
    ## Generete the LSE system for the consequents
    # AX = B with each of the patterns 
    # (w1 x)p11 + (w1 y)q11 + (w1 1)r11 + ..... 
    #   + (w16 x)p1_16 + (w16 y)q1_16 + (w16)r1_16 = J1 for output 1
    # A is (Pattern)x(Rules*(input+1)); 
    # X is (Rules*(input+1))x(outputs) and B is (Patterns)x(outputs)
    a <- t(matrix(matrix(c(object@X[pattern, ], 1), ncol=1) %*% 
      t(wNormalized), nrow=1, byrow=FALSE))
    b <- t(object@Y[pattern, , drop=FALSE])

    ##Rest of the Forward Pass
    ##Solve de LSE system by aproximations
    ## Recurrence formula for on-line learning
    S<-as.numeric(1/lamda)*(S-(S%*%a%*%t(a)%*%S)/as.numeric(lamda+t(a)%*%S%*%a))
    object@consequents <- object@consequents + 
      S%*%a%*%(t(b)-t(a)%*%object@consequents)
    
    ##Check if object is still valid
    validObject(object) 
  
    ## Obtain tita_5 output estimation
    tita_5 <- t(a)%*%object@consequents 
    object@errors <- c(object@errors, 
      sum((object@Y[pattern, , drop=FALSE]-tita_5)^2))

    ##Check convergence
    if(length(object@errors)!=0){
      if(abs(object@errors[length(object@errors)]) < tolerance){
        convergence <- TRUE
      }
    }

    ##Backward Propagation
    if(!convergence){
    ##For each input
    delta_eDelta_alpha <- mclapply(1:ncol(object@X), function(input, pattern){
      ##For each MembershipFunction
      lapply(1:length(object@premises[[input]]), function(MF){
        ##For each alpha, obtain delta_Error/delta_alpha
        parameters <- unlist(lapply(seq(along=object@premises[[input]][[MF]][]),
        function(alpha){
          ##Once we have defined the alpha we can calculate the 
          ##derivative_E/derivate_alpha. Once the alpha is definned for the 
          ##partial derivative of the Error
          ##For each output
          sum(unlist(lapply(1:ncol(object@Y),function(output){
    #Just to avoid over calculations, lets compute the 
    #derivate delta_w*/delta_alpha != 0 once for each consequent
    rulesWithAlpha <- which(object@rules[,input,drop=FALSE] == MF)
    delta_wDelta_aplha <- sapply(rulesWithAlpha, function(wj){
      if(ncol(object@rules)>1){
        #applies the prod(mf)
        # delta_wj/delta_alpha <- delta_MF/delta_alpha * MF!=f(alpha)  
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)* 
          prod(sapply((1:ncol(object@X))[-input], function(inPut){
            evaluateMF(object@premises[[inPut]][[object@rules[wj, inPut]]], 
            object@X[pattern,input])
          }))
      }else{
        #for the only mf in the rule
        derivateMF(object@premises[[input]][[MF]],object@X[pattern,input],alpha)
      }
    })

    ##Only for the consequents that participate in the output. Notice that 
    ##consequent are stored as object@consequents[,ouput] and not as vector 
    ##for each node
    sum(unlist(lapply(1:nrow(object@rules), function(consequent, rulesWithAlpha,
      delta_wDelta_aplha){
      # f <- c(tita_0, 1) %*% consequent for each node
      fConsequent <- as.numeric(c(object@X[pattern,],1) %*% 
        object@consequents[1:(ncol(object@X)+1) + 
          (consequent-1)*(ncol(object@X)+1),output])
    
      #The chain rule for the derivates of wi
      #(delta_wj/delta_alpha * wSum[pattern] - 
      #   w[consequent,pattern] * sum(delta_wj/delta_alpha))/wSum[pattern]^2

      acum <- 0
      #Check if wj=f(alpha)
      if(any(rulesWithAlpha %in% consequent)){
          # delta_wj/delta_alpha * wSum[pattern]
          acum <- delta_wDelta_aplha[rulesWithAlpha %in% consequent] * wSum
      }
      #w[consequent,pattern] * sum(delta_wj/delta_alpha)
      acum <- acum - w[consequent] * sum(delta_wDelta_aplha)
      acum <- acum / wSum^2

      return(fConsequent * acum)
    }, rulesWithAlpha, delta_wDelta_aplha))) * 
    #For the consequents that participate in the output
    # multiplied by the first part of delta_e/delta_alpha
      (object@Y[pattern,output]-tita_5[output])*(-2)
        })))#For each output, the sum of the outputs
    }))#For each alpha, obtain the vector of parameters of the corresponding MF
    names(parameters) <- names(object@premises[[input]][[MF]][])
    return(parameters)
  })#For each MembershipFunction
      },pattern)#For each input
      
      ##Once we have derivative_e/derivative_alpha we can calculate 
      ##delta_alpha = - eta * derivative_e/derivate_alpha 
      ##eta <- k / sqrt((derivative_e/derivate_alpha)^2)
         
      ## Update k
      # Rule 1 increase 10% if error undergoes four consecutive times
      if(length(object@errors)>=4){
        comparison <- cbind(last=4:2,first=3:1)
        comparison <- length(object@errors) + comparison - 4
        if(all(apply(comparison, 1, function(index){
          object@errors[index[1]]< object@errors[index[2]]}))){
            k <- 1.1 * k
        }
      } 
      # Rule 2 decrease 10% if error undergoes two consecutive combinations of 
      #one increase and one reduction
      if(length(object@errors)>=5){
        comparison <- cbind(last=c(5,3,3,1),first=c(4,4,2,2))
        comparison <- length(object@errors) + comparison - 5  
        if(all(apply(comparison, 1, function(index){
          object@errors[index[1]]< object@errors[index[2]]}))){
          k <- 0.9 * k  
        }
      } 
      
      ##Obtain eta
      eta <- k / sqrt(sum(unlist(delta_eDelta_alpha))^2)
      if(eta==Inf){eta <- k}

      ##Finally update alpha, delta_alpha = - eta * delta_e/delta_alpha 
      delta_alpha <- lapply(delta_eDelta_alpha, function(input){
        lapply(input,function(MF,input){
          return(- eta * MF[])
        },input)#For each MF
      })#For each input

      ##At last update the premises
      invisible(lapply(seq(along=object@premises), function(input){
        invisible(lapply(seq(along=object@premises[[input]]), function(MF){
          object@premises[[input]][[MF]][] <<- object@premises[[input]][[MF]][] + 
            delta_alpha[[input]][[MF]]
          return(NULL)
        }))#For each MembershipFunction
      }))#For each input

      }#!convergence

    
    ##Increase pattern and epoch
    pattern <- (pattern + 1)%%(nrow(object@X)+1)
    if(pattern==0){
      pattern <- 1
      epoch <- epoch + 1
      if(epoch < epochs) {
        print(paste("epoch: ", epoch))
      }
    }

    }#For each epoch and pattern


    ##Update ANFIS object
    object@fitted.values <- predict(object, object@X)
    object@residuals <- object@Y - object@fitted.values
    assign(nameObject,object,envir=parent.frame())

    ##training results
    return(list(error=object@errors, convergence=convergence, k=k, S=S))
})
