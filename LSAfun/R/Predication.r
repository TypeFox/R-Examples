##### Predication ##############################


# Print Method

print.Pred <- function(obj){
  cat("\nPredication Results\n\n")
  cat("Input: Predicate    ", obj$P,"    Argument:    ", obj$A, "\n\n")
  cat("Predication Vector ($PA):\n ",(obj$PA), "\n\n")
  cat("Predicate Vector without Argument ($P.Pred):\n ",(obj$P.Pred), "\n\n")
  cat("Used neighborhood words ($neighbors):\n", obj$neighbors, "\n")
  invisible(obj)
}

# Function

#' @export
Predication <- function(P,A,m,k,tvectors=tvectors,breakdown=FALSE,norm="none"){
  
  if(class(tvectors) == "matrix"){
    
    if(breakdown==TRUE){
      
      P <- breakdown(P)
      A <- breakdown(A)
      
    }
    
    tvectors_P   <- tvectors[P,]
    tvectors_A   <- tvectors[A,]
    nrow         <- nrow(tvectors)
    
    
    # m nearest types to P
    
    
    near.P <- neighbors(P,(m+1),tvectors=tvectors,
                         breakdown=breakdown)[2:(m+1)]   
    
    
    
    # k nearest to A
    
    near.PA <- multicos(A,names(near.P),tvectors=tvectors,
                        breakdown=breakdown)[1,]
    near.PA <- sort(near.PA,decreasing=T)[1:k]
    
    
    neighbors <- names(near.PA)
    tvectors_PA <- tvectors[names(near.PA),]
    
    
    
    if(k==1){
      
      if(norm=="none"){
        PA <- tvectors[P,]+tvectors[A,]+tvectors[names(near.PA),]
        P.Pred <- tvectors[P,]+tvectors[names(near.PA),]
      }
      
      if(norm=="all"){
        PA <- normalize(tvectors[P,]) + normalize(tvectors[A,]) + normalize(tvectors[names(near.PA),])
        P.Pred <- normalize(tvectors[P,]) + normalize(tvectors[names(near.PA),])
      }
      
      if(norm=="block"){
        PA <- normalize(tvectors[A,]) + normalize(tvectors[P,] + tvectors[names(near.PA),])
        P.Pred <- normalize(tvectors[P,] + tvectors[names(near.PA),])
      }
      
    }
    
    
    if(k >1){
      
      if(norm=="none"){
        PA <- tvectors[A,] + tvectors[P,] + colSums(tvectors[names(near.PA),])
        P.Pred <- tvectors[P,] + colSums(tvectors[names(near.PA),])
      }
      
      if(norm=="all"){        
        normPA <- tvectors[names(near.PA),]
        normPA <- t(apply(normPA,1,normalize))
        
        PA <- normalize(tvectors[A,]) + normalize(tvectors[P,]) + colSums(normPA)
        P.Pred <- normalize(tvectors[P,]) + colSums(normPA)
      }
      
      if(norm=="block"){
        PA <- normalize(tvectors[A,]) + normalize(tvectors[P,] + colSums(tvectors[names(near.PA),]))
        P.Pred <- normalize(tvectors[P,] + colSums(tvectors[names(near.PA),]))
      }
      
      
    }
    
    
    
    out <- list(PA=PA,P.Pred=P.Pred,neighbors=neighbors,P=P,A=A)
    
    class(out) <- "Pred"
    out
    
  }else{
    stop("tvectors must be a matrix!")
  } 
  
}
