###########################################################################
#  Original program realized by Hélène Crépeau (21-02-2011)
#  Updates realized by Sophie Baillargeon
############################################################################

petersen=function(mat, M, S=0)
{ 
  ### Validation des arguments
  mat <- as.matrix(mat)
  if (dim(mat)[2]!=2) stop("'mat' must have two columns")
#  if (any((mat%%1)!=0)||any(mat<0)) stop("'mat' must contain non-negative integers")
  if (!is.numeric(mat)||any(mat<0)) stop("'mat' must contain non-negative numerics")
  if (any(mat[,2]<mat[,1])) stop("each number in the second column of 'mat' must be greater or equal to the number in the first column on the same row") 
  
  valid.one(M,"numeric")
  
  valid.one(S,"numeric")
  ############
  
# gni : vector of the size of detected  groups 
  gni=mat[,2]
  
  mat_aggre=mat[gni>S,]
  G=dim(mat_aggre)[1]
  R=sum(mat_aggre[,1])
  C=sum(mat_aggre[,2])
  
  if(M<R) warning("")
  
  LP_T.hat=((M+1)*(C+1))/(R+1)-1
  var_LP=((M+1)*(C+1)*(M-R)*(C-R))/((R+1)^2*(R+2))
  se_LP=sqrt(var_LP)
  
  out<-list(G=G,R=R,C=C,LP_T.hat=LP_T.hat,se_LP=se_LP,mat_aggre=mat_aggre,call=match.call())
  class(out)<- "petersen"    
  return(out)
  
}


"print.petersen" <- function(x, ...) {
  cat("Given parameters:")
  cat("\nM: Total number of active collars during the census =",eval(x$call$M))
  S <- eval(x$call$S)
  if (is.null(S)) S <- 0  # If S was not given, it does not appear in the function call, but it takes the default value 0.
  cat("\nS: Minimum size of well aggregated groups =",S) 
  
  cat("\n\nDescriptive statistics on well aggregated groups:")
  cat("\nG: Number of well aggregated groups =",x$G) 
  cat("\nR: Total number of radio-collared animal observed in the well aggregated groups =",x$R)
  cat("\nC: Total number of animals observed in the well aggregated groups =",x$C)
  
  cat("\n\nEstimation:\n")
  tableau3 <- cbind(estimate=round(x$LP_T.hat,0),se=round(x$se_LP,1))
  rownames(tableau3) <- c("Total number of animals in a herd:")
  print.default(tableau3, print.gap = 2, quote = FALSE, right=TRUE, ...)
  
  cat("\nThis function calculates a Lincoln-Petersen estimator with Chapman's bias correction and the bias corrected standard error estimator of Seber and Wittes\n")
  
  invisible(x)
}



