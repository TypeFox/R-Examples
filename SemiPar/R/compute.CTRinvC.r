########## compute.CTRinvC.S #######################################

# Computes the inverse residual covariance matrix R
# for linear mixed models

# inputs: individual = id's for the "individiduals" on whom
#         the repeated mesaures are taken. This vector has
#         the same length as the response and is assumed to
#         be sorted into homogenous groupings. i.e.: 
#         1 1 1 2 2 3 3 3 3 3 is OK
#         1 2 3 1 2 3 1 3 3 3 is not
#         rho = AR(1) parameter
#         
#
# output: t(C) R.unscaled.inv C
#         C = [X Z]
#
# Last changed: 10/20/00


compute.CTRinvC <- function(X,Z,individual,rho=0)
{
  C <- cbind(X,Z)
  CTRinvC <- 0*t(C)%*%C

  # Calculate indices (matrix location) and time indices for each individiual

  if (length(individual)>1)
  {
    end.indiv <- (1:length(individual))[diff(individual)!=0]
    strt.indiv <- c(1,end.indiv+1)
    end.indiv <- c(end.indiv,length(individual))
    n.indiv <- length(end.indiv)
    indiv.inds <- NULL
    indiv.times <- NULL
    
    for (i in 1:n.indiv)
    {
       indiv.inds <- (strt.indiv[i]:end.indiv[i])
       indiv.times <- (1:length(indiv.inds))
       R.inv <- 
            solve(rho^abs(outer(indiv.times,indiv.times,"-")))
       if (ncol(R.inv)>1)
       {
          CTRinvC <- CTRinvC +
          t(C[indiv.inds,]) %*% R.inv %*% C[indiv.inds,]
       }
       else
       {
          CTRinvC <- CTRinvC +
               as.numeric(R.inv)*outer(C[indiv.inds,],C[indiv.inds,])
                  
       }
     }
  }
  else
  {
      CTRinvC <- t(C)%*%C
  }
    
  return(CTRinvC)
}

########## end compute.CTRinvC.S #######################################
