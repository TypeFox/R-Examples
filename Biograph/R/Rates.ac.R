Rates.ac <-
function(Stable) 
{
  #  Stable <- trans$Stable
  # Calculate transition rates (occurrence-exposure rates): rates
  #  Input: Stable.out, created in TAB_TRANS.r
  nage <- nrow(Stable)-2
  namstates <- unlist (unname(dimnames(Stable))[3])
  numstates <- length (namstates)
  namage <- unlist (unname(dimnames(Stable))[1])[1:nage]
  
  rates <- Stable[1:nage,4:(3+numstates),] # rates: origin = third dimension
  for (k in 1:numstates) rates[,k,] <- rates[,k,]/Stable[1:nage,2,]
  for (ix in 1:nage) diag(rates[ix,,]) <- 0
  dimnames(rates) <- list(age=namage,destination=namstates,origin=namstates)
  # Create Mac matrix (period-cohort): M = age,sex, destin, origin
  # Mac <- array (0,c(nage,numstates,numstates))
  Mac <- - aperm(rates,c(1,3,2)) #  (old)
  Mac <- rates
    # Replace infinite rate (when PY =0) by a large number
   Mac[Mac==Inf | Mac > 2 ] <- 2
   Mac[is.na(Mac)] =0
  for (ix in 1:nage)
  { # diag(Mac[ix,is,,]) <- mortality rate 
     diag(Mac[ix,,]) <- - apply(Mac[ix,,],2,sum) # column sum    
      # diagonal = total outmig rate + mortality rates
     # for (i in 1:numstates) {for (j in 1:numstates) if (is.na(Mac[ix,i,j])) Mac[ix,i,j] <- 0}
  }

   
   Mcum<- apply(Mac,c(2,3),cumsum) #Lambda.oe[51,,] = includes rates of 50-51 -> for P at age 51 = P[52,,]
   rownames(Mcum) <- c(rownames(Mac)[2:(length(rownames(Mac)))],"54")

 return (list (M = Mac,
               Mcum = Mcum))# origin = third dimension See Mac[,,1:2]
 }
