
cci.pest <- function(E, LastContact, Exitus, maxx)
{

NoPat <- length(E[,1]) # a number of patients
r <- ceiling(length(E[1,])/2) # a maximum number of disease remissions achieved by patients

# check whether a number of columns in the matrix E is even:
if ((length(E[1,])%%2)>0) E <- cbind(E,array(NA,c(NoPat,1))) # if the number of columns in the matrix E is not even, add a vector with NAs
  
# replace NAs in the matrix E with maximum follow-up and create a matrix R with failure causes:
R <- array(1,c(NoPat,2*r)) # allocation of a matrix with failure causes
for (i in 1:NoPat){ # do for each patient:
  for (j in 1:(2*r)){ # do for each event:
    if (is.na(E[i,j])){E[i,j] <- LastContact[i] # replace NA in E with maximum follow-up
      if (!(Exitus[i])){R[i,j] <- 0} # if there is NA in E and the patient is not dead, the exact ith failure time is not observed
      else {R[i,j] <- 2} }}} # if there is NA in E and the patient is dead, the exact failure time is known and the failure cause is death
# R - a matrix with failure causes:
# 0..the exact ith failure time is not observed
# 1..the exact failure time is observed and the failure cause is the achievement or the loss of disease remission
# 2..the exact failure time is observed and the failure cause is death
  
# allocate vectors with sums of cumulative incidence functions:
SumI <- array(0,maxx+1) # allocation of a vector with sum of cumulative incidence functions corresponding to achievements of disease remissions
SumI_star <- array(0,maxx+1) # allocation of a vector with sum of cumulative incidence functions corresponding to losses of disease remissions
SumI_death <- array(0,maxx+1) # allocation of a vector with sum of cumulative incidence functions corresponding to death in the 1st, 2nd etc. disease remission

# compute the sum of cumulative incidence functions corresponding to achievements of disease remissions:
for (i in seq(1,2*r,2)){	
  if (sum(R[,i]==1)>0) { # test if there is at least one patient with the achievement of the ith disease remission
    pIi <- cuminc(E[,i],R[,i]) # a list with estimates of cumulative incidence function (est) and variance (var) 
    pIi2 <- timepoints(pIi,seq(0,maxx,1)) # a list with estimates of cumulative incidence function (est) and variance (var) in each day
    Ii <- pIi2$est[1,] # the estimates of cumulative incidence function in each day
    Ii[is.na(Ii)] <- Ii[sum(!is.na(Ii))] # replace NAs with the last number before NAs
  } else { # if there is no patient with the achievement of the ith disease remission, cumulative incidence function is equal to a null curve
    Ii <- array(0,maxx+1)
  }
  SumI <- SumI+Ii
}

# compute the sum of cumulative incidence functions corresponding to losses of disease remissions:
for (i in seq(2,2*r,2)){
  if (sum(R[,i]==1)>0) { # test if there is at least one patient with the loss of the ith disease remission
    pIi_star <- cuminc(E[,i],R[,i]) # a list with estimates of cumulative incidence function (est) and variance (var) 
    pIi2_star <- timepoints(pIi_star,seq(0,maxx,1)) # a list with estimates of cumulative incidence function (est) and variance (var) in each day
    Ii_star <- pIi2_star$est[1,] # the estimates of cumulative incidence function in each day
    Ii_star[is.na(Ii_star)] <- Ii_star[sum(!is.na(Ii_star))] # replace NAs with the last number before NAs
  } else { # if there is no patient with the loss of the ith disease remission, cumulative incidence function is equal to a null curve
    Ii_star <- array(0,maxx+1)
  }
  SumI_star <- SumI_star+Ii_star
}

# compute the sum of cumulative incidence functions corresponding to death in the 1st, 2nd etc. disease remissions:
for (i in seq(1,2*r,2)){	
  if (sum(R[,i]==1)>0) { # test if there is at least one patient with the achievement of the ith disease remission
    Ri_death <- array(0,nrow(R)) # allocation of a vector with failure causes
    Ri_death[R[,i+1]==2] <- 1 # the event of interest (death in ith disease remission) (next command will replace deaths before ith disease remission with the value 2)
    Ri_death[R[,i]==2] <- 2 # death before the achievement of ith disease remission is a competitive risk
    Ri_death[R[,i+1]==1] <- 2 # ith loss of disease remission is a competitive risk
    if (sum(Ri_death==1)>0) { # test if there is at least one patient with death in ith disease remission
      Ei_death <- E[,i+1] # a vector with times to events
      pIi_death <- cuminc(Ei_death,Ri_death) # a list with estimates of cumulative incidence function (est) and variance (var) 
      pIi2_death <- timepoints(pIi_death,seq(0,maxx,1)) # a list with estimates of cumulative incidence function (est) and variance (var) in each day
      Ii_death <- pIi2_death$est[1,] # the estimates of cumulative incidence function in each day
      Ii_death[is.na(Ii_death)] <- Ii_death[sum(!is.na(Ii_death))] # replace NAs with the last number before NAs
    } else { # if there is no patient with death in ith disease remission, cumulative incidence function of deaths in ith disease remission is equal to a null curve
      Ii_death <- array(0,maxx+1)
    }
  } else { # if there is no patient with the achievement of the ith disease remission, cumulative incidence function of deaths in ith disease remission is equal to a null curve
      Ii_death <- array(0,maxx+1)
  }
  SumI_death <- SumI_death+Ii_death
}

cci.pest <- list(x=0:(maxx),y=SumI-SumI_star-SumI_death) # estimation of the current cumulative incidence function as SumI-SumI_star-SumI_death

}
