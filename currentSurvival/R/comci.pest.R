
comci.pest <- function(t, LastContact, Exitus, maxx)
{

NoPat <- length(t) # a number of patients

# replace NAs in t with maximum follow-up and compute a vector 'ind' with failure causes:
ind <- array(1,c(NoPat,1)) # allocation of a vector with failure causes
for (i in 1:NoPat){ # do for each patient:
  if (is.na(t[i])){t[i] <- LastContact[i] # replace NA in t with maximum follow-up
    if (!(Exitus[i])){ind[i] <- 0} # if there is NA in t and the patient is not dead, the exact ith failure time is not observed
    else {ind[i] <- 2} # if there is NA in t and patient is dead, the exact failure time is known and the failure cause is death
  }
}
# ind - a vector with failure causes:
# 0..the exact ith failure time is not observed
# 1..the exact failure time is observed and the failure cause is the achievement of the first disease remission
# 2..the exact failure time is observed and the failure cause is death

# estimate common cumulative incidence function:
if (sum(ind==1)>0) { # check whether there is at least one patient with the achievement of the first disease remission
  pomRi <- cuminc(t,ind) # create a list with estimates of cumulative incidence function (est) and variance (var); the achievement of the first disease remission is an event, death is a competitive risk 
  pomRi2 <- timepoints(pomRi,seq(0,maxx,1)) # a list with estimates of cumulative incidence function (est) and variance (var) in each day
  R1 <- pomRi2$est[1,] # the estimates of cumulative incidence function in each day
  R1[is.na(R1)] <- R1[sum(!is.na(R1))] # replace NAs with the last estimate before NAs
} else { # if there is no patient with the achievement of the first disease remission, cumulative incidence function is equal to a null curve
  R1 <- array(0,maxx+1)
}

comci.pest <- list(x=0:(maxx),y=R1)

}