
lfs.pest <- function(t, LastContact, Exitus, maxx)
{

NoPat <- length(t) # a number of patients

# replace NAs in t with maximum follow-up and compute a vector 'ind' with failure causes:
ind <- array(1,c(NoPat,1)) # allocation of a vector with failure causes
for (i in 1:NoPat){ # do for each patient:
  if (is.na(t[i])){t[i] <- LastContact[i] # replace NA in t with maximum follow-up
  if (!(Exitus[i])){ind[i] <- 0} }} # if there is NA in t and the patient is not dead, the exact ith failure time is not observed
# ind - a vector with failure causes:
# 0..the exact ith failure time is not observed
# 1..the exact failure time is observed and the failure cause is the loss of disease remission or the death
  
# estimate common leukaemia-free survival function:
S1star = summary(survfit(Surv(t,ind)~1)) # a survival function estimate corresponding to the loss of the first disease remission
S1_star_days = stretch(list(x=S1star$time,y=S1star$surv),max(LastContact)) # assignment of survival estimates to each day of the follow-up

lfs.pest <- list(x=0:maxx,y=S1_star_days$y[1:(maxx+1)])

}
