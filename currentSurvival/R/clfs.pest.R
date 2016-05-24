
clfs.pest <- function(E, LastContact, Exitus, maxx)
{

NoPat <- length(E[,1]) # a number of patients
r <- floor(length(E[1,])/2)+1 # a maximum number of disease remissions achieved by patients

# check whether a number of columns in the matrix E is odd:
if ((length(E[1,])%%2)==0) E <- cbind(E,array(NA,c(NoPat,1))) # if the number of columns in the matrix E is not odd, add a vector with NAs
  
# replace NAs in the matrix E with maximum follow-up and create a matrix R with failure causes:
R <- array(1,c(NoPat,2*r-1)) # allocation of a matrix with failure causes
for (i in 1:NoPat){ # do for each patient:
  for (j in 1:(2*r-1)){ # do for each event:
    if (is.na(E[i,j])){E[i,j] <- LastContact[i] # replace NA in E with maximum follow-up
      if (!(Exitus[i])){R[i,j] <- 0}}}} # if there is NA in E and the patient is not dead, the exact ith failure time is not observed
# R - a matrix with failure causes:
# 0..the exact ith failure time is not observed
# 1..the exact failure time is observed and the failure cause is the achievement or the loss of disease remission or death

# allocate vectors with sums of survival functions:
SumS = array(0,maxx+1) # allocation of a vector with sum of survival functions corresponding to achievements of disease remissions
SumS_star = array(0,maxx+1) # allocation of a vector with sum of survival functions corresponding to losses of disease remissions

# estimate survival function corresponding to the loss of the first disease remission:
S1star = summary(survfit(Surv(E[,1],R[,1])~1)) # a survival estimate corresponding to the loss of the first disease remission
S1_star_days = stretch(list(x=S1star$time,y=S1star$surv),max(LastContact)) # assignment of survival estimates to each day of the follow-up
S1_star = S1_star_days$y[1:(maxx+1)]

# compute the sum of survival functions corresponding to achievements of disease remissions:
for (i in seq(2,2*r-1,2)){
  Si = summary(survfit(Surv(E[,i],R[,i])~1)) # a survival estimate corresponding to the achievement of the ith disease remission
  Si_days = stretch(list(x=Si$time,y=Si$surv),max(LastContact)) # assignment of survival estimates to each day of the follow-up
  SumS=SumS+Si_days$y[1:(maxx+1)]
}

# compute the sum of survival functions corresponding to losses of disease remissions:
for (i in seq(3,2*r,2)){
  Si_star = summary(survfit(Surv(E[,i],R[,i])~1)) # a survival estimate corresponding to the loss of the ith disease remission
  Si_star_days = stretch(list(x=Si_star$time,y=Si_star$surv),max(LastContact)) # assignment of survival estimates to each day of the follow-up
  SumS_star=SumS_star+Si_star_days$y[1:(maxx+1)]
}

clfs.pest <- list(x=0:maxx,y=S1_star+SumS_star-SumS) # estimation of the current leukaemia-free survival function as S1_star+SumS_star-SumS

}

