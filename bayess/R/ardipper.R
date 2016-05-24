ardipper <- function(nsimu,n1,c2,c3,r2,q1){
# ACCEPT-REJECT ALGORITHM FOR THE OPEN POPULATION MODEL

echan=1:nsimu
for (i in 1:nsimu)
{
test <- T
while (test==T) {
y=rbinom(1,n1,q1)
if (runif(1) <= thresh(y,n1,c2,c3,r2,q1))
{
test <- F
echan[i] <- y
}
}
}
echan
}
