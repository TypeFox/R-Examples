mgpd_data <-
function(xdat,thr=rep(0,5)){
n               = dim(xdat)[1]
d               = dim(xdat)[2]
if(!(d %in% c(2:5))) print("Error message: invalid dimension.") else{
potdat          = NULL
id              = rep(0,n)
for(i in 1:d) potdat         = cbind(potdat,xdat[,i]-thr[i])
for(i in 1:d) id             = id+(potdat[,i]>0)
potdat          = potdat[id>0,]
}
potdat          = as.data.frame(potdat)
potdat
}
