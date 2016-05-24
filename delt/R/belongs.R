belongs<-function(obs,rec,coordi){
#Finds whether observation belongs to the rectangle
#
#obs is d-vector
#rec is 2*d-vector
#coordi is in 1:d
#
#Returns TRUE is obs is in rec, otherwise FALSE
#
#We need to check only whether coordi:s coordinate is in the
#interval
#
ans<-TRUE
if ((obs[coordi]<rec[2*coordi-1]) || (obs[coordi]>rec[2*coordi])) ans<-FALSE
return(ans)
}
