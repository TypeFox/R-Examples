sdd <-
function(x){
 #COMPUTES SD OF A dichotomous variable
 p<-mean(x)
 sdd<-sqrt((p*(1-p)))   
    return(sdd)}
