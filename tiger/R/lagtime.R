`lagtime` <-
function(x,y){
 if(var(x)==0 | var(y)==0){
    if(var(x)==0 & var(y)==0){
       return(0)
    } else {
       warning("unable to determine lagtime between constant and variable series")
       return(NA)
    }
 } 

 theRes<-ccf(x,y, plot=FALSE)
 return(theRes$lag[ min(which(theRes$acf==max(theRes$acf)))])
}

