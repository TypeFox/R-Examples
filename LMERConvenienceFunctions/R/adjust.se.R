adjust.se<-function(df,ncomp,alpha=0.05,two.tailed=TRUE){
        if(two.tailed){
             return(abs(qt(alpha/2/ncomp,df)))
        }else{
             return(abs(qt(alpha/ncomp,df)))
        }
}
