stab.fw <-
function(y,Gen,Env,times,Rep,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(X))X=NULL
   if(Rep==TRUE){
       dat=GetGEMean(y,Gen,Env,X)
       y1=dat$y
       if(is.null(X))X1=NULL
       else X1=dat[,-c(1:3)]
       Gen1=dat$Gen
       Env1=dat$Env
       result=fw_reg(y1,Gen1,Env1,times,X1,alpha)
   }
   else result=fw_reg(y,Gen,Env,times,X,alpha)


   return(result)
}
