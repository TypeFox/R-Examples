absspp <- function (dat, crit=2, disp=TRUE) { 
# ABSSPP Deletes rare variables/replicates
# Deletes variables with no variance or less than CRIT
# observations. If DISP is FALSE, messages are not displayed.
# SPPABS returns the columns of DAT that have been excluded.
count 	<- matrix(1,1,NROW(dat)) %*% as.matrix(dat>0)  # summation over the rows 

p	 	<-  (count<crit | diag(var(dat))==0)
sppabs 	<-  (1:NCOL(dat))[p]
if (any(p)){
   if (disp){
	  cat("\nvariables ignored:\n")
	  cat(sppabs)
	  cat("\n\n")
   }
    dat <- dat[,!p]

}
return(list(dat=dat, sppabs=sppabs))

}
