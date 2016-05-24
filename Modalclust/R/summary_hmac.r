summary.hmac<-function(object,...)
{
	z<-object
	n<-z$n.cluster
	levels<-z$level
	modes<-z$mode

	ans<-z["terms"]
	class(ans)<-"summary.hmac"


        ans$level<- unique(levels)
        ans$cluster<- unique(n)	
        ans$modes<-modes



        
 cat("At levels", ans$level, " ") 
 cat("there are ", ans$cluster, " clusters respectively \n","\n")
 cat("The modes of each level of clusters are", "\n","\n")
	for(i in 1:max(ans$level)){
      cat("Modes at level ", i,"\n")
   if(ncol(modes[[1]])>1)  print(data.frame(ans$modes[[i]]))
    else    cat(ans$modes[[i]],"\n")
	cat("\n")
    }
  	
}
