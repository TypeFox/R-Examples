##############################################################################
#
#     Sample Size From Age-Length Keys
#                           Fixed and Proportional Sampling
#                 Quinn and Deriso (1999: pp. 308-309)
##############################################################################

alkss<-function(x,lss=NULL,cv=NULL,allocate=1){
      if(is.null(x)) 
         stop ("No data")
      if(is.null(lss)) 
         stop ("length sample size is not specified") 
      if(is.null(cv)) 
         stop ("coefficient of variation is not specified") 
     
        nages<-ncol(x)-2
	alpha<-x[,2]/sum(x[,2])				
	x$sumA<-rowSums(x[,3:as.numeric(nages+2)])
	thetala<-x                                     
	for(i in 1:length(thetala[,1])){                                            
        for(j in 3:as.numeric(nages+2)){
           if(thetala$sumA[i]>0) thetala[i,j]<-thetala[i,j]/thetala$sumA[i]
           if(thetala$sumA[i]==0) thetala[i,j]<-0
        }
      }

	theta<-cbind(thetala,alpha)
	ta<-theta
	for(i in 1:length(ta[,1])){                                             
  	   for(j in 3:as.numeric(nages+2)){
            ta[i,j]<-ta[i,j]*ta$alpha[i]
         }
      }
	est_theta_a<-colSums(ta[,3:as.numeric(nages+2)])

	if(allocate==1){
 	  var1<-theta
        for(i in 1:length(var1[,1])){                                            
    		for(j in 3:as.numeric(nages+2)){
                var1[i,j]<-var1[i,j]*(1-var1[i,j])*var1$alpha[i]

            }
        }
       Sa<-colSums(var1[,3:as.numeric(nages+2)])
     }
	if(allocate==2){
 	   var1<-theta
  	   for(i in 1:length(var1[,1])){                                            
            for(j in 3:as.numeric(nages+2)){
                var1[i,j]<-var1[i,j]*(1-var1[i,j])*var1$alpha[i]^2
            }
         }
        Sa<-colSums(var1[,3:as.numeric(nages+2)])*length(var1[,2])
      }

	var2<-theta
	for(j in 3:as.numeric(nages+2)){                                        
         for(i in 1:length(var2[,1])){
            var2[i,j]<-var2$alpha[i]*(var2[i,j]-est_theta_a[j-2])^2
         }
         Ba<-colSums(var2[,3:as.numeric(nages+2)])

      }
	SS<-as.data.frame(round(Sa/(est_theta_a^2*cv^2-Ba/lss),0))
      names(SS)<-c(" ")
      dd<-c(paste("Age sample sizes for L=",lss,",cv=",cv,ifelse(allocate==1,
           ", Allocation: Proportional",", Allocation: Fixed")))
      SS<-list(dd,SS)
	names(SS)<-c("label","n")
      return(SS)
}