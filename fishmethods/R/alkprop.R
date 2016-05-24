alkprop<-function(x){
  if(is.null(x)) 
       stop ("No data")
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
	var1<-theta
	for(i in 1:length(var1[,1])){                                            
  	    for(j in 3:as.numeric(nages+2)){
            if(var1$sumA[i]>1) var1[i,j]<-(var1[i,j]*(1-var1[i,j])*var1$alpha[i]^2)/(var1$sumA[i]-1)
            if(var1$sumA[i]<=1) var1[i,j]<-0
         }
      }
	var2<-theta
	for(j in 3:as.numeric(nages+2)){                                        
   	    for(i in 1:length(var2[,1])){
            var2[i,j]<-(var2$alpha[i]*(var2[i,j]-est_theta_a[j-2])^2)/(sum(var2[,2]))

           }
      }
 	SE_theta_a<-sqrt(colSums(var1[,3:as.numeric(nages+2)])+ 
               colSums(var2[,3:as.numeric(nages+2)]))
	output<-as.data.frame(cbind(est_theta_a,SE_theta_a))
	output$CV<-ifelse(is.nan(SE_theta_a/est_theta_a),0,SE_theta_a/est_theta_a)
	names(output)<-c("prop","SE","CV")
	output<-list(output);names(output)<-c("results")
	return(output)
}