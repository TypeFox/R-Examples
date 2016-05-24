###################################################################
#---------------------------------------------------------------#
# empirical Response Function
#---------------------------------------------------------------#
erf<- function(theta,  data, whichItem, min=-3, max=3, Ncuts=12){
  
  breakPoints<-seq(min,max,by=(max-min)/Ncuts)
  centers<-breakPoints+(max-min)/(2*Ncuts)
  # remove last center
  centers<-centers[-length(centers)]
  
  probs<-unlist(lapply(split(data[,whichItem],cut(theta,breaks=breakPoints)),mean))
  group.N<-unlist(lapply(split(data[,whichItem],cut(theta,breaks=breakPoints)),length))
  se.p<- sqrt(( probs * (1-probs))/group.N) 
  #     for (i in 1:Ncuts) { 
  #         lines( c((1:Ncuts)[i],(1:Ncuts)[i]),c(probs[i],(probs[i]+1.96*se.p[i]) ),lwd=line.width) 
  #         }
  #     for (i in 1:Ncuts) { 
  #         lines( c((1:Ncuts)[i],(1:Ncuts)[i]),c(probs[i],(probs[i]-1.96*se.p[i]) ),lwd=line.width) 
  #         }         
  list(probs=probs, centers = centers, Ni=group.N, se.p = se.p)
}    


