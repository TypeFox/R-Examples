#Function: ghap.ancestral
#License: GPLv3 or later
#Modification date: 19 Jen 2015
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com
#Description: Assign the aplotype to one of the parental population 

ghap.ancestral<-function(
  hapstats.test,
  hapstats.parent1,
  hapstats.parent2,
  freq=0.05,
  prob.assign=0.60
){
  

  # Prepare test data
  hapstats.test$HAPLOtmp<-paste(hapstats.test$BLOCK,hapstats.test$ALLELE,sep="_")
  hapstats.test<-hapstats.test[,c("BLOCK","CHR","BP1","BP2","ALLELE","FREQ","HAPLOtmp")]
  colnames(hapstats.test)[6]<-"FREQ.TEST"
  
  # Prepare parent1 data
  hapstats.parent1$HAPLOtmp<-paste(hapstats.parent1$BLOCK,hapstats.parent1$ALLELE,sep="_")
  hapstats.parent1<-hapstats.parent1[,c("HAPLOtmp","FREQ")]
  colnames(hapstats.parent1)[2]<-"FREQ.PARENT1"
  
  # Prepare parent1 data
  hapstats.parent2$HAPLOtmp<-paste(hapstats.parent2$BLOCK,hapstats.parent2$ALLELE,sep="_")
  hapstats.parent2<-hapstats.parent2[,c("HAPLOtmp","FREQ")]
  colnames(hapstats.parent2)[2]<-"FREQ.PARENT2"
  
  # Merge data frames
  out<-merge(x = hapstats.test,y=hapstats.parent1,by.x = "HAPLOtmp",by.y="HAPLOtmp",all.x=TRUE)
  out<-merge(x = out,y=hapstats.parent2,by.x = "HAPLOtmp",by.y="HAPLOtmp",all.x=TRUE)
  out[is.na(out)==TRUE]<-0
  rm(hapstats.test,hapstats.parent1,hapstats.parent2)
  
  # Calculate probability
  out$PROB.PARENT1 <- out$FREQ.PARENT1/(out$FREQ.PARENT1+out$FREQ.PARENT2)
  out$PROB.PARENT2 <- out$FREQ.PARENT2/(out$FREQ.PARENT1+out$FREQ.PARENT2)
  out$PROB.PARENT1[is.na(out$PROB.PARENT1)==TRUE]<-0
  out$PROB.PARENT2[is.na(out$PROB.PARENT2)==TRUE]<-0
  out$PROB.PARENT1[out$FREQ.PARENT1 <= freq]<-0
  out$PROB.PARENT2[out$FREQ.PARENT2 <= freq]<-0
  
  # Assign the haplotype to one of the parental populations or set as unknow 
  out$ORIGIN <- "UNK"
  out$ORIGIN[out$PROB.PARENT1 > prob.assign] <- "PARENT1"
  out$ORIGIN[out$PROB.PARENT2 > prob.assign] <- "PARENT2"
  out <- out[,colnames(out) != "HAPLOtmp"]
  
  # Output
  return(out[order(out$BP1,out$BP2,out$ALLELE),])
  
}
