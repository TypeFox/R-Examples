########################################################
######################### DsDn #########################
########################################################
#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified


DsDn <- function(seq1,seq2,subModel=TRUE){
icode <-1

n1 <-dim(seq1)[1]
m1 <-dim(seq1)[2]
n2 <-dim(seq2)[1]
m2 <-dim(seq2)[2]

if(m1!=m2){
stop("m1!=m2")
}

cseq1 <- codonise64(seq1)
cseq2 <- codonise64(seq2)

cn1 <- dim(cseq1)[1]
cm1 <- dim(cseq1)[2]
cn2 <- dim(cseq2)[1]
cm2 <- dim(cseq2)[2]

N  <- getTable()
ns <- N$SynDif
na <- N$AsynDif
Ds <-0
Dn <-0

for(k in 1:cm1){
 csite1<-cseq1[,k]
 csite2<-cseq2[,k]
 # intersect(a,b) Indizes der Elemente in a die in b vorkommen

 if(length(sort(intersect(csite1,csite2)))==0){
   csite1 <- as.matrix(csite1)
   C <- counthaplotype(csite1)
   numHap1 <- C$numHap
   sizHap1 <- C$sizHap     
   seqHap1 <- C$seqHap
   csite2 <-  as.matrix(csite2)
   C<- counthaplotype(csite2)
   numHap2 <- C$numHap
   sizHap2 <- C$sizHap     
   seqHap2 <- C$seqHap

  if(subModel){ # Mit Substitutionsmodell
   D   <- i_codonsynnonsyn2(seqHap1,seqHap2,ns,na)
   cds <- D$cps
   cdn <- D$cpn

   sss1  <- decodonise64(seqHap1)
   sss2  <- decodonise64(seqHap2)
   m_num <- i_mutnum2(sss1,sss2)
    
   if((cds+cdn)==0){
     xds <-0
     xdn <-0 
   }else{
     xds<-(cds*m_num)/(cds+cdn)
     xdn<-(cdn*m_num)/(cds+cdn) 
   }  
  Ds<-Ds+xds
  Dn<-Dn+xdn
 }# Ende SubModel
 
 #if(!subModel){ # Ohne Substitutionsmodell
 #  seqHapall  <- rbind(seqHap1,seqHap2)
 #  D          <- synornonsyn(seqHapall)
 #  Ds         <- Ds + 1
 #  Dn         <- Dn + 1
 #}

 } 
}


return(list(Ds=Ds,Dn=Dn))

}

i_codonsynnonsyn2 <- function(csite1,csite2,ns,na){

   ns <- as.matrix(ns)
   na <- as.matrix(na)
   csite1 <- as.matrix(csite1)
   csite2 <- as.matrix(csite2)
   
   n1 <- dim(csite1)[1]
   m1 <- dim(csite1)[2]
   n2 <- dim(csite2)[1]
   m2 <- dim(csite2)[2]
cps <-0
cpn <-0
 
 for(i in 1:n1){
   for(j in 1:n2){
    x <- csite1[i]
    y <- csite2[j]
   xx <- ns[x,y]
   yy <- na[x,y]
   if(xx>0){cps<-cps+xx}
   if(yy>0){cpn<-cpn+yy}

   }
 }  

return(list(cps=cps,cpn=cpn))
}


i_mutnum2 <- function (seq1,seq2){

m_num <- 0
 for(k in 1:3){
  thissite1<-seq1[,k]
  thissite2<-seq2[,k]
   if(length(sort(intersect(thissite1,thissite2)))==0){
     m_num <- m_num +1
   }
 }
return(m_num) 
}
