########################################################
################## Counts haplotypes ###################
########################################################
# Copyright
# slightly modified from the PGE TOOLBOX.
#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified (unused)

counthaplotype <- function(seq){

 # die Funktion zaehlt die haploiden Datensaetze
 # zwei oder mehr gleiche Sequenzen untereinander --> haploid
 # uebernehme in seqHap nur die unterschiedlichen sequenzen

#  seqHap

n <- dim(seq)[1]
m <- dim(seq)[2]

if(n==1){

 numHap<-1
 sizHap<-1
 seqHap <-seq

 return(list(numHap=numHap,sizHap=sizHap,seqHap=seqHap))
}
  
      seq <- sortmatrix(seq) 
      seq <- as.matrix(seq)
      seqHap<-seq[1, ]
      curseq<-seq[1, ]
	    sizHap<-matrix(0,1,n)
	    numHap <-1
  
      for(i in 1:n){
	     if(sum(seq[i,]==curseq)!=m){
	      seqHap <-rbind(seqHap,seq[i,])
	      numHap <- numHap+1
	      curseq <- seq[i,]
	      }
       sizHap[numHap]<-sizHap[numHap]+1	
      }     
       sizHap<-sizHap[1:numHap]
       v <- sort(-sizHap,index.return=TRUE)
       sizHap <-  v$x
       sizHap <- -sizHap
       idx    <- v$ix
       seqHap <- as.matrix(seqHap)
       seqHap <- seqHap[idx, ]      
       seqHap <- as.matrix(seqHap)

return(list(numHap=numHap,sizHap=sizHap,seqHap=seqHap))

}
