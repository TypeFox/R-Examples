#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified
########################################################
######################## PsPn ##########################
########################################################

PsPn <- function(seq,subModel=TRUE){

icode <- 1
m <-dim(seq)[2]
n <-dim(seq)[1]

N <- getTable()

ns <- N$SynDif  # syndif
na <- N$AsynDif # asyndif

cseq <- codonise64(seq) #codiere sequenz in triplets
cn   <- dim(cseq)[1]
cm   <- dim(cseq)[2]

Ps <- 0
Pn <- 0


for(k in 1:cm){
   # Spaltenweise ueber die codonise64 Matrix
   csite  <- cseq[ ,k,drop=FALSE]
   # csite <- as.matrix(csite)
   H      <- counthaplotype(csite)
   numHap <- H$numHap
   sizHap <- H$sizHap
   seqHap <- H$seqHap
   
 if(numHap >1){ # wenn nicht alle gleich in der Spalte von codierter Triplet Matrix
   
#----------------------------------------------------
  if(subModel){  ##  Mit Substitutionsmodel
   
   C     <- i_codonsynnonsyn(seqHap,ns,na)
   cps   <- C$cps
   cpn   <- C$cpn
   sss   <- decodonise64(seqHap)
   m_num <- i_mutnum(sss)
     
    if((cps+cpn)==0){
      xps <- 0
      xpn <- 0
    }else{
      xps<-cps*m_num/(cps+cpn)
      xpn<-cpn*m_num/(cps+cpn)
    }  
   Ps<-Ps+xps
   Pn<-Pn+xpn  
  }# End SubModel  
#----------------------------------------------------

 #if(!subModel){ #  Ohne Substitutionsmodel
 #C  <- synornonsyn(seqHap)
 # if(C=="synonym")   {Ps <- Ps + 1}
 # if(C=="nonsynonym"){Pn <- Pn + 1}
 #} 
                           
 }# numHap >1
}# End for 
   
return(list(Ps=Ps,Pn=Pn))
}

# SubFunctions ###########
##########################
i_codonsynnonsyn <- function (csite,ns,na){

ns    <- as.matrix(ns)
na    <- as.matrix(na)
csite <- as.matrix(csite)
n     <- dim(csite)[1]
m     <- dim(csite)[2]
cps   <- 0
cpn   <- 0

#ns and na are a 64*64 Matrix

for (i in 1:(n-1)){
 for (j in (i+1):n){
    x <- csite[i]
    y <- csite[j]
    xx<-ns[x,y] #
    yy<-na[x,y]

     if(xx>0){
       cps <- cps + xx
     }
     if(yy>0){
       cpn <- cpn +yy
     }

 }
}
return(list(cps=cps,cpn=cpn))
}


i_mutnum <- function (seq){
#seq ist decodiert. gehe ueber alle drei Spalten

 m_num <- 0
for(k in 1:3){
 thissite <-seq[ ,k]
 x<-length(sort(unique(thissite)))-1
 m_num <- m_num +x
}

return(m_num)
}
