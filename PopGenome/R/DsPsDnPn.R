########################################################
###################### DsPsDnPn ########################
########################################################
#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

DsPsDnPn <- function (aln,vek1,vek2){

Ds<-0
Ps<-0
Dn<-0
Pn<-0
Ls<-0
Ln<-0

#SP <- splitten(aln,vek1,vek2)
#seq1 <- openfile2(SP$aln1)
#seq2 <- openfile2(SP$aln2)

aln1 <- aln[vek1,,drop=FALSE]
aln2 <- aln[vek2,,drop=FALSE]

seq1 <- aln1 #code(aln1) codiert wird in popgen
seq2 <- aln2 #code(aln2)

m<-dim(aln)[2]
n<-dim(aln)[1]

P1 <- PsPn(seq1) 

P2 <- PsPn(seq2)


Ps1 <- P1$Ps
Pn1 <- P1$Pn
Ps2 <- P2$Ps
Pn2 <- P2$Pn

Ps <- Ps1 + Ps2
Pn <- Pn1 + Pn2

D <- DsDn(seq1,seq2)
Ds <-D$Ds
Dn <-D$Dn 


return(list(Ds=Ds,Ps=Ps,Dn=Dn,Pn=Pn))
}
