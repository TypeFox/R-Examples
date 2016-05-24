########################################################
##################### decodonise64 #####################
########################################################
#Copyright: PGE Toolbox (Cai et. al) (Matlab)
#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

decodonise64 <- function (x){

x  <- as.matrix(x)

n  <- dim(x)[1]
m  <- dim(x)[2]
y  <- 5*matrix(1,n,(m*3))
T  <- codontable()
TABLE <- T$Protein
CODON <- T$Triplets
icode <- 1
stops <- which(TABLE[icode,]=="*")
# TABLE wir hier garnicht genutzt !!!

for(i in 1:n){
 for(j in 1:m){
  c<- CODON[x[i,j],]
  #print(c)
  y[i,(j-1)*3+1]<-c[1]
  y[i,(j-1)*3+2]<-c[2]
  y[i,(j-1)*3+3]<-c[3]
 }
}
return(y) 

}
