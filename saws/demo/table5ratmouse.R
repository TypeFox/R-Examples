cat("This program recreates the Rat and Mouse combined analysis from Table 5
of Fay, Freedman, Clifford, and Midthune (1997, Cancer Research 57: 39379-3988
but using the slightly different methods of Fay and Graubard (2001, Biometrics 1198-1206")


library(saws)
data(dietfat)
set<-as.factor(dietfat$SET)
n<-dietfat$N
length(n)
m<-dietfat$NTUM
n3<-dietfat$PN3
n6<-dietfat$PN6
sat<-dietfat$PZERO
mono<-dietfat$PMONO
restrict<-dietfat$RESTRICT
ln6<-n6
ln6[ln6>4]<-4
un6<- n6-4
un6[un6<0]<-0
## create matrix for clogist analysis
x<-matrix(c(restrict,n3,sat,mono,ln6,un6),ncol=6,dimnames=list(NULL,c("restrict","n3","sat","mono","ln6","un6")))
x[1,]

cout<-clogistCalc(n,m,x,set)

sm<-saws(cout,method="dm")
s1<-saws(cout,method="d1")
s2<-saws(cout,method="d2")
s3<-saws(cout,method="d3")
s4<-saws(cout,method="d4")
s5<-saws(cout,method="d5")

str(sm)
s1
s2
s3
s4
s5

n3PUFA.pvalues<-c(dm=sm$p.value[2],d1=s1$p.value[2],
    d2=s2$p.value[2],d3=s3$p.value[2],d4=s4$p.value[2],d5=s5$p.value[2])

cat("pvalues for n-3 PUFA, see Biometrics (2001, p. 1202)")
round(n3PUFA.pvalues,4)





