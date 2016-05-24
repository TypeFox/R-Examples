# binseqtest package depends on clinfun 
# so it is included automatically 
#require(clinfun)
#library(clinfun)
# to test the ph2simon function in the clinfun R package, we 
# reproduce Tables 1 and 2 of Simon (1989)

getResults<-function (x, ...) 
{
    # copied from print.ph2simon
    xout <- x$out
    nmax <- x$nmax
    n <- nrow(xout)
    nopt <- ((1:n)[xout[, 5] == min(xout[, 5])])[1]
    xopt <- xout[c(nopt, 1), ]
    dimnames(xopt)[[1]] <- c("Optimal", "Minimax")
    xopt
}


P0<-rep(c(.05,.1,.2,.3,.4,.5,.6,.7),each=3)
P1<-P0+.2
ALPHA<-rep(c(.1,.05,.05),8)
BETA<-rep(c(.1,.2,.1),8)

TABLE1<-matrix(NA,length(P0),16,dimnames=list(NULL,c("p0","p1","alpha","beta",
    paste("Opt",c("r1","n1","r","n","EN","PET")), 
    paste("mm",c("r1","n1","r","n","EN","PET"))) ) )
TABLE1[,1]<-P0
TABLE1[,2]<-P1
TABLE1[,3]<-ALPHA
TABLE1[,4]<-BETA
for (i in 1:length(P0)){
    out<-ph2simon(P0[i], P1[i], ALPHA[i], BETA[i])
    x<-getResults(out)
    TABLE1[i,5:10]<-x[1,]
    TABLE1[i,11:16]<-x[2,]
}
TABLE1

## Table 2 is same except P1=P0+.15
P1<-P0+.15

TABLE2<-matrix(NA,length(P0),16,dimnames=list(NULL,c("p0","p1","alpha","beta",
    paste("Opt",c("r1","n1","r","n","EN","PET")), 
    paste("mm",c("r1","n1","r","n","EN","PET"))) ) )
TABLE2[,1]<-P0
TABLE2[,2]<-P1
TABLE2[,3]<-ALPHA
TABLE2[,4]<-BETA
for (i in 1:length(P0)){
    out<-ph2simon(P0[i], P1[i], ALPHA[i], BETA[i])
    x<-getResults(out)
    TABLE2[i,5:10]<-x[1,]
    TABLE2[i,11:16]<-x[2,]
}
TABLE2


