`kendall` <-
function(data1,data2) {
n<-length(data1)
n2<-0;n1<-0;is<-0
n0<-n-1
for (j in 1:n0) { 
jj<-j+1
for (k in jj:n) {
a1<-data1[j]-data1[k]
a2<-data2[j]-data2[k]
aa<-a1*a2
if (! is.na(aa)) {
if (aa) {
n1<-n1+1
n2<-n2+1
if(aa > 0.0) is<-is+1
else is=is-1
} else { 
if (a1) n1<-n1+1
if (a2) n2<-n2+1
}
}
}
}
tau<-is/(sqrt(n1)*sqrt(n2))
#
#svar<-(4.0*n+10.0)/(9.0*n*(n-1.0))
#z<-tau/sqrt(svar)
#prob<-erfcc(abs(z)/1.4142136)
z<-is/sqrt(vark(data1,data2))
prob<- 2*pnorm(-abs(z))
return(list(stat=z,tau=tau,pvalue=prob))
}

