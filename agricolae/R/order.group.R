`order.group` <-
function(trt,means,N,MSerror,Tprob,std.err,parameter=1, snk=0, DFerror=NULL,alpha=NULL,sdtdif=NULL,vartau=NULL,console ) {
replications <- N
#N<-rep(1/mean(1/N),length(means))
n<-length(means)
z<-data.frame(trt,means,N=N,std.err,replications)
letras<-c(letters[1:26],LETTERS[1:26],1:9)
# ordena tratamientos
w<-z[order(z[,2],decreasing = TRUE), ]
M<-rep("",n)
k<-1
j<-1
k<-1
cambio<-n
cambio1<-0
chequeo=0
M[1]<-letras[k]
N<-w$N
while(j<n) {
chequeo<-chequeo+1
if (chequeo > n) break
for(i in j:n) {
nx<-abs(i-j)+1
if (nx==1) nx=2
if(snk ==1 ) Tprob <- qtukey(p=1-alpha,nmeans=nx, df=DFerror)
if(snk ==2 ) Tprob <- qtukey(p=(1-alpha)^(nx-1),nmeans=nx, df=DFerror)
if(snk ==7 ) {
  if(nx <= n-2) Tprob<- qtukey(p=(1-alpha)^(nx/n),nx,df=DFerror)
  if(nx > n-2) Tprob<- qtukey(p=(1-alpha),nx,df=DFerror)
}
if (!is.null(vartau))  {
ii0<-w[i, 1]
jj0<-w[j, 1]
if (snk == 3 | snk==4 ) sdtdif <- sqrt(vartau[ii0, ii0] + vartau[jj0, jj0] - 2 * vartau[ii0,jj0])
# DAU
if (snk == 5 | snk==6 ) sdtdif <- sqrt(vartau[ii0, jj0])
# PBIB and DAU  
if(snk == 3 | snk == 5 ) Tprob <- qt(p = (1 - alpha/2), df = DFerror)
if(snk == 4 | snk == 6 ) Tprob <- qtukey(p = (1 - alpha), nmeans=n,df = DFerror)
}
if(is.null(sdtdif))  minimo<-Tprob*sqrt(parameter*MSerror*(1/N[i]+1/N[j]))
if(!is.null(sdtdif)) minimo<-Tprob*sdtdif
s<-abs(w[i,2]-w[j,2])<=minimo
if(s) {
if(lastC(M[i]) != letras[k])M[i]<-paste(M[i],letras[k],sep="")
}
else {
k<-k+1
cambio<-i
cambio1<-0
ja<-j
for(jj in cambio:n) M[jj]<-paste(M[jj],sep="")
M[cambio]<-paste(M[cambio],letras[k],sep="")
for( v in ja:cambio) {
nx<-abs(v-cambio)+1
if(nx == 1)  nx=2
if(snk ==1 ) Tprob <- qtukey(p=1-alpha,nmeans=nx, df=DFerror)
if(snk ==2 ) Tprob <- qtukey(p=(1-alpha)^(nx-1),nmeans=nx, df=DFerror)
if(snk ==7 ) {
  if(nx <= n-2) Tprob<- qtukey(p=(1-alpha)^(nx/n),nx,df=DFerror)
  if(nx > n-2) Tprob<- qtukey(p=(1-alpha),nx,df=DFerror)
}

if (!is.null(vartau))  {
ii0<-w[i, 1]
jj0<-w[j, 1]
# PBIB
if (snk == 3 | snk==4 ) sdtdif <- sqrt(vartau[ii0, ii0] + vartau[jj0, jj0] - 2 * vartau[ii0,jj0])
# DAU
if (snk == 5 | snk==6 ) sdtdif <- sqrt(vartau[ii0, jj0])
# PBIB and DAU  
if(snk == 3 | snk == 5 ) Tprob <- qt(p = (1 - alpha/2), df = DFerror)
if(snk == 4 | snk == 6 ) Tprob <- qtukey(p = (1 - alpha), nmeans=n,df = DFerror)
}
if(is.null(sdtdif))  minimo<-Tprob*sqrt(parameter*MSerror*(1/N[i]+1/N[j]))
if(!is.null(sdtdif)) minimo<-Tprob*sdtdif
if(abs(w[v,2]-w[cambio,2])>minimo) {j<-j+1
cambio1<-1 
}
else break
}
break
}
}
if (cambio1 ==0 )j<-j+1
}
#-----
w <- data.frame(w, stat = M)
trt <- as.character(w$trt)
means <- as.numeric(w$means)
means1<-signif(means,4)
N <- as.numeric(w$N)
std.err <- as.numeric(w$std.err)
replications<-w$replications
cmax <- max(nchar(trt))
trt <- paste(trt, "                                      ")
trt <- substr(trt, 1, cmax)
output1<-cbind(trt,means)
rownames(output1)<-M
for (i in 1:n) {
	if(console)cat(rownames(output1)[i], "\t", output1[i,1], "\t", means1[i], "\n")
}
output <- data.frame(trt, means, M, N=replications, std.err)
invisible(output)
}
