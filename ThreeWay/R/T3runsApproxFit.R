T3runsApproxFit <-
function(X,n,m,p,maxa,maxb,maxc){

X=as.matrix(X)
PCsup3=pcasup3(X,n,m,p)		
A=PCsup3$A
B=PCsup3$B
C=PCsup3$C
rr1=n
if (n>m*p){
	rr1=m*p
	A=A[,1:rr1]
}
rr2=m
if (m>n*p){
	rr2=n*p
	B=B[,1:rr2]
}
rr3=p
if (p>m*n){
	rr3=m*n
	C=C[,1:rr3]
}
Z=permnew(t(A)%*%X,rr1,m,p)
Z=permnew(t(B)%*%Z,rr2,p,rr1)
Z=permnew(t(C)%*%Z,rr3,rr1,rr2)
HH=Z^2/sum(X^2)*100       # HH contains squared core values, divided by SS(X), multiplied by 100
if (maxa>rr1){
	maxa=rr1
}
if (maxb>rr2){
	maxb=rr2
}
if (maxc>rr3){
	maxc=rr3
}
conv=1e-6
o=0
for (r1 in 1:maxa){
	for (r2 in 1:maxb){
		for (r3 in 1:maxc){
			if (((r1*r2)>=r3) & ((r1*r3)>=r2) & ((r2*r3)>=r1)){
				o=o+1
				H=permnew(HH[1:r1,],r1,rr2,rr3)
				H=permnew(H[1:r2,],r2,rr3,r1)
				H=H[1:r3,]
				fp=sum(H)
				if (o==1){
					out=matrix(c(r1,r2,r3,fp),1,4)
				} else{
					out=rbind(out,c(r1,r2,r3,fp))
				}
			}
		}
	}
}
if (nrow(out)>1){
	s=colSums(t(out[,1:3]))
} else{
	s=sum(out[,1:3])
}
out=cbind(out,s)
out=out[ord(s)$a,]
if (is.matrix(out)==FALSE){
	out=as.matrix(t(out))
}
rownames(out)=rep("T3",length=dim(out)[1])
colnames(out)=c("P","Q","R","Fit (%) ","P+Q+R")
return(out)
}
