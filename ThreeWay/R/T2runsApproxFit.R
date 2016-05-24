T2runsApproxFit <-
function(X,n,m,p,maxa,maxb,maxc,model){

X=as.matrix(X)
PCsup2=pcasup2(X,n,m,p,model)		
A=PCsup2$A
B=PCsup2$B
C=PCsup2$C
rr1=n
rr2=m
rr3=p
if (model==1){
	Y=permnew(X,n,m,p)
	Y=permnew(Y,m,p,n)
	C=eigen(Y%*%t(Y))$vectors
}
if (model==2){
	Y=permnew(X,n,m,p)
	B=eigen(Y%*%t(Y))$vectors
}
if (model==3){
	A=eigen(X%*%t(X))$vectors
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
if (model==1){
	for (r1 in 1:maxa){
		for (r2 in 1:maxb){
			o=o+1
			H=permnew(HH[1:r1,],r1,rr2,rr3)
			H=permnew(H[1:r2,],r2,rr3,r1)
			H=H[1:rr3,]
			fp=sum(H)
			if (o==1){
				out=matrix(c(r1,r2,rr3,fp),1,4)
			} else{
				out=rbind(out,c(r1,r2,rr3,fp))
			}
		}
	}
}
if (model==2){
	for (r1 in 1:maxa){
		for (r3 in 1:maxc){
			o=o+1
				H=permnew(HH[1:r1,],r1,rr2,rr3)
				H=permnew(H[1:rr2,],rr2,rr3,r1)
				H=H[1:r3,]
			fp=sum(H)
			if (o==1){
				out=matrix(c(r1,rr2,r3,fp),1,4)
			} else{
				out=rbind(out,c(r1,rr2,r3,fp))
			}
		}
	}
}
if (model==3){
	for (r2 in 1:maxb){
		for (r3 in 1:maxc){
			o=o+1
			H=permnew(HH[1:rr1,],rr1,rr2,rr3)
			H=permnew(H[1:r2,],r2,rr3,rr1)
			H=H[1:r3,]
			fp=sum(H)
			if (o==1){
				out=matrix(c(rr1,r2,r3,fp),1,4)
			} else{
				out=rbind(out,c(rr1,r2,r3,fp))
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
