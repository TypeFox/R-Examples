normvari <-
function(A){		

r=ncol(A)
m=nrow(A)
d=SUM(A)$row  		# check on near zero rows, which will not be weighted
for (i in 1:m){ 
	if(d[i]<1e-16){		
		d[i]=1
	}
}
A=A/((d^.5)%*%matrix(1,1,r))
VARIM=varim(A)
B=((d^.5)%*%matrix(1,1,r))*VARIM$B
out=list()
out$B=B
out$T=VARIM$T
out$f=VARIM$f
return(out)
}
