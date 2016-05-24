varim <-
function(A){

conv=.000001
T=matrix(0,ncol(A),ncol(A))
r=ncol(A)
m=nrow(A)
for (i in 1:r){
	T[i,i]=1
}
B=A
f=sum(((A*A)-matrix(1,m,1)%*%colSums(A*A)/m)^2)
fold=f-2*conv*f
if (f==0){
	fold=-conv
}
iter=0

while ((f-fold)>(f*conv)){
	fold=f
	iter=iter+1
	for (i in 1:(r-1)){
		for (j in (i+1):r){
			x=B[,i]
			y=B[,j]
			xx=T[,i]
			yy=T[,j]
			u=x^2-y^2
			v=2*x*y
			u=u-matrix(1,m,1)*sum(u)/m
			v=v-matrix(1,m,1)*sum(v)/m  
			a=2*sum(u*v)
			b=sum(u^2)-sum(v^2)
			c=(a^2+b^2)^.5
			if (a>=0){
				sign=1
			}
			if (a<0){
				sign=-1
			}
			o=1
			if (c<.00000000001){
			    cos=1
				sin=0
				o=0
			}
			if (c>=.00000000001){
				vvv=-sign*((b+c)/(2*c))^.5
				sin=(.5-.5*vvv)^.5
				cos=(.5+.5*vvv)^.5
			}
			v=cos*x-sin*y
			w=cos*y+sin*x
			vv=cos*xx-sin*yy
			ww=cos*yy+sin*xx
			if (o==1){
			if (vvv>=0){	     # prevent permutation of columns
				B[,i]=v
				B[,j]=w
				T[,i]=vv
				T[,j]=ww
			}
			if (vvv<0){
				B[,i]=w
				B[,j]=v
				T[,i]=ww
				T[,j]=vv
			}
			}
		}
	}
	f=sum(((B*B)-matrix(1,nrow(B),1)%*%colSums(B*B)/nrow(B))^2)
}
out=list()
out$f=f
out$B=B
out$T=T
return(out)
}
