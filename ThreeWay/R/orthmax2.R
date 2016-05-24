orthmax2 <-
function(A1,A2,gam1,gam2,conv){				

narg=nargs()
if (narg<5){
	conv=1e-6
}
m1=nrow(A1)
r=ncol(A1)
m2=nrow(A2)
r2=ncol(A2)
if (r2!=r){
	stop("Error! Column orders of A and B unequal!")
}

T=diag(r)
B1=A1
B2=A2
f=sum(B1^4)-gam1/m1*sum((colSums(B1^2))^2)+sum(B2^4)-gam2/m2*sum((colSums(B2^2))^2)

if (r>1){
	fold=f-2*conv*abs(f)
	if (f==0){
		fold=-conv
	}
	iter=0
	while (f-fold>abs(f)*conv){
		fold=f
		iter=iter+1
		for (i in 1:(r-1)){
			for (j in (i+1):r){
		
				# Jennrich & Clarkson
				xx=T[,i]
				yy=T[,j]
				a=0
				b=0
			
				# for A1
				x=B1[,i]
				y=B1[,j]
				x2=x^2
				y2=y^2
				a=a+(-gam1/m1)*(.25*(sum(x2-y2))^2 - (sum(x*y))^2) + .25*sum(x2^2 + y2^2 - 6*x2*y2)
				b=b+(-gam1/m1)*sum(x*y)*sum(x2-y2) + sum((x^3)*y - x*(y^3))

				# for A2
				x=B2[,i]
				y=B2[,j]
				x2=x^2
				y2=y^2
				a=a+(-gam2/m2)*(.25*(sum(x2-y2))^2 - (sum(x*y))^2) + .25*sum(x2^2 + y2^2 - 6*x2*y2)
				b=b+(-gam2/m2)*sum(x*y)*sum(x2-y2) + sum((x^3)*y - x*(y^3))      
				theta=0.25*atan2(b,a)
				cs=cos(theta)
				sn=sin(theta)
				x=B1[,i]
				y=B1[,j]
				B1[,i]=cs*x+sn*y
				B1[,j]=cs*y-sn*x
				x=B2[,i]
				y=B2[,j]
				B2[,i]=cs*x+sn*y      
				B2[,j]=cs*y-sn*x
				T[,i]=cs*xx+sn*yy
				T[,j]=cs*yy-sn*xx
			}
		}	
		f=sum(B1^4)-gam1/m1*sum((colSums(B1^2))^2)+sum(B2^4)-gam2/m2*sum((colSums(B2^2))^2)
	}
}
out=list()
out$B1=B1
out$B2=B2
out$T=T
out$f=f
return(out)
}
