makeamat <-
function(xmat,m){
	n=length(xmat)/m
#sort data -- ascending
	sm=1e-10
	amat1=matrix(0,nrow=n*(n-1)/2,ncol=n)
	comp=matrix(0,nrow=n,ncol=n)
	nr=0
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			bigger=1;l=0
			while(l<m){
				l=l+1
				if(xmat[i,l]>xmat[j,l]+sm){bigger=0}
			}
			if(bigger==1){
				nr=nr+1
				amat1[nr,j]=1;amat1[nr,i]=-1
				comp[i,j]=1
			}
		}
	}
	amat2=amat1[1:nr,]
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			if(comp[i,j]==1){
				if(j<n){
				for(k in (j+1):n){
					if(comp[j,k]==1&comp[i,k]==1){
						comp[i,k]=2
					}
				}}
			}
		}
	}
	dump=1:nr<0
	id=0
	for(i in 1:(n-1)){
		for(j in (i+1):n){
			if(comp[i,j]>sm){
				id=id+1
				if(comp[i,j]==2){dump[id]=TRUE}
			}
		}
	}
	amat3=amat2[!dump,]
	ans=new.env()
	ans$amat=amat3
	ans
}
