statvc3levelt1 <-
function(y,cluster,m1,m2,n2,p,weight){

mat=matrix(0,p,p)
w=0

if(weight=="observation")
{
for(i in 1:m1)
	{
	if(m2[i]>1)
		{
		for(j in 1:(m2[i]-1))
			{
		for(k in (j+1):m2[i])	
			{
			yclj=y[cluster[[i]][[j]],,drop=F]
			yclk=y[cluster[[i]][[k]],,drop=F]
			wcur=(n2[[i]][j]+n2[[i]][k])/2	
			w=w+wcur*n2[[i]][j]*n2[[i]][k]
			sumjk=as.matrix(apply(yclj,2,sum)) %*% t(as.matrix(apply(yclk,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2
			}
			}
		}   # end if(m2[i]>1)
	}
}


else if(weight=="pair")
{
for(i in 1:m1)
	{
	if(m2[i]>1)
		{
		for(j in 1:(m2[i]-1))
			{
		for(k in (j+1):m2[i])	
			{
			yclj=y[cluster[[i]][[j]],,drop=F]
			yclk=y[cluster[[i]][[k]],,drop=F]			
			wcur=1	
			w=w+wcur*n2[[i]][j]*n2[[i]][k]
			sumjk=as.matrix(apply(yclj,2,sum)) %*% t(as.matrix(apply(yclk,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2
			}
			}
		}   # end if(m2[i]>1)
	}
}

else if(weight=="cluster")
{
for(i in 1:m1)
	{
	if(m2[i]>1)
		{
		for(j in 1:(m2[i]-1))
			{
		for(k in (j+1):m2[i])	
			{
			yclj=y[cluster[[i]][[j]],,drop=F]
			yclk=y[cluster[[i]][[k]],,drop=F]	
			wcur=(n2[[i]][j]*n2[[i]][k])		
			w=w+wcur*n2[[i]][j]*n2[[i]][k]
			sumjk=as.matrix(apply(yclj,2,sum)) %*% t(as.matrix(apply(yclk,2,sum)))
			mat=mat+wcur*(sumjk+t(sumjk))/2
			}
			}
		}   # end if(m2[i]>1)
	}
}



list(mat,w)

}
