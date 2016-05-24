"srswor1" <-
function(n,N)
{j=0
 s=numeric(N)
 for(k in 1:N) if(runif(1)<(n-j)/(N-k+1)) {j=j+1;s[k]=1;}
 s
}

