unor <-
function(a,binsize,out=c("symb","num"))
{#discretizacion 1R. It has only one restriction that does not
#allow a discretization in more than 20 intervals.
#This can be fixed by considering the maximum number 
# of intervals as a parameterof the function. 
    b=a[order(a[,1]),]
    n=dim(a)[1]
    cumcut=c(1,binsize)
    mayor=moda(b[cumcut[1]:cumcut[2],2])
    j=binsize
    while(j<n){j=j+1;if(b[j,2]!=mayor)break}
    j1=j-1
#    print(j1)
    while(j1<n){j1=j1+1;if(b[j1,1]!=b[j1-1,1])break}
    cumcut[2]=j1-1
#   print(cumcut)
    for(jj in 2:20)
    {maxclass= moda(b[(cumcut[jj]+1):(cumcut[jj]+binsize),2])[1]
     mayor=c(mayor,maxclass)
     ind=cumcut[jj]+binsize
     while(ind<n){ind=ind+1;if(b[ind,2]!=maxclass)break}
     tempo=ind-1
     while(tempo<n){tempo=tempo+1;if(b[tempo,1]!=b[tempo-1,1])break}
     cumcut=c(cumcut,tempo-1)
     if(tempo>n)break
    }
 #   print(cumcut)
    nm=length(mayor)
#merging
    finalcut=cumcut[((1:(nm-1))[abs(diff(mayor))>=1])+1] 
a1=as.vector(a[,1])
cutpoints=b[finalcut,1]
cutpoints=c(-Inf,cutpoints,Inf)
if(out=="num")
a1=cut(a1,cutpoints,labels=F)
else{a1=cut(a1,cutpoints)}
a1
}
