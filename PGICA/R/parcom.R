parcom <-
function(subj,w,files,m,n,fegfe,ind,X)
{
    #	added following line to avoid package notes
#	X=get("X");
    w1=t(w[[subj]])
    winv=solve(w1)
    load(files[subj])
    x1=t(X)
    
    gradn=matrix(0,m*m,n)
    for(k in 1:m){
      for(j in 1:m){
        gradn[((k-1)*m+1)+j-1,]=x1[,j]*fegfe[,k]+t(winv)[j,k]
      }
    }
    grad=apply(gradn,1,sum)
    mind=floor(n/ind)
    hes=0
    for(i in 1:mind){
    hes=hes+gradn[,((i-1)*ind+1):(i*ind)]%*%t(gradn[,((i-1)*ind+1):(i*ind)])}
  if((mind*ind)<n){
    hes=hes+gradn[,(mind*ind+1):n]%*%t(gradn[,(mind*ind+1):n])}
  res=list(solve(hes)%*%grad)
  return(res)
}
