pairup=function(x,type="less") {
  x=as.matrix(x)
  n=dim(x)[1]
  i=rep(1:n,rep(n,n))
  j=rep(1:n,n)
  c1=apply(x,2,function(y){rep(y,rep(length(y),length(y)))})
  c2=apply(x,2,function(y){rep(y,length(y))})
  ans=cbind(c1,c2)
  ans=switch(type, less=ans[(i<j), ], leq=ans[i<=j, ], neq=ans)
  ans
}
