dm <-
function (x,y,diff=FALSE) { #produces dominance or difference matrix
n_x=length(x)
n_y=length(y)
#dx <- matrix(nrow=(n_x+1), ncol=(n_y+1))
dx <- matrix(nrow=n_x, ncol=n_y)
for (i in 1:n_x)
  {for (j in 1:n_y) 
  {dx[i,j]<--sign(y[j]-x[i])
   if(diff==TRUE){dx[i,j]<--c(y[j]-x[i])}}} 
#dx[,(n_y+1)] <- rowSums(dx)
#dx[(n_x+1),] <- colSums(dx)
#rownames(dx)<-c(x,"di.")
#colnames(dx)<-c(y,"d.j")
rownames(dx)<-x
colnames(dx)<-y
return(dx)}

