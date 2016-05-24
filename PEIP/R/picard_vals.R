picard_vals <-
function(U,sm,d)
  {
k=length(sm);
utd = vector(length=k)
utd_norm = vector(length=k)
for(i in 1:k)
  {
    utd[i] =t(U[,i]) %*%   d;
    utd_norm[i]=utd[i]/sm[i];
}

return(list(utd=utd,utd_norm=utd_norm) )
}
