suv <-
function(u,v,irow,jcol){
###  computes (u%*%t(v))[i,j]= sum(u[i,]*v[j,]) for all pairs in irow,jcol
  dd=dim(u)
  nnrow=as.integer(dd[1])
  nncol=as.integer(nrow(v))
  nrank=dd[2]
  storage.mode(u)="double"
  storage.mode(v)="double"
  storage.mode(irow)="integer"
  storage.mode(jcol)="integer"
  nomega=as.integer(length(jcol))
  .Fortran("suv",
           nnrow,nncol,nrank,u,v,irow,jcol,nomega,
           r=double(nomega),
           PACKAGE="softImpute"
           )$r
}

