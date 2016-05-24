na2Incomplete=function(from){
    d=as.integer(dim(from))
  rows=rep(1:d[1],d[2])
  cols=rep(1:d[2],rep(d[1],d[2]))
  nisna=!is.na(from)  
  new("Incomplete",as(sparseMatrix(i=rows[nisna],j=cols[nisna],x=from[nisna],dims=d),"dgCMatrix"))
  }
setAs("matrix","Incomplete",na2Incomplete)
sparse2Incomplete=function(from) new("Incomplete",as(from,"dgCMatrix"))
setAs("sparseMatrix","Incomplete",sparse2Incomplete)

 
Incomplete=function(i,j,x){
    new("Incomplete",as(sparseMatrix(i=i,j=j,x=x),"dgCMatrix"))
      }

  
