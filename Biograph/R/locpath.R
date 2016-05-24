locpath <-
function(Bdata) 
  {loc <- which (colnames(Bdata)=="path",arr.ind=TRUE)
  return (loc)}
