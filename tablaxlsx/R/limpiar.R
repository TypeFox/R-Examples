limpiar <-
function(x,valores=NA,filas=TRUE,columnas=TRUE){
  if(is.null(dim(x)) | length(dim(x))!=2) return(NULL)
  if(is.null(valores)) return(NULL)
  if(min(dim(x))==0) return(NULL)	
  res1=as.matrix(apply(x,2,"%in%", valores))
  lfilas=!apply(res1,1,all)
  lcolumnas=!apply(res1,2,all)
  if(!filas) lfilas=rep(TRUE,nrow(x))
  if(!columnas) lcolumnas=rep(TRUE,ncol(x))
  return(list(Filas=lfilas,Columnas=lcolumnas))
}
