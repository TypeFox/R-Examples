reducir <-
function(x,valores=NA,filas=TRUE,columnas=TRUE){
  a1=limpiar(x,valores,filas,columnas)
  nomfilas=rownames(x)
  nomcolumnas=colnames(x)
  nomdim=names(dimnames(x))
  x2=x[a1$Filas,a1$Columnas]
  at1=attributes(x)
  if(!is.null(attr(x,"cabColumna"))) attr(x2,"cabColumna")=attr(x,"cabColumna")
  if(!is.null(attr(x,"cabFila"))) attr(x2,"cabFila")=attr(x,"cabFila ")
  #x2=t(t(x2))
  if(!is.data.frame(x2)) dim(x2)=c(sum(a1$Filas),sum(a1$Columnas))
  if(!is.null(nomfilas)) rownames(x2)=nomfilas[a1$Filas]
  if(!is.null(nomcolumnas))colnames(x2)=nomcolumnas[a1$Columnas]
  if(!is.null(nomdim)) names(dimnames(x2))=nomdim
  if(!inherits(x,"ftable") & 
     (is.null(attr(x,"cabFila")) | is.null(attr(x,"cabColumna")))){
    return(x2)
  }  
  if(inherits(x,"ftable")){
    cabf=cabe(attr(x,"row.vars"))
    cabc=t(cabe(attr(x,"col.vars")))
  }else{
    cabf=attr(x,"cabFila")
    cabc=attr(x,"cabColumna")
  } 
  if(!is.null(cabf)){
    cabf1=cabf[a1$Filas,]
    if(is.null(dim(cabf1))) {
      cabf1=as.matrix(cabf1)
      colnames(cabf1)=colnames(cabf)
    }
    cabf2=quitarrep(cabf1)
    rs1=apply(cabf2,2,spanrow)
    attr(cabf1,"rowspan")=rs1
  }
  if(!is.null(cabc)){
    cabc1=cabc[,a1$Columnas]
    if(is.null(dim(cabc1))) {
      cabc1=t(as.matrix(cabc1))
      rownames(cabc1)=rownames(cabc)
    }
    cabc2=quitarrep(t(cabc1))
    rs2=t(apply(cabc2,2,spanrow))
    attr(cabc1,"colspan")=rs2  
  }
  attr(x2,"cabFila")=cabf1
  attr(x2,"cabColumna")=cabc1
  return(x2)
}
