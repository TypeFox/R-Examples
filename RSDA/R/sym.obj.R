sym.obj <-
function(sym.data,number.sym.obj) {
  if((number.sym.obj>sym.data$N)||(number.sym.obj<=0)) 
    stop("number.sym.obj out of range") 
  sym.obj<-list(M=sym.data$M,var.types=sym.data$sym.var.types,var.length=
                sym.data$sym.var.length,var.names=sym.data$sym.var.names,
                obj.data.vector=sym.data$data[number.sym.obj,])
  return(sym.obj)
}
