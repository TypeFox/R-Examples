sym.var <-
function(sym.data,number.sym.var) {
  if((number.sym.var>sym.data$M)||(number.sym.var<=0)) 
    stop("number.sym.var out of range") 
  pos<-sym.data$sym.var.starts[number.sym.var]
  adv<-sym.data$sym.var.length[number.sym.var]
  sym.var<-list(N=sym.data$N,var.name=sym.data$sym.var.names[number.sym.var],
                var.type=sym.data$sym.var.types[number.sym.var],
                obj.names=sym.data$sym.obj.names,
                var.data.vector=sym.data$meta[,pos:(pos+adv-1)])
  return(sym.var)
}
