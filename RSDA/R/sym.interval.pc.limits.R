sym.interval.pc.limits <-
function(sym.data,prin.curve,num.vertex,lambda,var.ord) {
  num.vars<-sym.data$M 
  num.ind<-sym.data$N
  
  res<-as.data.frame(prin.curve)
  res$lambda<-lambda
  
  sym.indiv<-rep("X",sum(num.vertex))
  
  start<-1
  finish<- num.vertex[1]
  sym.indiv[start:finish]<-sym.data$sym.obj.names[1]
  
  for (i in 2:num.ind)
  {
    previous<-num.vertex[i-1]
    start<-start+previous
    finish<-num.vertex[i]+finish
    sym.indiv[start:finish]<-sym.data$sym.obj.names[i]
  }
  
  res$symindiv<-sym.indiv
  var.type<-rep('$I',num.vars+1)
  variables<-rep("X",num.vars)
  
  for(i in 1:num.vars){
    variables[var.ord[i]]<-paste0("prin_surface_" , as.character(i))
  }
  colnames(res)[1:num.vars]<-variables
  variables<-c(variables[var.ord],"lambda")
  
  sym.res<-classic.to.sym(dataTable = res,
                          concept = c("symindiv"),
                          variables = variables,
                          variables.types = var.type)
  return(sym.res)
}
