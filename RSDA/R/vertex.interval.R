vertex.interval <-
function(sym.data) {
  if ((sym.data$sym.var.types[1] != "$I")) {
    stop("Variables have to be continuos or Interval")
  }
  else{
    nn <- sym.data$N
    mm <- sym.data$M 
    num.vertex<-rep(-1,nn)
    vertex <- matrix(0, 1, mm)
    vertex <- as.data.frame(vertex)
    colnames(vertex) <- sym.data$sym.var.names
    sym.text<-"as.matrix(sym.data$data["
    for(i in 1:nn){
      current.row<-as.character(i)
      previous<-"1:2"
      command<-paste0(sym.text,current.row,",",previous,"])")
      for(j in 2:mm){
        col.current.min<-2*j-1
        col.current.max<-2*j
        nxt.grid<-paste0(as.character(col.current.min),":",as.character(col.current.max))
        command<-paste0(command,",",sym.text,current.row,",",nxt.grid,"])")
      }
      command<-paste0("expand.grid(",command,")")
      aux<-eval(parse(text=command))
      aux<-sqldf("select distinct * from aux")
      num.vertex[i]<-dim(aux)[1]
      colnames(aux) <- sym.data$sym.var.names
      vertex<-rbind(vertex,aux)
    }
    num.vertexf<-dim(vertex)[1]
    return(list(vertex = vertex[2:num.vertexf,], num.vertex = num.vertex))
  }
}
