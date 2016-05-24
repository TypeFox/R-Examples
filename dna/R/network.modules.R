network.modules=function(s,m,epsilon,plot=FALSE,interactive=FALSE,...){
 p=as.integer(nrow(s))
 if (is.null(colnames(s))){
  colnames(s)=paste("Gene",1:p)
  rownames(s)=paste("Gene",1:p)
 }
 out=.C("rgmd",as.double(s),module=integer(p),as.integer(m),as.double(epsilon),p,PACKAGE="dna") 
 out$module=as.factor(out$module)
 names(out$module)=rownames(s)
 if (plot==TRUE){
  require(igraph)
  if (sum(out$module!=0)>0){
   graph.genenames=names(out$module)[out$module!=0]
   graph.s=s[out$module!=0,out$module!=0]
   Is=abs(graph.s)>=epsilon
   rs=row(graph.s)
   cs=col(graph.s)
   vIs=c(Is)
   vrs=c(rs)
   vcs=c(cs)
   orc=rs<cs
   ex=vrs[vIs&orc]
   ey=vcs[vIs&orc]
   edges=rbind(ex,ey)
   g=graph.empty(directed=FALSE)
   g=add.vertices(g,length(graph.genenames),names=graph.genenames)
   g=add.edges(g,edges)
   if (interactive==TRUE)
    tkplot(g,vertex.label=V(g)$names,...)
   else
    plot(g,vertex.label=V(g)$names,...)
  }
  else{
   cat("No plot created since there are no modules in this network.\n")
  }
 }
 new("modules",module=out$module)
}

