
#given a correlations thingy as produced by brainwaver,
#produce a list of adjacency matrices fiddled to have a preferred number of edges
#Actually this is not quite possible, but come as close as choose.thresh.nbedges
#will allow us. A functional parameter allows us to say things like produce the graphs
#with n log n edges where n is the number of nodes

correlations.to.adjacencies<-function(correlations, edge.func)
{
  n.regions<-attr(correlations, "num.time.series")
  tmp<-attr(correlations,'names')
  n.levels<-length(tmp)/3  
  level.names<-tmp[1:n.levels]
  
  proc.length<-attr(correlations, "proc.length")

  
  if(exists("edge.func",mode="function")){
  favourite.edges<-edge.func(n.regions)
  }else{
  if(is.integer(edge.func)){
	favourite.edges<-edge.func
	}else{
	stop('Error the number of edges is not mentioned')
	}}

  adjmats=list()

  for(j in 1:n.levels){
    sup=choose.thresh.nbedges(correlations[[j]],nb.edges=favourite.edges, proc.length=proc.length, num.levels=j)
    adjacencies<-const.adj.mat(correlations[[j]], sup=sup, proc.length=proc.length, num.levels=j)
    cat(file=stderr(),"level", j, "threshold", sup, "edges",sum(adjacencies)/2,"\n")
    tmplist<-list(adjacencies)
    names(tmplist)<-level.names[[j]]
    adjmats<-c(adjmats,tmplist)
  }

  adjmats
}

#given a certain number of data points, what's the best number of wavelet levels
#to process?
ideal.wavelet.levels<-function(brain)
{
  number.of.samples<-length(brain[[1]])
  floor(log(number.of.samples,2))-5
}

#just your usual x^2+y^2+z^2
distance<-function(x,y,z)
{
  square<-function(x) x*x
  sqrt(square(x[2]-x[1])+square(y[2]-y[1])+square(z[2]-z[1]))
}

