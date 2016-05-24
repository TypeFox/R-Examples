plotL.graph <-
function(rl, rlist,nr,species,links,...)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    element <- extract.graph(rl,rlist,nr)
    element1 <- cluster.id(element)
    plot_graph(element1,species=species,links)
  }
