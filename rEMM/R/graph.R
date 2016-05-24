
## igraph coercions
as.igraph <-  function(object) UseMethod("as.igraph")

as.igraph.TRACDS <-  function(object) smc_as.igraph(object@tracds_d$mm)

    
## graph coercions
as.graph <-  function(object) UseMethod("as.graph")

as.graph.TRACDS <-  function(object) {
    if(!.installed("graph")) stop ("Package graph needed! Please install from Bioconductor.")
    smc_as.graph(object@tracds_d$mm)
}


