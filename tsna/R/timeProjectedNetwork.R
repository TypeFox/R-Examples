# create a crude version of discritized time projected graph from nD object

timeProjectedNetwork<-function(nd,start = NULL, end = NULL, time.increment = NULL, 
                               onsets = NULL, termini = NULL, ...){
  # slice up the dynamic network into a list of static networks
  nets<-get.networks(nd,start = start, end = end, time.increment = time.increment, 
                     onsets = onsets, termini = termini, 
                     retain.all.vertices=TRUE, # make sure to keep network size consistent so ids won't scramble
                     ...)
  if(is.bipartite(nd)){
    warning('input network is bipartite and projected network will not be')
  }
  if(is.hyper(nd)){
    stop('time projected network not supported for hypergraphic networks')
  }
  # init a new network of appropriate size
  projected<-network.initialize(length(nets)*network.size(nd),loops = has.loops(nd))
  # TODO: copy network-level attributes
  
  # loop over edges to create within-slice edges
  for(s in seq_along(nets)){
    el<-as.matrix.network.edgelist(nets[[s]])
    # multiply the vertex indices by s to give them for the next slice
    el<-el+(network.size(nd)*(s-1))
    if(nrow(el)>0){
      # create the edges and also assign the edge type attribute
      # and any other attributes when edges are created
      add.edges(projected,tail=el[,1],head=el[,2],
                names.eval = 'edge.type',
                vals.eval = 'within_slice')
      # figure out the eids of edges we just added to use in adding attributes later
      eids<-(projected$gal$mnext-nrow(el)):(projected$gal$mnext-1) 
      # handle un-directed edges by adding arc in each direction, with same attribute handling
      if(!is.directed(nd)){
        add.edges(projected,tail=el[,2],head=el[,1],
                  names.eval = 'edge.type',vals.eval = 'within_slice')
        eids<-c(eids,(projected$gal$mnext-nrow(el)):(projected$gal$mnext-1)) 
      }
      eattrnames<-list.edge.attributes(nets[[s]])
      evals<-lapply(eattrnames,function(attrname){get.edge.attribute(nets[[s]],attrname,unlist=FALSE)})
      # if edges were added in both directions, the vector of eids should be twice as long, need to copy values
      if(!is.directed(nd)){
        evals<-lapply(evals,function(vals){
          c(vals,vals)
        })
      }
      # to work around issue with set edge attribute being unable to distinguish when it should be setting multiple vs single attributes
      # need to unnest values if setting single
      if (length(eattrnames)==1){
        evals<-unlist(evals,recursive = FALSE)
      }
      set.edge.attribute(projected,eattrnames,evals,e = eids)
    }
    
  }
 
  # copy vertex attributes from slice networks
  for(attrname in list.vertex.attributes(nets[[1]])){
    # this gets all of the vertex attributes for each timepoint, and assigns them all at once 
    set.vertex.attribute(projected,attrname=attrname,
                         value=unlist(lapply(nets,get.vertex.attribute,attrname=attrname,unlist=FALSE),recursive=FALSE)
    )
  }
  
  # loop again to create between-slice vertex self ties
  # TODO: omit ties as indicated by vertex activity?
  if (length(nets)>1){
    for(s in 1:(length(nets)-1)){
      tails<-1:network.size(nd)+network.size(nd)*(s-1)
      heads<-1:network.size(nd)+network.size(nd)*s
      add.edges(projected,tail=tails, head=heads, 
                names.eval = 'edge.type',vals.eval = 'identity_arc' ) 
    }
  }
  return(projected)
}