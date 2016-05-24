#  File networkDynamic/R/utilities.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012,2013 the statnet development team

################################################
# utilities.R
# Author: Zack W. Almquist, Pavel , Li Wang lxwang@uw.edu, skyebend@uw.edu
#
#
################################################



# see https://statnet.csde.washington.edu/trac/wiki/NetworkDynamicConverterFunctions
# for full specs on this constructor function

networkDynamic <- function(base.net=NULL,edge.toggles=NULL,vertex.toggles=NULL,
                  edge.spells=NULL,vertex.spells=NULL,edge.changes=NULL,vertex.changes=NULL,
                  network.list=NULL,onsets=NULL,termini=NULL,vertex.pid=NULL,start=NULL,end=NULL,net.obs.period=NULL,verbose=TRUE,create.TEAs=FALSE,edge.TEA.names=NULL,vertex.TEA.names=NULL,...) {
  
  
  if (!is.null(start) && !is.null(end)) {
    if (start > end) stop ("start must be less than end")
  }
  if (!is.null(start) || !is.null(end)) {
    if (!is.null(net.obs.period)) stop("net.obs.period can not be specified with start and end arguments")
  }
  if (!is.null(base.net)) {
    if (!is.network(base.net)) stop("base.net must be either NULL or a network object")
  }
  
  if (is.null(network.list) & (!is.null(onsets) || !is.null(termini))) {
    warning("Onsets and termini arguments only used when a network.list is provided.")
  }
  
  # ------- network.list specified ---------
  if (!is.null(network.list)) {
    # ---- check parameters ----
    if (!is.null(vertex.toggles) || !is.null(vertex.spells) ||!is.null(vertex.changes) ||
      !is.null(edge.toggles) || !is.null(edge.spells) ||!is.null(edge.changes)) {
      stop("Vertex or edge activation arguments can not be used along with a network.list")
    }
    
    if (!is.null(onsets)) {
      if (!(is.numeric(onsets) && (length(onsets) == length(network.list)))) {
        stop("Onsets and termini must be NULL or numeric, and the same length as network.list")
      }
    }
    if (!is.null(termini)) {
      if (!(is.numeric(termini) && (length(termini) == length(network.list)))) {
        stop("Onsets and termini must be NULL or numeric, and the same length as network.list")
      }
    }
    
    net.sizes = sapply(network.list, function(y){if (is.network(y)) network.size(y) else NA})
    if (sum(sapply(net.sizes, is.na)) > 0) {
      stop("All elements of network.list must be network objects")
    }
    
    # check network attributes in case they are not the same on all networks
    # this only checks non-user attributes
    net.attribs <- net.attribs <- c('directed', 'hyper', 'loops', 'multiple', 'bipartite')
    temp <- do.call(rbind, lapply(network.list, 
                                 function(x) {
                                    sapply(net.attribs, function(a) get.network.attribute(x, a))}))
    if (!all(sapply(seq_len(ncol(temp)), function(x) length(unique(temp[,x])))==1))
      warning("Networks in network.list have different network properties. Only the first network's attributes are retained.")
    
    if (length(unique(net.sizes)) != 1) {
      if (is.null(vertex.pid)) {
        stop("vertex.pid must be specified when network sizes of network.list are different")
      }
    }
    
    if (!is.null(vertex.pid)) {
      if (vertex.pid %in% list.vertex.attributes(network.list[[1]])) {
        # placeholder
      } else {
        stop("vertex.pid must be present in the network vertex attributes")
      }
    }
    
    # observation period
    if (!is.null(start) && !is.null(onsets)) {
      stop("only one of start and onsets should be specified")
    }
    
    if (is.null(net.obs.period)) {
      # end is not used. if start is specified, assume end to be start + length
      
      if (is.null(start) && is.null(onsets)) {
        start <- 0;
        if (verbose){
          cat("Neither start or onsets specified, assuming start=0\n")
        }
      }
      # only start
      if (is.null(onsets)) {
        net.obs.period <- list(observations=list(c(start,start+length(network.list))),mode="discrete", time.increment=1,time.unit="step")
        onsets<-seq(from=start, length=length(network.list))
        termini<-seq(from=start, length=length(network.list)) + 1
        if (verbose){
          cat("Onsets and termini not specified, assuming each network in network.list should have a discrete spell of length 1\n")
        }
      } else {
        # only onsets
        if (is.null(termini)) stop("onsets and termini must be specified together")
        obs <- lapply(1:length(termini), function(i) c(onsets[i], termini[i]))
        net.obs.period <- list(observations=obs, mode='continuous',time.increment=NA,time.unit='unknown')
      }
        
    } 
    copyTEAs<-create.TEAs #flag to control if we process TEAs, which is slow
    if (copyTEAs){
      # construct the list of vertex attributes that may become TEAs
      TEAvertAttrs<-unique(unlist(lapply(network.list,list.vertex.attributes)))
      TEAvertAttrs<-TEAvertAttrs[!TEAvertAttrs%in%c('na','vertex.names',vertex.pid)]
      # also for edges
      TEAedgeAttrs<-unique(unlist(lapply(network.list,list.edge.attributes)))
      TEAedgeAttrs<-TEAedgeAttrs[!TEAedgeAttrs%in%c('na')]
      # and for networks
      TEAnetAttrs<-unique(unlist(lapply(network.list,list.network.attributes)))
      TEAnetAttrs<-TEAnetAttrs[!TEAnetAttrs%in%c('directed', 'hyper', 'loops', 'multiple', 'bipartite','mnext','net.obs.period','vertex.pid','edge.pid','n')]
      if (verbose){
        cat("Dynamic attributes (TEAs) will be created for the following attributes of the networks on network.list:\n")
        if (length(TEAnetAttrs)){
          cat("\tNetwork: ",TEAnetAttrs,"\n")
        }
        if (length(TEAvertAttrs)){
          cat("\tVertices: ",TEAvertAttrs,"\n")
        }
        if (length(TEAedgeAttrs)){
          cat("\tEdges: ",TEAedgeAttrs,"\n")
        }
      }
    }
    
    # ---- vertex.pid present ----
    if (!is.null(vertex.pid)) {
      # build base network vertices <-> pids
      base.net.pids <- NULL      
      for (i in seq_along(network.list)) {
        if (vertex.pid %in% list.vertex.attributes(network.list[[i]])) {
          net.pids <- get.vertex.attribute(network.list[[i]], vertex.pid)
          if (!is.unique.list(net.pids)) 
            stop("vertex.pid attribute must be unique for each vertex.")
          base.net.pids <- c(base.net.pids, net.pids)
        } else {
          stop("vertex.pid must be present in the network vertex attributes for each network")
        }
      }  
      base.net.pids <- sort(unique(base.net.pids))
      
      # initialize, copy network attributes
      vattrs=character(0)
      nattrs=character(0)
      if (is.null(base.net)){ 
        base.net <- network.list[[1]]
        if (verbose){
          cat("Argument base.net not specified, using first element of network.list instead\n")
        }
        # attributes will not be copied from base net
      } else {
        # get the list of attributes to copy from 
        vattrs<-list.vertex.attributes(base.net)
        if (copyTEAs){
          # don't create a TEA if attribute of the same name exists only in base net
          TEAvertAttrs<-setdiff(TEAvertAttrs,vattrs)
        }
        
        nattrs<-list.network.attributes(base.net)
        nattrs<-setdiff(nattrs,c("bipartite","directed","hyper","loops","mnext","multiple","n"))
      }
      # todo: should we warn here that used first network as base net?
      # I guess the reason to init instead of copy is in case the pids imply a different network size?
      out.net <- network.initialize(length(base.net.pids), directed = base.net%n%"directed", 
                                   hyper = base.net%n%"hyper", loops = base.net%n%"loops", 
                                   multiple = base.net%n%"multiple", bipartite = base.net%n%"bipartite")
      if (verbose){
        cat(paste("Initialized network of size",network.size(out.net),"inferred from number of unique vertex.pids\n"))
      }
      # copy vertex attributes (only occurs if real base net specified)
      
      for(attr in vattrs){
        set.vertex.attribute(out.net,attr,get.vertex.attribute(base.net,attr))
      }
      # copy network attributes EXCEPT network properties (if real base.net specified)
      for (attr in nattrs){
        set.network.attribute(out.net,attr,get.network.attribute(base.net,attr))
      }
    
      # store the vertex.pids
      set.vertex.attribute(out.net, vertex.pid, base.net.pids)
      # store the name of the vertex.pid for later use
      set.network.attribute(out.net,'vertex.pid',vertex.pid)
      
      base.net <- out.net
        
      # get combined edge list, indexed by the vertices in base.net
      for (i in seq_along(network.list)) {
        edgelist <- as.edgelist(network.list[[i]])
        net.pids <- get.vertex.attribute(network.list[[i]], vertex.pid)
        # convert the network vertex indices to base.net indices, using pids
        edges <- apply(edgelist, c(1,2), function(x) {
          which(base.net.pids == net.pids[x])[1]
        })
        
        # activate the vertices
        vertices<-match(net.pids, base.net.pids)
        os<-rep(onsets[i], length(vertices))
        ts<-rep(termini[i], length(vertices))
        activate.vertices(base.net,onset=os,terminus=ts,v=vertices)
        
        if (copyTEAs){
          # copy any non-standard, vertex, attributes into TEAs
          for(attr in TEAvertAttrs){
            activate.vertex.attribute(base.net,attr,get.vertex.attribute(network.list[[i]],attr,unlist=FALSE),onset=onsets[i],terminus=termini[i],v=get.vertex.id(base.net,net.pids))
          }
        }
        
        # activate the edges
        for (e in seq_len(nrow(edges))) {
          t <- edges[e,1]
          h <- edges[e,2]
          # add edge if necessary
          if (length(get.edgeIDs(base.net, t, h)) == 0) {
            add.edge(base.net, tail=t, head=h)
          }
          eid<-get.edgeIDs(base.net,t,h)[1]
          activate.edges(base.net, e = eid,  onset=onsets[i], terminus=termini[i])
          if (copyTEAs){
            #  copy any non-standard edge attributes into TEAs
            for(attr in TEAedgeAttrs){
              activate.edge.attribute(base.net,attr,get.edge.attribute(network.list[[i]],attrname=attr,unlist=FALSE)[get.edgeIDs(network.list[[i]],v=edgelist[e,1],alter=edgelist[e,2])],e=eid,onset=onsets[i],terminus=termini[i])
            }
          }
          
        } # end edge copy loop
        
        if (copyTEAs){
          # copy any non-standard network attributes into TEAs
          for(attr in TEAnetAttrs){
            activate.network.attribute(base.net,attr,get.network.attribute(network.list[[i]],attr,unlist=FALSE),onset=onsets[i],terminus=termini[i])
          }
        }
        
      } # end network loop
      
      
      
      
      
    } else {
      # ---- no vertex.pid, all networks same size ----
      
      # initialize, copy network attributes
      if (is.null(base.net)) base.net <- network.list[[1]]
      if (verbose){
        cat("Argument base.net not specified, using first element of network.list instead\n")
      }
      out.net <- network.initialize(network.size(base.net), directed = base.net%n%"directed", 
                                   hyper = base.net%n%"hyper", loops = base.net%n%"loops", 
                                   multiple = base.net%n%"multiple", bipartite = base.net%n%"bipartite")
      # copy vertex attributes  These may be overwritten by TEAs if they appear on later networks
      vattrs<-list.vertex.attributes(base.net)
      for(attr in vattrs){
        set.vertex.attribute(out.net,attr,get.vertex.attribute(base.net,attr))
      }
      # copy network attributes
      nattrs<-list.network.attributes(base.net)
      nattrs<-setdiff(nattrs,c("bipartite","directed","hyper","loops","mnext","multiple","n"))
      for (attr in nattrs){
        set.network.attribute(out.net,attr,get.network.attribute(base.net,attr))
      }
      base.net <- out.net
      
      # get combined edge list, indexed by the vertices in base.net
      for (i in seq_along(network.list)) {
        edges <- as.edgelist(network.list[[i]])
        
        # all vertices are assumed to be active
        activate.vertices(base.net, onset=onsets[i], terminus=termini[i])
        
        if(copyTEAs){
          # copy any vertex attributes into a vertex TEA
          for(attr in TEAvertAttrs){
            activate.vertex.attribute(base.net,attr,get.vertex.attribute(network.list[[i]],attr,unlist=FALSE),onset=onsets[i],terminus=termini[i])
          }
        }
        
        # activate the edges
        for (e in seq_len(nrow(edges))) {
          t <- edges[e,1]
          h <- edges[e,2]
          # add edge if necessary
          if (length(get.edgeIDs(base.net, t, h)) == 0) {
            add.edge(base.net, tail=t, head=h)
          }
          eid <- get.edgeIDs(base.net,t,h)[1]
          activate.edges(base.net, e = eid, onset=onsets[i], terminus=termini[i]) 
          if (copyTEAs){
            # copy any edge attributes into an edge TEA
            for(attr in TEAedgeAttrs){
              activate.edge.attribute(base.net,attr,get.edge.attribute(network.list[[i]],attrname=attr,unlist=FALSE)[get.edgeIDs(network.list[[i]],v=t,alter=h)],e=eid,onset=onsets[i],terminus=termini[i])
            }
          }
        }
        if (copyTEAs){
          # copy any non-standard network attributes into TEAs
          for(attr in TEAnetAttrs){
            activate.network.attribute(base.net,attr,get.network.attribute(network.list[[i]],attr,unlist=FALSE),onset=onsets[i],terminus=termini[i])
          }
        }
        
      }
      
      
      
    }
    ## check that net.obs.period has appropriate structure
    .check.net.obs.period(net.obs.period)
    if (verbose){
      cat("Created net.obs.period to describe network\n")
      cat(.print.net.obs.period(net.obs.period))
    }
    set.network.attribute(base.net, "net.obs.period", net.obs.period)
    
  } else {  # end network list block
    #todo: if vertex.pid is not null, should we store it? (as a way of setting vertex pid)
    
    # ---------- edge or vertex timings specified  -----------
    # ---- check parameters ----
    if (!is.null(network.list)) {
      stop("Vertex or edge activation arguments can not be used along with a network.list")
    }
    vertex.args <- list(vertex.toggles, vertex.changes, vertex.spells)
    vertex.which <- which(!sapply(vertex.args, is.null))
    if (length(vertex.which) > 1) 
      stop("Only one of vertex.toggles, vertex.spells and vertex.changes should be specified.")
    # pick out the non-null argument, if it is given
    vertex.data <- (if (length(vertex.which)==1) vertex.args[[vertex.which]] else NULL)
  
    
    if (!is.null(vertex.data)){
      # try to convert it if it is not a matrix or a data.frame
      # (don't want to convert data.frame to matrix if we can help it because may force numerics into chars)
      if(!is.data.frame(vertex.data) & !is.matrix(vertex.data)){
        vertex.data <- as.matrix(vertex.data)
      }
      
      if (!is.null(vertex.changes)) {
        if (ncol(vertex.data) < 3) stop("vertex.changes requires 3 columns: time, vertex.id, direction")
        cnames<-colnames(vertex.data,do.NULL=FALSE) # avoid problem when no names exist
        cnames[1:3]<- c("time", "vertex.id", "direction")
        colnames(vertex.data)<-cnames
        if(!is.numeric(vertex.data[,'time'])){
          stop("vertex.changes requires the time column to be numeric")
        }
        if(!is.numeric(vertex.data[,'vertex.id'])){
          stop("vertex.changes requires the vertex.id column to be numeric")
        }
        if(!is.numeric(vertex.data[,'direction'])){
          stop("vertex.changes requires the direction column to be numeric")
        }
      }
      if (!is.null(vertex.toggles) ) {
        if (ncol(vertex.data) < 2) stop("vertex.toggles requires 2 columns: time, vertex.id")
        cnames<-colnames(vertex.data,do.NULL=FALSE)
        cnames[1:2] <- c("time", 'vertex.id')
        colnames(vertex.data)<-cnames
        if(!is.numeric(vertex.data[,'time'])){
          stop("vertex.toggles requires the time column to be numeric")
        }
        if(!is.numeric(vertex.data[,'vertex.id'])){
          stop("vertex.toggles requires the vertex.id column to be numeric")
        }
      }
      
      if (!is.null(vertex.spells)) {
        if (ncol(vertex.data) < 3) stop("vertex.spells requires 3 columns: onset, terminus, vertex.id") 
        cnames<-colnames(vertex.data,do.NULL=FALSE)
        cnames[1:3] <- c('onset', 'terminus', 'vertex.id')
        colnames(vertex.data)<-cnames
        
        if(!is.numeric(vertex.data[,'vertex.id'])){
          stop("vertex.spells requires the vertex.id column to be numeric")
        }
        if(!is.numeric(vertex.data[,'onset'])){
          stop("vertex.spells requires the onset time column to be numeric")
        }
        if(!is.numeric(vertex.data[,'terminus'])){
          stop("vertex.spells requires the terminus time column to be numeric")
        }
      }
      
    }
        
    edge.args <- list(edge.toggles, edge.changes, edge.spells)
    edge.which <- which(!sapply(edge.args, is.null))
    if (length(edge.which) > 1)
      stop("Only one of edge.toggles, edge.spells and edge.changes should be specified.")
    edge.data <- (if (length(edge.which)==1) edge.args[[edge.which]] else NULL)
    
    if (!is.null(edge.data)) {
      
      if(!is.data.frame(edge.data) & !is.matrix(edge.data)){
        # don't want to cast to matrix if it is data frame 'cause may mangle numeric columns
        edge.data <- as.matrix(edge.data)
      }
      if (!is.null(edge.changes)) {
        if (ncol(edge.data) < 4) stop("edge.changes requires 4 columns: time, tail, head, direction")
        cnames<-colnames(edge.data,do.NULL=FALSE)
        cnames[1:4] <- c('time', 'tail', 'head', 'direction')
        colnames(edge.data)<-cnames
        if (!is.numeric(edge.data[,'time'])){
          stop('the time column of the edge.changes argument to networkDynamic must be numeric')
        }
        if (!is.numeric(edge.data[,'tail'])){
          stop('the tail column of the edge.changes argument to networkDynamic must be a numeric vertex id')
        }
        if (!is.numeric(edge.data[,'head'])){
          stop('the head column of the edge.changes argument to networkDynamic must be a numeric vertex id')
        }
        if (!is.numeric(edge.data[,'direction'])){
          stop('the direction column of the edge.changes argument to networkDynamic must be numeric')
        }
      }
      if (!is.null(edge.toggles)) {
        if (ncol(edge.toggles) < 3) stop("edge.toggles requires 3 columns: time, tail, head")
        cnames<-colnames(edge.data,do.NULL=FALSE)
        cnames[1:3] <- c('time', 'tail', 'head')
        colnames(edge.data)<-cnames
        if (!is.numeric(edge.data[,'time'])){
          stop('the time column of the edge.toggles argument to networkDynamic must be numeric')
        }
        if (!is.numeric(edge.data[,'tail'])){
          stop('the tail column of the edge.toggles argument to networkDynamic must be a numeric vertex id')
        }
        if (!is.numeric(edge.data[,'head'])){
          stop('the head column of the edge.toggles argument to networkDynamic must be a numeric vertex id')
        }
      }
      if (!is.null(edge.spells)) {
        if (ncol(edge.spells) < 4) stop("edge.spells requires 4 columns: onset, terminus, tail, head")
        cnames<-colnames(edge.data,do.NULL=FALSE)
        cnames[1:4] <- c('onset', 'terminus', 'tail', 'head')
        colnames(edge.data)<-cnames
        if (!is.numeric(edge.data[,'onset'])){
          stop('the onset time column of the edge.spells argument to networkDynamic must be numeric')
        }
        if (!is.numeric(edge.data[,'terminus'])){
          stop('the terminus time column of the edge.spells argument to networkDynamic must be numeric')
        }
        if (!is.numeric(edge.data[,'tail'])){
          stop('the tail column of the edge.spells argument to networkDynamic must be a numeric vertex id')
        }
        if (!is.numeric(edge.data[,'head'])){
          stop('the head column of the edge.spells argument to networkDynamic must be a numeric vertex id')
        }
      }
    } # end edge data preformatting
      
  
    
    if (is.null(edge.data) && is.null(vertex.data)) warning('neither edge or vertex data were included for network construction')
    
    # ---- initialize base.net ----
    # fill in base network if it is not given
    max.vertex <- max(vertex.data[,'vertex.id'], edge.data[,'tail'], edge.data[,'head'],0)
    if (is.null(base.net)){
      if (verbose){
        cat(paste("Initializing base.net of size",max.vertex,"imputed from maximum vertex id in edge records\n"))
      }
      base.net <- network.initialize(max.vertex)
    } 
    
    
    # strict construction for now
    if (max.vertex > network.size(base.net)) stop("base.net network size is smaller than size implied by vertex.ids in vertex or edge argument")
    
    # remove any activity from base.net (for now)
    delete.vertex.activity(base.net)
    if (network.edgecount(base.net) > 0){
      delete.edge.activity(base.net)
      cat("Edge activity in base.net was ignored\n")
    }
    
    
    # ---- vertex data ----
    if (!is.null(vertex.data)){
      # sort by time
      vertex.data <- vertex.data[order(vertex.data[,1,drop=FALSE]), ,drop=FALSE]
      
      # initialize
      if (!is.null(vertex.toggles)) activate.vertices(base.net, onset=-Inf, terminus=Inf)
      
      if (!is.null(vertex.spells)) { # doing vertex spell data
        activate.vertices(base.net, v=vertex.data[,'vertex.id'], onset=vertex.data[,'onset'], terminus=vertex.data[,'terminus'])
        
        # vertex TEA stuff
        if (create.TEAs){
          
          # if vertex.TEA.names is missing, try to read in col names from vertex.spells
          if (is.null(vertex.TEA.names)){
            vertex.TEA.names<-colnames(vertex.spells)
            if(!is.null(vertex.TEA.names)){
              # remove the first four names because they are onset, terminus, etc.
              vertex.TEA.names<-tail(vertex.TEA.names,-3)
            }
          }
          
          # check that the length of vertex TEA names matches extra columns given
          if (ncol(vertex.spells)-3 != length(vertex.TEA.names)){
            stop('the vector of vertex.TEA.names must match the number of remaining columns in vertex.spells')
          }
          vids<-vertex.spells[,3] # the ids of the vertices in each row
          uniqueVids<-unique(vids)
          # loop for each attribute
          for (attrIndex in seq_along(vertex.TEA.names)){
            
            # construct the edge TEA directly using list operations
            # because the API methods are orders of magnitude too slow
            # NOTE: this makes the assumption that values always change
            # in other words, adjacent spells with the same value will not be merged,
            # where they would be using the activate. methods
            toBeTEA<-lapply(uniqueVids,function(vid){
              # find index of all spells that need to associated to this id
              rows<-which(vids==vid)
              # set sort order
              rows<-rows[order(vertex.spells[rows,1])]
              # construct values list
              vals<-as.list(vertex.spells[rows,3+attrIndex])
              spls<-as.matrix(vertex.spells[rows,1:2,drop=FALSE])
              dimnames(spls)<-NULL
              # veryify spell matrix consistency
              if(!all(chkspellmat(spls))){
                warning("vertex spell data induced an invalid spell matrix for vertex TEA '",vertex.TEA.names[attrIndex],"' for vertex id ", vid)
              }
              
              list(vals,spls)
            })
            set.vertex.attribute(base.net,attrname=paste(vertex.TEA.names[attrIndex],'active',sep='.'),value=toBeTEA, v=uniqueVids )
          }
          
          # debug
          if (verbose & length(vertex.TEA.names)>0){
            cat('Activated TEA vertex attributes: ',paste(vertex.TEA.names,collapse=', '))
          }
        }
        
      } else {
        for (i in seq_len(nrow(vertex.data))) {
          # todo: maybe do this in try catch so we can give appropriate line numbers for errors?
          at <- vertex.data[i,'time']
          v <- vertex.data[i,'vertex.id']
          change.activate <- (if (is.null(vertex.changes)) !is.active(base.net, at=at, v=v) else vertex.data[i,'direction']==1) 
          if (change.activate) {
            activate.vertices(base.net, v=v, onset=at, terminus=Inf)
          } else {
            deactivate.vertices(base.net, v=v, onset=at, terminus=Inf)
          }
        }
      }
    }
    
    # ---- edge data ----
    if (!is.null(edge.data)) {        
      # sort by onset time
      edge.data <- edge.data[order(edge.data[,1]), ,drop=FALSE]
      
      # initialize
      #if (is.null(edge.spells)) activate.edges(base.net, onset=-Inf, terminus=Inf)
      
      #if we are in the spells case, 
      if (!is.null(edge.spells) ){
        # if no edges exist yet we can avoid actually looping
        if (network.edgecount(base.net) ==0){
          dyads<-unique(edge.data[,3:4,drop=FALSE])
          tails<-as.list(dyads[,1])
          heads<-as.list(dyads[,2])
          add.edges(base.net,tail=tails,head=heads)
          eids<-sapply(seq_len(nrow(edge.data)),function(i){get.edgeIDs(base.net,v=edge.data[i,3],alter=edge.data[i,4])})
          if (length(eids)>0){
            activate.edges(base.net,onset=edge.data[,1],terminus=edge.data[,2],e=eids)
          }
        } else {  # there are pre-existing edges, so need to check for them while looping
          for (i in seq_len(nrow(edge.data))){
            eid<-get.edgeIDs(base.net,v=edge.data[i,3],alter=edge.data[i,4])
            if (length(eid)==0){
              # if the edge doesn't exist, add it and get the new edge id
              add.edge(base.net,tail=edge.data[i,3],head=edge.data[i,4])
              eid<-get.edgeIDs(base.net,v=edge.data[i,3],alter=edge.data[i,4])
            } 
            activate.edges(base.net,onset=edge.data[i,1],terminus=edge.data[i,2],e=eid)
          }
        }
  
      } else {  # we are processing changes or toggles
        
        # if no edges exist in base.net, can at least pre-create the edges
        if (network.edgecount(base.net) ==0){
          dyads<-unique(edge.data[,2:3,drop=FALSE])
          tails<-as.list(dyads[,1])
          heads<-as.list(dyads[,2])
          add.edges(base.net,tail=tails,head=heads)
          
          
          # try to pre-generate the activity matrix for each edge
          # and set it directly, bypasing the activity methods
		      spls<-lapply(seq_len(nrow(dyads)),function(e){
            
		        spellsFromChanges(edge.data,dyads[e,1],dyads[e,2],strict=FALSE)
            # strict=FALSE means it will ignore any out-of-order changes
		      })
		      base.net<-set.edge.attribute(base.net,'active',spls)
          
          
        } else { # some pre-existing edges, so have to loop to avoid hurting our heads
           #TODO: this could be optimized much more, see version in tergm
          # if there are pre-existing edges, assume edges in base.net to be active initially
          if (!is.null(edge.toggles)) activate.edges(base.net, onset=-Inf, terminus=Inf)
        
          for (i in seq_len(nrow(edge.data))) {
            t <- edge.data[i,2] 
            h <- edge.data[i,3]
            e <- get.edgeIDs(base.net, t, h)
            # add edge if not present in the base.net (as inactive?)
            # TODO: problem here with directed vs undirected networks?
            # TODO: how to handle multiplex/duplicate edge case?
            if (length(e) == 0) {
              add.edge(base.net, t, h)
              e <- get.edgeIDs(base.net, t, h)
              if (!is.null(edge.toggles)) deactivate.edges(base.net, e=e, onset=-Inf, terminus=Inf)
            }
            at <- edge.data[i,'time']
            change.activate <- (if (!is.null(edge.toggles)) !is.active(base.net, at=at, e=e) else edge.data[i,'direction']==1) 
            if (change.activate) {
              activate.edges(base.net, e=e, onset=at, terminus=Inf)
            } else {
              deactivate.edges(base.net, e=e, onset=at, terminus=Inf)
            }
            
          }
        }
      } # end of non-spell edge creation
      
      # begin edge TEA stuff
      # only do TEAs for edge spells
      if (!is.null(edge.spells) & create.TEAs){
        
        # if edge.TEA.names is missing, try to read in col names from edge.spells
        if (is.null(edge.TEA.names)){
          edge.TEA.names<-colnames(edge.spells)
          if(!is.null(edge.TEA.names)){
            # remove the first four names because they are onset, terminus, etc.
            edge.TEA.names<-tail(edge.TEA.names,-4)
          }
        }
        
        # check that the length of edge TEA names matches extra columns given
        if (ncol(edge.spells)-4 != length(edge.TEA.names)){
          stop('the vector of edge.TEA.names must match the number of remaining columns in edge.spells')
        }
        # get the vector of eids corresponding to the edges that have been created
        eids<-get.dyads.eids(base.net,tails=edge.spells[,3],heads=edge.spells[,4])
        uniqueEids<-unique(eids)
        # loop for each attribute
        for (attrIndex in seq_along(edge.TEA.names)){
          # construct the edge TEA directly using list operations because the API methods 
          # are orders of magnitude slower for this use case
          # NOTE: this makes the assumption that values always change
          # in other words, exactly adjacent spells with the same value will not be merged,
          # where they would be using the activate. methods
          toBeTEA<-lapply(uniqueEids,function(eid){
            # grab all of the spells that need to associated to this id
            rows<-which(eids==eid)
            # set sort order for spells
            rows<-rows[order(edge.spells[rows,1])]
            # construct values list
            vals<-as.list(edge.spells[rows,4+attrIndex])
            spls<-as.matrix(edge.spells[rows,1:2,drop=FALSE])
            dimnames(spls)<-NULL
            # veryify spell matrix consistency
            if(!all(chkspellmat(spls))){
              warning("edge spell data induced an invalid spell matrix for edge TEA '",edge.TEA.names[attrIndex],"' for edge id ", eid)
            }
            list(vals,spls)
          })
          set.edge.attribute(base.net,attrname=paste(edge.TEA.names[attrIndex],'active',sep='.'),value=toBeTEA,e =uniqueEids )
          
        }
        
        # debug
        if (verbose & length(edge.TEA.names)>0){
          cat('Activated TEA edge attributes: ',paste(edge.TEA.names,collapse=', '))
        }
      }
      
    } # end edge data
    
    # if the net.obs.period was not provided, create one
    if (is.null(net.obs.period)) {
      # observation.period and censoring
      # default to max and min of time range
      if(is.null(start) | is.null(end)){
        timeRange<-range(get.change.times(base.net,ignore.inf=FALSE))
      }
      if (is.null(start)){
        start <- timeRange[1]
      }
      if (is.null(end)){
        end <- timeRange[2]
      }
      
      # if we are in the toggles case, need
      net.obs.period <- list(observations=list(c(start, end)))
      if (!is.null(edge.spells) || !is.null(vertex.spells)) {
        
        net.obs.period$mode <- 'continuous'
        net.obs.period$time.increment<-NA
        net.obs.period$time.unit<-"unknown"
      } else {
        net.obs.period$mode <- 'discrete'
        net.obs.period$time.increment<-1
        net.obs.period$time.unit<-"step"
      }
    } else {
      if (!is.null(start)) stop("start and end arguments should not be specified with net.obs.period argument")
    }
    # verify that net.obs.period has good structure
    .check.net.obs.period(net.obs.period)
    if (verbose){
      cat("Created net.obs.period to describe network\n")
      cat(.print.net.obs.period(net.obs.period))
    }
    set.network.attribute(base.net, "net.obs.period", net.obs.period)
    
  } # end non-network.list part
  
  # if only base net is specified, set.nD.class on it and return. 
  return(set.nD.class(base.net))
  
}

# draft of internal function for retriving batch set of eids associated with dyads

get.dyads.eids<-function(net, tails, heads){
  if(length(tails)!=length(heads)){
    stop('the length of the tails and heads parameters must be the same to define the set of dyads to check')
  }
  sapply(seq_len(length(tails)),function(e){
    eids<-get.edgeIDs(net,v=tails[e],heads[e])
    if (length(eids)>1){
      warning('get.dyads.eids found multiple edges for given dyad (multiplex network), only smallest eid returned')
      eids<-min(eids)
    }
    if (length(eids)<1){
      eids<-NA
    }
    return(eids)
    })
}

# Get activity functions

# wrapper functions to return activity matrices of edges and vertices
get.edge.activity <- function(x, e=seq_along(x$mel), as.spellList=FALSE,active.default=TRUE) {
  if(length(x$mel)>0) 
    if((min(e) < 1) || (max(e) > x%n%"mnext"-1)) 
      stop("Illegal edge in get.edge.activity.\n")
  
  if (as.spellList) {
    return(as.data.frame.networkDynamic(x, e=e,active.default=active.default,start=-Inf,end=Inf))
  } else {
    spls<-get.edge.attribute(x$mel, "active", unlist=FALSE)
    if (active.default & length(e)>0){
      # hard to distinguish between edges that are missing one ones with no activity
      # figure out which e correspond to deleted edges and get rid of 'em
      eExist<-which(!sapply(x$mel,is.null))
      # replace spells that exist and are null
      spls[intersect(eExist,which(sapply(spls,is.null)))]<-list(matrix(c(-Inf,Inf),ncol=2))
    }
    # also remove any 'null spells' (Inf,Inf)
    spls<-.removeNullSpells(spls)
    #trim to e
    spls<-spls[e]
    return(spls)
  }
  
}

get.vertex.activity <- function(x, v=seq_len(network.size(x)), as.spellList=FALSE,active.default=TRUE) {
  if((min(v,Inf) < 1) || (max(v,-Inf) > network.size(x))) 
    stop("Illegal vertex in get.vertex.activity.\n")  
  
  if (as.spellList) {
    return(get.vertex.spelllist(x, v=v,active.default=active.default))
  } else {
    vam=get.vertex.attribute(x, "active", unlist=FALSE)
    vam<-vam[v]
    if(length(vam)>0){
      # if active default, need to add spells for vertices with no spells
      if (active.default ){
        vam[is.na(vam)]<-list(matrix(c(-Inf,Inf),ncol=2))
      }
      # also 'remove any' 'null spells' (Inf,Inf)
      # don't actually remove them, because then would lose order of v
      vam<-.removeNullSpells(vam)
    } else { # no vertices case
      vam<-list()
    }
  }
  return(vam)
}



# tail and head are nodeIDs
as.data.frame.networkDynamic<-function(x, row.names = NULL, optional = FALSE,e=seq_along(x$mel), start=NULL, end=NULL, active.default=TRUE,...){
  if(is.null(start) && !is.null(attr(x,"start"))) start <- attr(x,"start")
  if(is.null(end) && !is.null(attr(x,"end"))) end <- attr(x,"end")
  
  # deprecation warning because we want to use net.obs.period instead of invisible attrs
  if( !is.null(attr(x,"start")) |  !is.null(attr(x,"end"))){
    .Deprecated(msg='Indicating the network activity bounds for censoring with attrs for "start" and "end" attached to a network object has been deprecated. \nPlease use the network-level attribute "net.obs.period" instead.  \nSee ?net.obs.period for more information' )
  }

  # check if net.obs.period is present
  # use max and min for censoring bounds if not included.
  if(!is.null(x%n%'net.obs.period')){
    if(is.null(start)){
      start<-min(unlist((x%n%'net.obs.period')$observations))
    }
    if(is.null(end)){
      end<-max(unlist((x%n%'net.obs.period')$observations))
    }
  }
  
  # if start and end are still unset, set them to Inf so censoring and spell removal will work correctly
  if(is.null(start)){
    start<- -Inf
  }
  if(is.null(end)){
    end<- Inf
  }
  
  tm<-lapply(e,function(y){
    edge<-x$mel[[y]]
    if(is.null(edge)) NULL else{
      active<-edge$atl$active
      if (!is.null(active)) {
        ac<-matrix(rep(cbind(edge$outl,edge$inl,y),nrow(active)),ncol=3,byrow=TRUE)
        cbind(active,ac)
      } else {
        # activity attribute is missing, so use active default 
        if(active.default){
          # create -Inf Inf spell here
          matrix(c(-Inf,Inf,edge$outl,edge$inl,y),ncol=5,byrow=TRUE)
        }
      }
    }
  })
  out <- do.call(rbind,tm)
  
  # determine which spells are active and subset them
  if(!is.null(out) && nrow(out)>0 ){
    spls.active<-sapply(seq_len(nrow(out)),function(s){
      # odd syntax c(out[s,1],out[s,1]) is needed to force evaluation or something
      # otherwise the dataframe subset is *extremely* expensive
      spells.overlap(c(start,end),c(out[s,1],out[s,2]))
    })
    out<-out[spls.active,,drop=FALSE]
  }
  
  if (is.null(out)) {
    out <- data.frame(onset=numeric(), terminus=numeric(), tail=numeric(), head=numeric(),edge.id=numeric())
    #warning("Network does not have any edge activity")
  } else {
    colnames(out)<-c("onset","terminus","tail","head","edge.id")
  }
  out<-data.frame(out)
  
  # do censoring
  out$onset.censored <- out$onset < start | out$onset==-Inf
  out$terminus.censored <- out$terminus > end | out$terminus==Inf
  
  if(!is.null(start)) out$onset[out$onset.censored] <- start
  
  if(!is.null(end)) out$terminus[out$terminus.censored] <- end
  
  out$duration <- out$terminus-out$onset  
  
  # have to permute columns to put edge.id at end
  out<-out[,c(1,2,3,4,6,7,8,5)]

  # sort output by eid, onset,terminus
  out<-out[order(out[,8],out[,1],out[,2]),]
  
  return(out)
}

# return [onset, terminus, vertex]
get.vertex.spelllist = function (x, v=seq.int(network.size(x)), start=NULL, end=NULL,active.default=TRUE) {
  # todo: need to replace this with reading net.obs.period attribute
  if(is.null(start) && !is.null(attr(x,"start"))) start <- attr(x,"start")
  if(is.null(end) && !is.null(attr(x,"end"))) end <- attr(x,"end")
  
  # deprecation warning because we want to use net.obs.period instead of invisible attrs
  if( !is.null(attr(x,"start")) |  !is.null(attr(x,"end"))){
    .Deprecated(msg='Indicating the network activity bounds for censoring with attrs for "start" and "end" attached to a network object has been deprecated. \nPlease use the network-level attribute "net.obs.period" instead.  \nSee ?net.obs.period for more information' )
  }
  
  # check if net.obs.period is present
  # use max and min for censoring bounds if not included.
  if(!is.null(x%n%'net.obs.period')){
    if(is.null(start)){
      start<-min(unlist((x%n%'net.obs.period')$observations))
    }
    if(is.null(end)){
      end<-max(unlist((x%n%'net.obs.period')$observations))
    }
  }
  
  # grab the list of spells
  node.list <- get.vertex.activity(x,v=v,active.default=active.default)
  # deal with NA or nullcaused by vertices w/o spells defined
  if (active.default & length(node.list)>0){
    node.list[is.na(node.list)]<-list(insert.spell(NULL,onset=-Inf,terminus=Inf)) # always active
    node.list[is.null(node.list)]<-list(insert.spell(NULL,onset=-Inf,terminus=Inf)) # always active
  } 
  # don't include spell if null or na
  v<-v[!is.na(node.list)]
  node.list<-node.list[!is.na(node.list)]
  v<-v[!sapply(node.list,is.null)]
  node.list<-node.list[!sapply(node.list,is.null)]
  
  out <- lapply(seq_along(v), function(i){cbind(node.list[[i]][,1], node.list[[i]][,2],v[i])})
  out <- do.call(rbind, out)
  
  if (is.null(out)) {
    out <- data.frame(onset=numeric(), terminus=numeric(), vertex.id=numeric())
  } else {
    colnames(out) <- c("onset", "terminus", "vertex.id")
  }
  
  out <- data.frame(out)
  
  # figure out censoring
  out$onset.censored <- out$onset==-Inf
  out$terminus.censored <- out$terminus==Inf
  
  if(!is.null(start)) out$onset[out$onset.censored] <- start
  
  if(!is.null(end)) out$terminus[out$terminus.censored] <- end
  # figure out duration
  out$duration <- out$terminus-out$onset  
  
  # impose sort order
  #out <- out[order(out$onset, out$vertex),]
  out <- out[order(out$vertex.id, out$onset,out$terminus),]
  out
}

################
### end networkDynamic-> other formats
################

print.networkDynamic <- function(x, ...){
  cat("NetworkDynamic properties:\n")
  times <- get.change.times(x,ignore.inf=FALSE)
  if (length(times)==0){
    cat("  network contains no time information")
  } else {
    cat("  distinct change times:", length(times), "\n")
    maxrange<-range(times)
    cat("  maximal time range:", maxrange[1], "until ",maxrange[2],"\n")
  }
  # TEAs
  ntea <-list.network.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
  if (network.size(x)>0){
    vtea <-list.vertex.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
    etea <-list.edge.attributes.active(x,onset=-Inf,terminus=Inf,dynamic.only=TRUE)
  } else {
    vtea<-character(0)
    etea<-character(0)
  }
  if (length(ntea)+length(vtea)+length(etea)>0){
    cat("\n Dynamic (TEA) attributes:\n")
    if (length(ntea)>0){
      cat("  Network TEAs:")
      cat(paste("   ",ntea,"\n"))
    }
    if (length(vtea)>0){
      cat("  Vertex TEAs:")
      cat(paste("   ",vtea,"\n"))
    }
    if (length(etea)>0){
      cat("  Edge TEAs:")
      cat(paste("   ",etea,"\n"))
    }
  }
  
  # report on spell ranges for tea attributes?
  
  if ('net.obs.period'%in%list.network.attributes(x)){
    cat("\nIncludes optional net.obs.period attribute:\n")
    .print.net.obs.period(get.network.attribute(x,'net.obs.period',unlist=FALSE))
  }
  # print standard network stuff
  # todo: should we override print.network so we can supress printing of net.obs.period and info on TEAs?
  cat("\n")
  NextMethod("print")
}


##############
### as.networkDynamic
### converts various objects to networkDynamic
##############

is.networkDynamic <- function(x){
  "networkDynamic" %in% class(x)
}

as.networkDynamic <- function(object,...){
  UseMethod("as.networkDynamic")
}

# doesn't do anything. Returns object as is
as.networkDynamic.networkDynamic <- function(object,...){
  return(object)
}

# only sets the networkDynamic class
as.networkDynamic.network<- function(object,...){
  set.nD.class(object)
  return(object)
}

# remove networkDynamic class but leave object unchanged
as.network.networkDynamic<-function(x,...){
  if(is.networkDynamic(x)){
    class(x)<-class(x)[class(x)!='networkDynamic']
  }
  return(x)
}


# if x is an nD object, return x
# otherwise, modify it in its parent frame
set.nD.class <- function(x){
  if(!is.networkDynamic(x)) {
    xn <- substitute(x)
    class(x) <- c("networkDynamic", class(x))
    if(.validLHS(xn, parent.frame()))
      on.exit(eval.parent(call('<-',xn, x)))
    return(invisible(x))
  }
  return(x)
}

as.edgelist <- function(nw, attrname = NULL, as.sna.edgelist = FALSE,...){
  el <- as.matrix.network.edgelist(nw, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...)
  if(!is.directed(nw)) el[,1:2] <- cbind(pmin(el[,1],el[,2]),pmax(el[,1],el[,2]))
  el
}


# given a vector or list, return true if it is unique
is.unique.list <- function(x) {
  length(x) == length(unique(x))
}



## returns the min and max times of vertex and edge timings
mintime <- function(vertex.data, edge.data) {
  if ('time' %in% colnames(vertex.data)) t1 = 'time' else t1 = 'onset'
  if ('time' %in% colnames(edge.data)) t2 = 'time' else t2 = 'onset'
  if ((length(vertex.data[,t1]) + length(edge.data[,t2]))==0){
    return(-Inf)
  } else {
    return(min(vertex.data[,t1], edge.data[,t2]))
  }
}

maxtime <- function(vertex.data, edge.data) {
  if ('time' %in% colnames(vertex.data)) t1 = 'time' else t1 = 'terminus'
  if ('time' %in% colnames(edge.data)) t2 = 'time' else t2 = 'terminus'
  if ((length(vertex.data[,t1]) + length(edge.data[,t2]))==0){
    return(Inf)
  } else {
    return(max(vertex.data[,t1], edge.data[,t2]))
  }
}

## checks that net.obs.period object has appropriate structure

.check.net.obs.period <- function(x){
  if (!is.list(x)){
    stop("net.obs.period must be a list object")
  } 
  elements<-names(x)
  if(!'observations'%in%elements){
    stop('net.obs.period must contain a sub-list named "observations" containing at least one spell')
  }
  if (!all(is.numeric(unlist(x$observations)))){
    stop("all elements of the 'observations' component of net.obs.period must be numeric")
  }
  if (!all(sapply(x$observations,length)==2)){
    stop("all elements of the 'observations' component of net.obs.period must be a vector of length 2")
  }
  
  # must be possible to compute a max and min bound
  tryCatch({max(unlist(x$observations))},
           warning=function(w){
             stop("unable to compute a max and min value for the 'observations' component of net.obs.period", 
                  call. = FALSE)
           },
           error=function(e){
             stop("unable to compute a max and min value for the 'observations' component of net.obs.period", 
                  call. = FALSE)
           }
  )         
  
  if(!'mode'%in%elements){
    stop('net.obs.period must contain an element named "mode" indicating the time model')
  }
  if(!'time.increment'%in%elements){
    stop('net.obs.period must contain a element named "time.increment" which gives the natural time increment (but it can have the value NA if does not apply)')
  }
  if(!'time.unit'%in%elements){
    stop('net.obs.period must contain a element named "time.unit" which names the time unit for the network object')
  }
}

# internal function to pretty-print net.obs.period attribute
.print.net.obs.period<-function(nop){
  cat(" Network observation period info:\n")
  cat(paste("  Number of observation spells:",length(nop$observations),"\n"))
  maxrange<-range(nop$observations)
  cat(paste("  Maximal time range observed:",maxrange[1],"until",maxrange[2],"\n"))
  cat(paste("  Temporal mode:",nop$mode,"\n"))
  cat(paste("  Time unit:",nop$time.unit,"\n"))
  cat(paste("  Suggested time increment:",nop$time.increment,"\n"))
}

# function to adjust all of the activity spells in a network by a specified amount
adjust.activity <-function(nd,offset=0,factor=1){
  # check args
  if (!is.networkDynamic(nd)){
    stop("adjust.activity requires a networkDynamic object as its first argument")
  }
  if (!is.numeric(offset)){
    stop("the offset argument must be a positive or negative numeric value giving the amount of the time adjustment")
  }
  
  if (!is.numeric(factor)){
    stop("the factor argument must be a positive or negative numeric value giving the amount of the time should be multiplied by")
  }
  xn <- substitute(nd)   # needed for proper assignment in calling environment
  
  # change all the vertex spells
  nd$val<-lapply(nd$val, function(v){
    if (!is.null(v$active)){
      v$active<-(v$active+offset)*factor
    }
    v
  })
  
  # change all the vertex teas
  vattributes <- list.vertex.attributes(nd)
  vattributes.active <- grep(".active", vattributes, value = TRUE)
  for (attrname in vattributes.active){
    nd$val<-lapply(nd$val, function(v){
      if (!is.null(v[[attrname]])){
        v[[attrname]][[2]]<-(v[[attrname]][[2]]+offset)*factor
      }
      v
    })
  }
  
  # change all the edge spells
  nd$mel<-lapply(nd$mel, function(e){
    if (!is.null(e$atl$active)){
      e$atl$active<-(e$atl$active+offset)*factor
    }
    e
  })
  
  # change all the edge teas
  nd$mel<-lapply(nd$mel, function(e){
    # get the list of attribute names for the edge
    eattributes.active <- grep(".active", names(e$atl), value = TRUE)
    for(attrname in eattributes.active)
    if (!is.null(e$atl[[attrname]])){
      e$atl[[attrname]][[2]]<-(e$atl[[attrname]][[2]]+offset)*factor
    }
    e
  })
  
  # change the network teas
  nattributes <- list.network.attributes(nd)
  nattributes.active <- grep(".active", nattributes, value = TRUE)
  for (attrname in nattributes.active){
      nd$gal[[attrname]][[2]]<-(nd$gal[[attrname]][[2]]+offset)*factor
  }
  
  # change the net.obs.period
  obs<-nd%n%'net.obs.period'
  if(!is.null(obs)){
    # transform the observation spells
    obs$observations<-lapply(obs$observations,function(observation){
      (observation+offset)*factor
    })
    # also transform the time.increment
    obs$time.increment<-obs$time.increment*factor
    nd%n%'net.obs.period'<-obs
  }
  if(.validLHS(xn, parent.frame()))
    on.exit(eval.parent(call('<-',xn, nd)))
  invisible(nd)
}


# function to create edge spell matrix for a single dyad (tail,head)
# from matrix of changes [time, tail, head, direction]
spellsFromChanges<-function(changes,tail,head,strict=TRUE){
  dyadRows<-changes[,2]==tail & changes[,3]==head
  changes<-changes[dyadRows,,drop=FALSE]
  # remove any row names that may exist
  dimnames(changes)<-NULL
  # if there are only three columns, we are actually dealing with toggles
  # assume that the first toggle is an activation
  if (ncol(changes)==3){
    # append 1 and 0 to indicate alternating activation and deactivation
    changes<-cbind(changes,1:nrow(changes)%%2)
  } else {
    # check we don't have unbalenced activations or deactivations
    badRows<-findRep(changes[,4])
    if(length(badRows)>0){
      if (strict){
        stop('encountered unbalanced changes for dyad ',tail,' ',head, ' at times ',paste(changes[badRows,1]))
      } else {
        # painfully find and remove offending spell row(s)
        changes<-changes[-badRows,,drop=FALSE]
      }
    }
  }
  # if there is an odd number of rows, spells will be unbalenced
  # so need to pad beginning or end with inf
  if (nrow(changes)%%2>0){
    if (changes[1,4]==1){ # if first spell was an activation
      # last spell should be open interval
      changes<-rbind(changes,c(Inf,tail,head,0))
    } else {  # first spell was deactivation, so first spell should be left-open interval activation
      changes<-rbind(c(-Inf,tail,head,1),changes)
    }
  }
  # bind onsets and termini into a spell matrix
  splmat<-cbind(changes[changes[,4]==1,1,drop=FALSE],changes[changes[,4]==0,1,drop=FALSE])
  return(splmat)
}


# find the index of values which are repeats of the previous values
findRep<-function(x){
  if (length(x)>1){
   # find cases where vector value matches its value offset by 1
   return(which(x[1:(length(x)-1)]==x[2:length(x)])+1)
  } else {
   return(numeric(0))
  }        
}

