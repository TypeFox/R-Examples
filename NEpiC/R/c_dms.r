c_dms <-
  function(network, gene2weight, d=1,r=0.1)
  {
    ###########################################################################################################
    ######################################## Starting dense module searching #########################################
    ####################################################################################################################
    #==================================================================================================================#
       
    cat("genes used: ", length(gene2weight[,1]), "\n", sep="")
    
    rawG <- graph.data.frame(network,directed=F)
    g.weight <- as.numeric(gene2weight[,2])
    intG <- integGM(rawG,as.character(gene2weight[,1]),g.weight)
    GWPI <- simplify(intG)
    
    cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
    dm.result <- globalQueryJ(GWPI,search_r=d,r)
    
    #================================================ extract genesets ================================================#
    cat("extracting modules...\n", sep="")
    genesets <- list()
    for(k in 1:length(dm.result))
    {
      node = names(dm.result[k])
      g = dm.result[[k]]
      genesets[[node]] <- V(g)$name
    }
    seed.genes <- names(genesets)
    
    #================================== clean genesets by removing identical records ==================================#
    cat("removing identical modules...\n", sep="")
    identical.idx <- list()
    for(k in 1:length(genesets))
    {
      tmp.idx <- c(k);
      for(kt in 1:length(genesets))
      {
        if(kt==k)next()
        genesk <- genesets[[k]]
        genest <- genesets[[kt]]
        if(length(genesk)!=length(genest))next()
        overlap = intersect(genesk, genest)
        if(length(overlap)==length(genest))
        {
          tmp.idx <- c(tmp.idx, kt)
        }
      }
      if(length(tmp.idx)>1)
      {
        tmp.idx <- sort(tmp.idx)
        identical.idx[[seed.genes[k]]] <- tmp.idx
        #cat(k, ".", sep="")
      }
    }
    
    toremove.idx <- c()
    for(k in 1:length(identical.idx))
    {
      tmp.idx <- identical.idx[[k]]
      toremove.idx <- c(toremove.idx, tmp.idx[-1])
    }
    toremove.idx <- unique(toremove.idx)
    genesets.clear <- genesets[-toremove.idx]
    
    #================================================= random network =================================================#
    cat("permutation on random network...\n", sep="")
    genesets.length <- c()
    for(k in 1:length(genesets.clear))
    {
      genes <- genesets.clear[[k]]
      genesets.length <- c(genesets.length, length(genes))
    }
    genesets.length <- unique(genesets.length)
    
    genes.idx <- seq(1, length(V(GWPI)$name))
    graph.g.weight = data.frame(GWPIene=V(GWPI)$name, gain.weight=V(GWPI)$weight)
    genesets.length.null.dis <- list()
    length.max = max(genesets.length)+5
    for(k in 5:length.max)
    {
      l.zperm <- c()
      for(j in 1:100000)
      {
        idx.pseudo=sample(genes.idx, size=k)
        l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo, 2])/sqrt(length(idx.pseudo)));
      }
      genesets.length.null.dis[[as.character(k)]] = l.zperm
      cat(k, ".", sep="");
    }
    
    genesets.length.null.stat <- list()
    for(k in 5:length.max)
    {
      l.zperm <- genesets.length.null.dis[[as.character(k)]]
      k.mean <- mean(l.zperm)
      k.sd <- sd(l.zperm)
      genesets.length.null.stat[[as.character(k)]] = c(k.mean, k.sd)
    }
    
    ################################################### Normalization ##################################################
    #==================================================================================================================#
    zim <- data.frame(gene=names(genesets.clear), Zm=-9, Zn=-9, zcount=-9)
    for(k in 1:length(genesets.clear))
    {
      genes <- genesets.clear[[k]]
      match(genes, graph.g.weight[,1]) -> idx
      idx <- idx[!is.na(idx)]
      zim[k,2] <- sum(graph.g.weight[idx, 2])/sqrt(length(idx)) 
      
      tmp <- genesets.length.null.stat[[as.character(length(idx))]]
      zim[k, 3]=(zim[k,2]-tmp[1])/tmp[2]
      zim[k, 4]=sum(genesets.length.null.dis[[as.character(length(idx))]]>=zim[k,2]) 
    }
    zom = zim[order(zim[,3], decreasing=T), ]
    #==================================================================================================================#
    
    #================================================== save results ==================================================#
    res.list <- list()
    res.list[["GWPI"]]                        = GWPI
    res.list[["graph.g.weight"]] 		          = graph.g.weight
    res.list[["genesets.clear"]] 		          = genesets.clear
    res.list[["genesets.length.null.dis"]] 	  = genesets.length.null.dis
    res.list[["genesets.length.null.stat"]] 	= genesets.length.null.stat
    res.list[["zi.matrix"]]                   = zim
    res.list[["zi.ordered"]]                  = zom
    save(res.list, file="RESULT.list.RData")
    cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
    return(res.list);
    ########################################### End of dense module searching ###########################################
  }

