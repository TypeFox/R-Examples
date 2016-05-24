plot.rcox <- function(x,y,...){


##   if (!("package:Rgraphviz" %in% search())){
##     if("Rgraphviz" %in% rownames(installed.packages())){
##       require(Rgraphviz)
##     } else {
##       cat("The Rgraphviz package (from Bioconductor) must be installed to display the models\n")
##       return()
##     }
##   }

  
##   if (!("package:Rgraphviz" %in% search())){
##     cat("The Rgraphviz package (from Bioconductor) must be installed...\n")
##     return()
##   }
  
##     if("Rgraphviz" %in% installed.packages()){
##       library(Rgraphviz)
##     } else {
##       cat("The Rgraphviz package (from Bioconductor) must be loaded to display the models\n")
##       return()
##     }
##   }

  if (!require("Rgraphviz")){
    cat("The Rgraphviz package (from Bioconductor) must be installed to display the models\n")
    return()
  }

##   if(!("Rgraphviz" %in% rownames(installed.packages()))){
##     cat("The Rgraphviz package (from Bioconductor) must be installed to display the models\n")
##     return()
##   }

    

  
  m2 <- x


  eccList <- getSlot(m2,'ecc')
  vccList <- getSlot(m2,'vcc')

  coef <- coef(x)

  if (is.null(coef))
    coef <- 1:c(length(eccList)+length(vccList))
  
  #l <- sapply(vccList, length)
  
  o <- order(coef[1:length(vccList)])
  vccColors <- heat.colors(length(vccList))
  vccColors <- vccColors[o]
  
  V <- unlist(vccList)
  V <- V[order(V)]
  vertexColors <- NULL
  for (i in 1:length(vccList)){
    tmp <- vccList[[i]]
    #print(tmp)
    if (length(tmp)==1){
       vcolor <- "white"    
    } else {
       vcolor <- vccColors[i] 
    }
    d <- c(vstr = rep(vcolor, length(tmp)))
    names(d) <- tmp
    vertexColors <- c(vertexColors, d)    
  }

  nAttrs <- list()
  nAttrs$fillcolor <- vertexColors

  edL <- vector("list", length=length(V))
  names(edL) <- V
  nv <- 1:length(V)
  names(nv) <- V

  ed <- unlist(eccList,recursive=FALSE)

  for (i in 1:length(V)){
    idx <- sapply(ed, function(x) is.element(V[i],x))
    e <- setdiff (unlist(ed[idx]),V[i])
    edL[[V[i]]] <- list(edges=nv[e])
  }
  edL<- edL[sapply(edL,length)>0]

  ##G <- new("graphNEL", nodes=V, edgeL=edL)
  G <- new("graphNEL", nodes=V, edgeL=edL)

  ##G <- new("graphNEL", nodes=V)
  
  eccColors<-topo.colors(length(eccList))

  edgeColors <- NULL
  if (length(eccList)>0){
    o <- order(coef[-(1:length(vccList))])
    eccColors <- eccColors[o]

    for (i in 1:length(eccList)){
      tmp <- eccList[[i]]; ltmp <- length(tmp)
      for (j in 1:ltmp){
        ee <- tmp[[j]]
        ee <- ee[order(ee)]
        #G <- addEdge(ee[1], ee[2], G, weight=1)
        estr <- paste(ee[1],"~",ee[2],sep='')
        if (ltmp > 1){
          ecolor <- eccColors[i]
          d <- c(estr = ecolor)
          names(d) <- estr
          edgeColors <- c(edgeColors, d)                  
      }
      }
    }
  }

  #print(edgeColors)
  if (!is.null(edgeColors))
    eAttrs <- list(color=edgeColors)
  else
    eAttrs <- list()
  
  plot(G, "neato", nodeAttrs = nAttrs, edgeAttrs = eAttrs)
  #return(G)
}
