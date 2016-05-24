
plot.gexf <- function(x, EdgeType = c("curve", "line"), output.dir = NULL, ...){  
#   if(!is.null(gexf.object$positions)){
#     library(sna)
#     nNodes <- nrow(gexf.object$nodes)
#     links <- matrix(rep(0, nNodes*nNodes), ncol = nNodes)
#     relations <- gexf.object$edges[,c(3,4)]
#     
#     for(edge in 1:ncol(relations)){
#       links[(relations[edge,]$target), (relations[edge,]$source)] <- 1
#     }
#     
#     positions <- gplot.layout.kamadakawai(links, layout.par=list())
#     positions <- cbind(positions, 0) # needs a z axis
#   }
#   
#   nodecolors <- data.frame(r = rep(255, nNodes),
#                            g = rep(255, nNodes),
#                            b = rep(255, nNodes),
#                            a = rep(255, nNodes))
#   
#   gexf.object <-  write.gexf(nodes=gexf.object$nodes,
#                              edges=gexf.object$edges[,c("source","target")],
#                              nodesVizAtt=list(
#                                color=nodecolors,
#                                position=positions
#                              ))
#   
  if(length(unlist(x$atts.definitions))){
    html <- readLines(system.file("sigmajs/index_att.html", package="rgexf"), warn=FALSE)
  } else {
    html <- readLines(system.file("sigmajs/index.html", package="rgexf"), warn=FALSE)    
  }

  html <- gsub("EdgeTypePar", EdgeType[1], html)
  
  s <- Rhttpd$new()
  s$start(listen='127.0.0.1')
  
  # html
  my.app <- function(env) {  
    res <- Response$new()
    res$write(paste(html, collapse="\n"))
    res$finish()
  }
  s$add(app=my.app, name='plot')
  
  # graph
  graph <- function(env){
    res <- Response$new()
    res$write(x$graph)
    res$finish()
  }
  s$add(app=graph, name='data')
  
  # jquery
  jquery <- function(env){
    res <- Response$new()
    res$write(paste(readLines(system.file("sigmajs/jquery.min.js", package="rgexf"), warn=FALSE), collapse="\n "))
    res$finish()
  }
  s$add(app=jquery, name='jquery.js')
  
  
  # sigmajs
  sigmajs <- function(env){
    res <- Response$new()
    res$write(paste(readLines(system.file("sigmajs/sigma.min.js", package="rgexf"), warn=FALSE) , collapse="\n "))
    res$finish()
  }
  s$add(app=sigmajs, name='sigmajs')
  
  
  # parseGexf
  parseGexf <- function(env){
    res <- Response$new()
    res$write(paste(readLines(system.file("sigmajs/sigma.parseGexf.js", package="rgexf"), warn=FALSE), collapse="\n "))
    res$finish()
  }
  s$add(app=parseGexf, name='sigmaparseGexfjs')
  s$browse('plot')
  
  if(length(output.dir)){
    
    wd <- getwd()
    setwd(output.dir)
    
    html2 <- gsub("custom", "js", html)
    html2 <- gsub("js\">", ".js\">", html2)
    html2 <- gsub("/js/data", "graph.gexf", html2)
    
    writeLines(html2, "index.html")
    writeLines(x$graph, "graph.gexf")
    
    if(!file.exists("js")) dir.create("js")
    
    writeLines(paste(readLines(system.file("sigmajs/jquery.min.js", package="rgexf"), warn=FALSE), collapse="\n "),
               file.path("js", "jquery.js"))
    writeLines(paste(readLines(system.file("sigmajs/sigma.min.js", package="rgexf"), warn=FALSE), collapse="\n "),
               file.path("js", "sigma.js"))
    writeLines(paste(readLines(system.file("sigmajs/sigma.parseGexf.js", package="rgexf"), warn=FALSE), collapse="\n "),
               file.path("js", "sigmaparseGexf.js"))
    
    setwd(wd)
  }
  
}
