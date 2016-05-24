### #####################################################
###
### Plot iModels
### "D"iscrete variables are "D"ots (=gray)
### "C"ontinuous variabesl are "C"ircles (=transparent)
###
### #####################################################

iplot.iModel <- function(x,...){
  ig <- ugList(x$glist, result="igraph")
  V(ig)$label <- V(ig)$name
  V(ig)$size  <- 50
  ig$cex      <-  4
  V(ig)$label.cex <- 1.2

  switch(class(x)[1],
         "dModel"={
           V(ig)$color <- "grey"
         },
         "cModel"={
           V(ig)$color <- "white"
         },
         "mModel"={
           V(ig)$color <- "white"
           disc.idx <- match(x$datainfo$disc.names, V(ig)$name) #-1
           V(ig)[disc.idx]$color <- "grey"
         })

  ig$layout <- layout.lgl
  plot(ig)
  return(invisible(x))
}

plot.iModel <- function(x,...){
  uG <- ugList(x$glist)
  switch(class(x)[1],
         "dModel"={
           fillv <- rep("lightgray", length(x$varNames))
           names(fillv) <- x$varNames
         },
         "cModel"={
           fillv <- rep("transparent", length(x$varNames))
           names(fillv) <- x$varNames
         },
         "mModel"={
           dv <- x$datainfo$disc.names
           cv <- x$datainfo$cont.names
           fillv  <- c(rep("lightgray", length(dv)), rep("transparent", length(cv)))
           names(fillv) <- c(dv,cv)
         })

  plot(uG, nodeAttrs=list(fillcolor=fillv))
}
