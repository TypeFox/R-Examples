################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for combining other endnodes conditional on factor attributes: 
# class "utility.endnode.cond"
# ==============================================================================


# constructor:
# ------------

utility.endnode.cond.create <- function(name.node,          # character(1)
                                        attrib.levels,      # data.frame
                                        nodes,              # list of nodes
                                        utility      = TRUE,
                                        required     = FALSE,
                                        col          = "black",
                                        shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
  if ( !is.data.frame(attrib.levels) )
  {
    cat("*** Warning: attrib.levels must be a data frame","\n")
    check.ok <- F
  }
  if ( length(names(attrib.levels)) != length(unique(names(attrib.levels))) )
  {
    cat("*** Warning: cColumn names of attrib.levels must be different","\n")
    check.ok <- F
  }
  if ( nrow(attrib.levels) != length(nodes) )
  {
    cat("*** Warning: Number of rows of attrib.levels not equal to",
        "number of nodes provided:",nrow(attrib.levels),length(nodes),"\n")
    check.ok <- F
  }
  if ( length(nodes) < 1 )
  {
    cat("*** Warning: No nodes provided","\n")
    check.ok <- F
  }
  for ( i in 1:length(nodes) )
  {
    if ( nodes[[i]]$utility != utility )
    {
      funtype   <- "utility"; if ( !utility )            funtype   <- "value"
      funtype.i <- "utility"; if ( !nodes[[i]]$utility ) funtype.i <- "value"
      cat("***Warning: incompatible function types: new node is of type",
          funtype,"node",nodes[[i]]$name," is of type",funtype.i,"\n")
      check.ok <- F
    }
  }
  if ( ! check.ok )
  {
    cat("*** Warning: Node \"",name.node,"\" could not be constructed","\n",
        sep="")
    return(NA)
  }
  
  # construct class:
  
  node <- list()
  node$name          <- name.node
  node$description   <- "utility/value conditional combination end node"
  node$type          <- "endnode"
  node$attrib.levels <- attrib.levels
  for ( i in 1:ncol(attrib.levels) ) 
  {
    node$attrib.levels[,i] <- as.character(node$attrib.levels[,i])
  }
  node$attrib        <- names(attrib.levels)
  for ( i in 1:length(nodes) )
  {
    node$attrib     <- c(node$attrib,nodes[[i]]$attrib)
  }
  node$attrib        <- unique(node$attrib)
  node$nodes         <- nodes
  node$required      <- required
  node$utility       <- utility
  node$col           <- col
  node$shift.levels  <- shift.levels
  class(node)        <- "utility.endnode.cond" 
  
  # print return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.cond <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update conditional nodes:
  
  n <- node
  for ( i in 1:length(n$nodes) )
  {
    n$nodes[[i]] <- updatepar(n$nodes[[i]],par)
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.cond <- function(x,
                                          attrib,   # data.frame
                                          par=NA,
                                          ...)
{
  node <- x
  
  # check availability of attributes:
  
  if ( ! is.data.frame(attrib) )
  {
    warning("Node \"",node$name,"\": attrib must be a data frame",sep="")
    return(NA)
  }
  ind <- match(node$attrib,names(attrib))
  if(sum(is.na(ind))>0)
  {
    ind.na <- is.na(ind)
    warning("Node \"",node$name,"\": attribute(s) \"",
            paste(node$attrib[ind.na],collapse=","),"\" not found",sep="")
    return(rep(NA,nrow(attrib)))
  }
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # select rows compatible with conditioning attributes:
  
  u    <- rep(NA,nrow(attrib))
  calc <- rep(FALSE,nrow(attrib))
  for ( i in 1:ncol(n$attrib.levels) )  # evaluate NAs
  {
    calc <- calc | is.na(attrib[,names(n$attrib.levels)[i]])
  }
  while( TRUE )
  {
    # identify first row that has not yet been evaluated:
    
    startind <- match(FALSE,calc)
    
    # break if all were evaluated:
    
    if ( is.na(startind) ) break
    
    # find rows with the same attribute combinations:
    
    ind.attrib <- as.character(attrib[startind,names(n$attrib.levels)[1]]) == 
      as.character(attrib[,names(n$attrib.levels)[1]])
    if ( ncol(n$attrib.levels) > 1 )
    {
      for ( i in 2:ncol(n$attrib.levels) )
      {
        ind.attrib <- 
          ind.attrib & 
          ( as.character(attrib[startind,names(n$attrib.levels)[i]]) == 
              as.character(attrib[,names(n$attrib.levels)[i]]) )
      }
    }
    ind.attrib <- which(ind.attrib)
    
    # find corresponding node:
    
    ind.node <- as.character(attrib[startind,names(n$attrib.levels)[1]]) ==
      as.character(n$attrib.levels[,names(n$attrib.levels)[1]])
    if ( ncol(n$attrib.levels) > 1 )
    {
      for ( i in 2:ncol(n$attrib.levels) )
      {
        ind.node <- 
          ind.node & 
          ( as.character(n$attrib.levels[,names(n$attrib.levels)[i]]) == 
              as.character(attrib[startind,names(n$attrib.levels)[i]]) )
      }
    }
    ind.node <- which(ind.node)
    
    # evaluate node for all attribute rows with same conditional values:
    
    if ( length(ind.node) > 0 )
    {      
      u[ind.attrib] <- evaluate(n$nodes[[ind.node[1]]],attrib[ind.attrib,])
      if ( length(ind.node) > 1 )
      {
        cat("*** Warning: multiple combinations of the same",
            "attribute levels in node",n$name,"\n")
      }
    }            
    calc[ind.attrib] <- T
  }      
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.endnode.cond <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.cond <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("attribute/node combinations:","\n")
  nodes.names <- character(0)
  for ( i in 1:length(node$nodes) ) nodes.names[i] <- node$nodes[[i]]$name
  print(cbind(node$attrib.levels,node=nodes.names))
  for ( i in 1:length(node$nodes) ) 
  {
    cat("**","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.endnode.cond <-
  function(x,
           par       = NA,
           col       = utility.calc.colors(),
           gridlines = c(0.2,0.4,0.6,0.8),
           main      = "",
           cex.main  = 1,
           nodes     = x$name,
           ...)
  {
    node <- x
    if ( is.na(nodes[1]) | ! is.na(match(node$name,nodes)) )
    {
      nrow <- floor(sqrt(length(node$nodes)))
      ncol <- floor(length(node$nodes)/nrow+0.999)
      par.def <- par(no.readonly=T)
      par(mfrow=c(nrow,ncol),mar=c(4.3,3.8,2.8,0.8),oma=c(0,0,2,0)) 
      for ( i in 1:length(node$nodes) )             # c(bottom, left, top, right)
      {
        title <- main
        for ( j in 1:ncol(node$attrib.levels) )
        {
          title <- paste(title," ",colnames(node$attrib.levels)[j],"=",
                         as.character(node$attrib.levels[i,j]),sep="")
        }
        plot(node$nodes[[i]],par=par,col=col,gridlines=gridlines,main=title,cex.main=cex.main,...)
      }
      mtext(node$name,outer=TRUE,cex=cex.main)
      par(par.def)
    }
    if ( length(node$nodes) > 0 )
    {
      for ( i in 1:length(node$nodes) )
      {
        if ( is.na(nodes[1]) | !is.na(match(node$nodes[[i]]$name,nodes)) )
        {
          plot(node$nodes[[i]],
               par=par,
               col=col,
               gridlines=gridlines,
               cex.main=cex.main,
               ...)
        }
      }
    }
  }

