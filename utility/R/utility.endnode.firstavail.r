################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for getting the results of the first node of a list that can 
# successfully be evaluated
# class "utility.endnode.firstavail"
# ==============================================================================


# constructor:
# ------------

utility.endnode.firstavail.create <- function(name.node,          # character(1)
                                              nodes,              # list of nodes
                                              utility      = TRUE,
                                              required     = FALSE,
                                              col          = "black",
                                              shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
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
  node$description   <- "utility/value endnode to evaluate first available subnode"
  node$type          <- "endnode"
  node$attrib        <- character(0)
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
  class(node)        <- "utility.endnode.firstavail" 
  
  # print return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.firstavail <- function(x,par=NA,...)
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

evaluate.utility.endnode.firstavail <- function(x,
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
  
  # update parameters:
  
  n <- updatepar(node,par)

  # evaluate nodes:
  
  u <- rep(NA,nrow(attrib))
  for ( i in 1:nrow(attrib) )
  {
    for ( j in 1:length(node$nodes) )
    {
      u[i] <- evaluate(node$nodes[[j]],attrib[i,])
      if ( !is.na(u[i]) ) break
    }
  }
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.endnode.firstavail <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.firstavail <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("nodes:","\n")
  for ( i in 1:length(node$nodes) ) cat("  ",node$nodes[[i]]$name,"\n")
  for ( i in 1:length(node$nodes) ) 
  {
    cat("**","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.endnode.firstavail <-
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
        title <- paste(main,i,":",node$nodes[[i]]$name)
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

