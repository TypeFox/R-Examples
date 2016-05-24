################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for discrete factor attributes: 
# class "utility.endnode.discrete"
# ==============================================================================


# constructor:
# ------------

utility.endnode.discrete.create <- function(name.node,          # character(1)
                                            attrib.levels,      # data.frame
                                            u,                  # numeric(n)
                                            names.u      = rep(NA,length(u)),
                                            utility      = TRUE,
                                            required     = FALSE,
                                            col          = "black",
                                            shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
  if ( !is.data.frame(attrib.levels) )
  {
    cat("*** Warning: Attrib.levels must be a data frame","\n")
    check.ok <- F
  }
  if ( length(names(attrib.levels)) != length(unique(names(attrib.levels))) )
  {
    cat("*** Warning: Column names of attrib.levels must be different","\n")
    check.ok <- F
  }
  if ( nrow(attrib.levels) != length(u) )
  {
    cat("*** Warning: Number of rows of attrib.levels not equal to",
        "number of elements of u:",nrow(attrib.levels),length(u),"\n")
    check.ok <- F
  }
  if ( length(names.u) != length(u) )
  {
    cat("*** Warning: Number of elements of names.u not equal",
        "number of elements of u:",length(names.u),length(u),"\n")
    check.ok <- F
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
  node$description   <- "utility/value discrete attribute end node"
  node$type          <- "endnode"
  node$attrib.levels <- attrib.levels
  for ( i in 1:ncol(attrib.levels) ) 
  {
    node$attrib.levels[,i] <- as.character(node$attrib.levels[,i])
  }
  node$attrib       <- names(attrib.levels)
  node$u            <- u
  node$names.u      <- names.u
  node$required     <- required
  node$utility      <- utility
  node$col          <- col
  node$shift.levels <- shift.levels
  class(node)       <- "utility.endnode.discrete" 
  
  # print and return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.discrete <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update adequate values in interpolation list:
  
  n <- node
  for ( i in 1:length(n$u) )
  {
    if ( ! is.na(n$names.u[i]) )
    {
      ind <- which(n$names.u[i] == names(par) )
      if ( length(ind) > 1 )
      {
        warning("Node \"",node$name,"\": multiple occurrences of parameter",
                names(par)[ind[1]])
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$u[i] <- par[ind]
      }
    } 
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.discrete <- function(x,
                                              attrib,   # data.frame
                                              par = NA,
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
  
  # check levels of attributes:
  
  for ( i in 1:ncol(node$attrib.levels) )
  {
    n <- names(node$attrib.levels)[i]                       # attribute name
    l <- unique(node$attrib.levels[,i]); l <- l[!is.na(l)]  # defined levels
    a <- unique(attrib[,n]); a <- a[!is.na(a)]              # requested levels
    ind.na <- is.na(match(a,l))
    if ( sum(ind.na) > 0 )
    {
      warning("Node \"",node$name,"\": unknown attribute level(s): \"",paste(a[ind.na],collapse=","),
              "\" of attribute \"",n,"\"",sep="")
    }
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
    
    # find corresponding value:
    
    ind.u <- as.character(attrib[startind,names(n$attrib.levels)[1]]) ==
      as.character(n$attrib.levels[,names(n$attrib.levels)[1]])
    if ( ncol(n$attrib.levels) > 1 )
    {
      for ( i in 2:ncol(n$attrib.levels) )
      {
        ind.u <- 
          ind.u & 
          ( as.character(n$attrib.levels[,names(n$attrib.levels)[i]]) == 
              as.character(attrib[startind,names(n$attrib.levels)[i]]) )
      }
    }
    ind.u <- which(ind.u)
    
    # evaluate node for all attribute rows with same conditional values:
    
    if ( length(ind.u) == 1 )
    {
      u[ind.attrib] <- n$u[ind.u]
    }
    else
    {
      if ( length(ind.u) > 1 )
      {
        warning("Node \"",node$name,"\": multiple combinations of the same",
                "attribute levels in node \"",n$name,"\"",sep="")
      }
    }            
    calc[ind.attrib] <- T
  }      
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.endnode.discrete <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.discrete <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  cat("attribute(s):   ",paste(node$attrib,collapse=","),"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("attribute/value combinations:","\n")
  names.u <- ifelse(is.na(node$names.u),"",node$names.u)
  print(cbind(node$attrib.levels,u=node$u,names.u=names.u))
}


# plot:
# -----

plot.utility.endnode.discrete <- 
  function(x,
           par       = NA,
           col       = utility.calc.colors(),
           gridlines = c(0.2,0.4,0.6,0.8),
           main      = "",
           cex.main  = 1,
           ...)
  {
    # plot frame:
    
    node <- x
    length = 101
    n <- updatepar(node,par)
    title <- main; if ( nchar(title) == 0 ) title <- n$name
    funtype <- "utility"; if ( !n$utility ) funtype <- "value"
    plot(numeric(0),numeric(0),type="l",
         xlim=c(0,1),ylim=c(0,1),
         xlab=paste(n$attrib,collapse=","),ylab=funtype,main=title,
         xaxs="i",yaxs="i",yaxt="n",xaxt="n",cex.main=cex.main,...)
    
    # colored bar along y axis:
    
    if ( length(col)>1 & !node$utility )
    {
      num.grid = 100
      endpoints <- seq(0,1,length.out=num.grid+1)+1/(2*num.grid)
      midpoints <- 0.5*(endpoints[-1]+endpoints[-length(endpoints)])
      cols <- utility.get.colors(midpoints,col)
      for ( i in 1:(num.grid-1) )
      {
        lines(-0.01*c(1,1),endpoints[c(i,i+1)],col=cols[i],lwd=3,lend=2,xpd=TRUE)
      }
    }
    
    # axes (should overly colored bar):
    
    labels=character(length(n$u))
    for ( i in 1:length(n$u) )
    {
      labels[i] <- paste(as.character(n$attrib.levels[i,]),collapse=",")
    }
    axis(side=1,at=((1:length(n$u))-0.5)/length(n$u),labels=labels)
    axis(side=2)
    
    # plot gridlines:
    
    if ( !node$utility )
    {
      if ( !is.na(gridlines[1]) )
      {
        for ( level in gridlines ) abline(h=level,lty="dashed")
      }
    }
    
    # plot points:
    
    color <- "black"; if(length(col)==1) color <- col
    points(((1:length(n$u))-0.5)/length(n$u),n$u,pch=19,col=color,xpd=TRUE)
  }

