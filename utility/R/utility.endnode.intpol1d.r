################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for single-attribute interpolation: 
# class "utility.endnode.intpol1d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.intpol1d.create <- function(name.node,    # character(1)
                                            name.attrib,  # character(1)
                                            range,        # numeric(2)
                                            x,            # numeric(n)
                                            u,            # numeric(n)
                                            names.x     = rep(NA,length(x)),
                                            names.u     = rep(NA,length(u)),
                                            utility     = TRUE,
                                            required    = FALSE,
                                            col         = "black",
                                            shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
  if ( length(x) != length(u) )
  {
    cat("*** Warning: x and u of different length:",
        length(x),length(u))
    check.ok <- F
  }
  if ( length(names.x) != length(names.u) )
  {
    cat("*** Warning: names.x and names.u of different length:",
        length(names.x),length(names.u),"\n")
    check.ok <- F
  }
  if ( length(x) != length(names.x) )
  {
    cat("*** Warning: x and names.x of different length:",
        length(x),length(names.x),"\n")
    check.ok <- F
  }
  if ( range[1] >= range[2] )
  {
    cat("*** Warning: Minimum of range not smaller than maximum:",
        range[1],range[2],"\n")
    check.ok <- F
  }
  if ( sum(x[-1]-x[-length(x)] > 0) != length(x)-1 &
         sum(x[-1]-x[-length(x)] < 0) != length(x)-1 )
  {
    cat("*** Warning: x values in interpolation node must either be","\n",
        "strictly increasing or strictly decreasing","\n")
    check.ok <- F
  } 
  if ( ! check.ok )
  {
    cat("*** Warning: node \"",name.node,"\" could not be constructed","\n",
        sep="")
    return(NA)
  }
  
  # construct class:
  
  node <- list()
  node$name         <- name.node
  node$description  <- "utility/value 1d interpolation end node" 
  node$type         <- "endnode"
  node$attrib       <- name.attrib
  node$range        <- range
  node$x            <- x
  node$u            <- u
  node$names.x      <- names.x
  node$names.u      <- names.u
  node$required     <- required
  node$utility      <- utility
  node$col          <- col
  node$shift.levels <- shift.levels
  class(node)       <- "utility.endnode.intpol1d" 
  
  # print and return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.intpol1d <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update adequate values in interpolation list:
  
  n <- node
  for ( i in 1:length(n$x) )
  {
    if ( ! is.na(n$names.x[i]) )
    {
      ind <- which(n$names.x[i] == names(par) )
      if ( length(ind) > 1 )
      {
        warning("Node \"",node$name,"\": multiple occurrences of parameter \"",
                names(par)[ind[1]],"\"",sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$x[i] <- par[ind]
      }
    } 
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

evaluate.utility.endnode.intpol1d <- function(x,
                                              attrib,   # data.frame, numeric
                                              par = NA,
                                              ...)
{
  node <- x
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # extract attributes:
  
  if ( is.data.frame(attrib) | is.matrix(attrib) )
  {
    if ( length(which(colnames(attrib)==n$attrib)) != 1 )
    {
      warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
      return(rep(NA,nrow(attrib)))
    }
    a <- attrib[,n$attrib]
  }
  else
  {
    if ( ! is.vector(attrib) )
    {
      warning("Node \"",node$name,"\": unknown format of attribute \"",n$attrib,"\"",sep="")
      return(NA)
    }
    if ( length(names(attrib)) == 0 )
    {
      a <- attrib
    }
    else
    {
      ind <- which(names(attrib)==n$attrib)
      if ( length(ind) != 1 )
      {
        if ( length(ind) > 1)
        {
          warning("Node \"",node$name,"\": multiple occurrences of attribute \"",
                  n$attrib,"\"",sep="")
        }
        else
        {
          warning("Node \"",node$name,"\": attribute \"",n$attrib,"\" not found",sep="")
        }
        return(NA)
      }
      a <- attrib[ind]
    }
  }
  
  # evaluate results:
  
  if ( !is.numeric(a) )
  {
    if ( is.factor(a) ) a <- as.numeric(as.character(a))
    else                a <- as.numeric(a)
  }
  
  u <- approx(x=n$x,y=n$u,xout=a)$y
  ind.out.of.range <- (a < n$range[1]) | (a > n$range[2])
  u <- ifelse(ind.out.of.range,NA,u)
  if ( sum(ind.out.of.range,na.rm=T) > 0 )
  {
    ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
    warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib,"\" out of range: ",
            paste(a[ind.not.na],collapse=","),sep="")
  }
  
  # return results:
  
  return(u)
}


# print:
# ------

print.utility.endnode.intpol1d <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}

# summary:
# --------

summary.utility.endnode.intpol1d <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  cat("attribute:      ",node$attrib,"\n")
  cat("attribute range:",node$range[1],"-",node$range[2],"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("data pairs:","\n")
  names.x <- ifelse(is.na(node$names.x),"",node$names.x)
  names.u <- ifelse(is.na(node$names.u),"",node$names.u)
  print(data.frame(names.x=names.x,x=node$x,u=node$u,names.u=names.u))
}


# plot:
# -----

plot.utility.endnode.intpol1d <- 
  function(x,
           par       = NA,
           col       = utility.calc.colors(),
           gridlines = c(0.2,0.4,0.6,0.8),
           main      = "",
           cex.main  = 1,
           ...)
  {
    node <- x
    n <- updatepar(node,par)
    utility.endnode.plot1d(node      = n,
                           col       = col,
                           gridlines = gridlines,
                           main      = main,
                           cex.main  = cex.main,
                           ...)
    points(n$x,n$u,cex=1.5,xpd=TRUE) 
  }

