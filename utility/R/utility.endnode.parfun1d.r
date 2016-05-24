################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for 1d (single attribute) parametric function: 
# class "utility.endnode.parfun1d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.parfun1d.create <- function(name.node,    # character(1)
                                            name.attrib,  # character(1)
                                            range,        # numeric(2)
                                            name.fun,     # name of f(a,par)
                                            par,          # numeric(n)
                                            names.par    = rep(NA,length(par)),
                                            utility      = TRUE,
                                            required     = FALSE,
                                            col          = "black",
                                            shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
  if ( length(par) != length(names.par) )
  {
    cat("*** Warning: par and names.par of different length:",
        length(par),length(names.par),"\n")
    check.ok <- F
  }
  if ( range[1] >= range[2] )
  {
    cat("*** Warning: Minimum of range not smaller than maximum:",
        range[1],range[2],"\n")
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
  node$name         <- name.node
  node$description  <- "utility/value 1d parametric function end node"
  node$type         <- "endnode"
  node$attrib       <- name.attrib
  node$range        <- range
  node$name.fun     <- name.fun
  node$par          <- par
  node$names.par    <- names.par
  node$required     <- required
  node$utility      <- utility
  node$col          <- col
  node$shift.levels <- shift.levels
  class(node)       <- "utility.endnode.parfun1d" 
  
  # print and return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.parfun1d <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update adequate values in interpolation list:
  
  n <- node
  for ( i in 1:length(n$par) )
  {
    if ( ! is.na(n$names.par[i]) )
    {
      ind <- which(n$names.par[i] == names(par) )
      if ( length(ind) > 1 )
      {
        warning("Node \"",node$name,"\": multiple occurrences of parameter \"",
                names(par)[ind[1]],"\"",sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$par[i] <- par[ind]
      }
    } 
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.parfun1d <- function(x,
                                              attrib,   # data.frame, numeric
                                              par = NA,
                                              ...)
{
  node <- x
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # extract attributes:
  
  if ( is.data.frame(attrib) )
  {
    if ( length(which(names(attrib)==n$attrib)) != 1 )
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
  
  u <- do.call(n$name.fun,list(a,n$par))
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
# -----

print.utility.endnode.parfun1d <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.parfun1d <- function(object,...)
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
  cat("function:       ",node$name.fun,"\n")
  cat("parameters:","\n")
  names.par <- ifelse(is.na(node$names.par),"",node$names.par)
  print(data.frame(names.par=names.par,par=node$par))
}


# plot:
# -----

plot.utility.endnode.parfun1d <- 
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
  }


# ==============================================================================
# simple parametric utility functions
# ==============================================================================


utility.fun.exp <- function(attrib,par)   # par[1]:  absolute risk aversion
{                                         # par[2]:  minimum of attribute range (default=0)
  # par[3]:  maximum of attribute range (default=1)
  atrans <- attrib
  if ( length(par) >= 3 ) atrans <- (attrib-par[2])/(par[3]-par[2])
  if ( par[1] == 0 ) return(atrans)
  return((1-exp(-atrans*par[1]))/(1-exp(-par[1])))
}

