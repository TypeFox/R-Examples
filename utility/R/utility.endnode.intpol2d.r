################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# endnode for interpolation based on isolines of two attributes: 
# class "utility.endnode.intpol2d"
# ==============================================================================


# constructor:
# ------------

utility.endnode.intpol2d.create <- function(name.node,   # character(1)
                                            name.attrib, # character(2)
                                            ranges,      # list(2) of numeric(2)
                                            isolines,    # list(n) of list of
                                            # x, y, and, optionally
                                            # names.x, names.y
                                            u,           # numeric(n)
                                            names.u      = rep(NA,length(u)),
                                            lead         = 0,
                                            utility      = TRUE,
                                            required     = FALSE,
                                            col          = "black",
                                            shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T 
  if ( length(name.attrib) != 2 )
  {
    cat("*** Warning: name.attrib must be of length 2","\n")
    check.ok <- F
  }
  if ( length(ranges) != 2 )
  {
    cat("*** Warning: ranges must be a list of two ranges","\n")
    check.ok <- F
  }
  else
  {
    if ( length(ranges[[1]]) != 2 )
    {
      cat("*** Warning: ranges[[1]] must contain two elements","\n")
      check.ok <- F
    }
    else
    {
      if ( ranges[[1]][1] >= ranges[[1]][2] )
      {
        cat("*** Warning: Minimum of range not smaller than maximum:",
            ranges[[1]][1],ranges[[1]][2],"\n")
        check.ok <- F
      }
    }
    if ( length(ranges[[2]]) != 2 )
    {
      cat("*** Warning: ranges[[2]] must contain two elements","\n")
      check.ok <- F
    }
    else
    {
      if ( ranges[[2]][1] >= ranges[[2]][2] )
      {
        cat("*** Warning: Minimum of range not smaller than maximum:",
            ranges[[2]][1],ranges[[2]][2],"\n")
        check.ok <- F
      }
    }
  } 
  if ( length(isolines) < 2 )
  {
    cat("*** Warning: at least two isolines are required","\n")
    check.ok <- F
  } 
  if ( length(isolines) != length(u) )
  {
    cat("*** Warning: isolines and u are of different length:",
        length(isolines),length(u),"\n")
    check.ok <- F
  }
  for ( i in 1:length(isolines) )
  {
    len.x <- length(isolines[[i]]$x) 
    if ( len.x < 2 )
    {
      cat("*** Warning: element x of isoline[[",i,"]] ",
          "must be of length > 1","\n",sep="")
      check.ok <- F
    }
    if ( len.x != length(isolines[[i]]$y)  )
    {
      cat("*** Warning: x and y in isoline[[",i,"]] ",
          "have different lengths:",
          len.x," ",length(isolines[[i]]$y),"\n",
          sep="")
      check.ok <- F
    }
    if ( length(isolines[[i]]$names.x) == 0 ) isolines[[i]]$names.x <- rep(NA,len.x) 
    if ( len.x != length(isolines[[i]]$names.x) ) 
    {
      cat("*** Warning: x and names.x in isoline[[",i,"]] ",
          "have different lengths:",
          len.x," ",length(isolines[[i]]$names.x),"\n",
          sep="")
      check.ok <- F
    }
    if ( length(isolines[[i]]$names.y) == 0 ) isolines[[i]]$names.y <- rep(NA,len.x) 
    if ( len.x != length(isolines[[i]]$names.y) ) 
    {
      cat("*** Warning: y and names.y in isoline[[",i,"]] ",
          "have different lengths:",
          len.x," ",length(isolines[[i]]$names.y),"\n",
          sep="")
      check.ok <- F
    }      
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
  node$description  <- "utility/value 2d interpolation end node"
  node$type         <- "endnode"
  node$attrib       <- name.attrib
  node$ranges       <- ranges
  node$isolines     <- isolines
  node$u            <- u
  node$names.u      <- names.u
  node$lead         <- lead
  node$required     <- required
  node$utility      <- utility
  node$col          <- col
  node$shift.levels <- shift.levels
  class(node)       <- "utility.endnode.intpol2d" 
  
  # print and return class
  
  #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
  return(node)
}


# update parameter values:
# ------------------------

updatepar.utility.endnode.intpol2d <- function(x,par=NA,...)
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
                names(par)[ind[1]],sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$u[i] <- par[ind]
      }
    } 
    for ( j in 1:length(n$isolines[[i]]$x) )
    {
      if ( ! is.na(n$isolines[[i]]$names.x[j]) )
      {
        ind <- which(n$isolines[[i]]$names.x[j] == names(par) )
        if ( length(ind) > 1 )
        {
          warning("Node \"",node$name,"\": multiple occurrences of parameter",
                  names(par)[ind[1]],sep="")
          ind <- ind[1]
        }
        if ( length(ind) == 1 )
        {
          n$isolines[[i]]$x[j] <- par[ind]
        }
      }
      if ( ! is.na(n$isolines[[i]]$names.y[j]) )
      {
        ind <- which(n$isolines[[i]]$names.y[j] == names(par) )
        if ( length(ind) > 1 )
        {
          warning("Node \"",node$name,"\": multiple occurrences of parameter",
                  names(par)[ind[1]],sep="")
          ind <- ind[1]
        }
        if ( length(ind) == 1 )
        {
          n$isolines[[i]]$y[j] <- par[ind]
        }
      }
    } 
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.endnode.intpol2d <- function(x,
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
    ind <- match(n$attrib,colnames(attrib))
    if ( sum(ifelse(is.na(ind),1,0)) > 0 )
    {
      warning("Node \"",node$name,"\": attribute(s) \"",
              paste(n$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
      return(rep(NA,nrow(attrib)))
    }
    a <- attrib[,ind]
  }
  else
  {
    if ( ! is.vector(attrib) )
    {
      warning("Node \"",node$name,"\": unknown format of attribute(s) \"",n$attrib,"\"",sep="")
      return(NA)
    }
    if ( length(names(attrib)) == 0 )
    {
      if ( length(attrib) == 2 )
        a <- as.matrix(attrib,nrow=1)
    }
    else
    {
      ind <- match(n$attrib,names(attrib))
      if ( sum(ifelse(is.na(ind),1,0)) > 0 )
      {
        warning("Node \"",node$name,"\": attribute(s) \"",
                paste(n$attrib[is.na(ind)],collapse=","),"\" not found",sep="")
        return(rep(NA,nrow(attrib)))
      }
      a <- as.matrix(attrib[ind],nrow=1)
    }
  }
  
  # evaluate results:
  
  if ( is.data.frame(a) )
  {
    if ( !is.numeric(a[,1]) )
    {
      if ( is.factor(a[,1]) ) a[,1] <- as.numeric(as.character(a[,1]))
      else                    a[,1] <- as.numeric(a[,1])
    }
    if ( !is.numeric(a[,2]) )
    {
      if ( is.factor(a[,2]) ) a[,2] <- as.numeric(as.character(a[,2]))
      else                    a[,2] <- as.numeric(a[,2])
    }
  }
  else
  {
    if ( !is.numeric(a) )
    {
      if ( is.factor(a) ) a <- as.numeric(as.character(a))
      else                a <- as.numeric(a)
    }
  }
  
  ind <- order(n$u)
  u <- utility.intpol2d(xy=a,isolines=n$isolines[ind],
                        levels=n$u[ind],lead=n$lead)
  
  ind.out.of.range <- (a[,1]<n$range[[1]][1])|(a[,1]>n$range[[1]][2])
  u <- ifelse(ind.out.of.range,NA,u)
  if ( sum(ind.out.of.range,na.rm=T) > 0 )
  {
    ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
    warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib[1],"\" out of range: ",
            paste(a[ind.not.na,1],collapse=","),sep="")
  }
  
  ind.out.of.range <- (a[,2]<n$range[[2]][1])|(a[,2]>n$range[[2]][2])
  u <- ifelse(ind.out.of.range,NA,u)
  if ( sum(ind.out.of.range,na.rm=T) > 0 )
  {
    ind.not.na <- ifelse(is.na(ind.out.of.range),F,ind.out.of.range)
    warning("Node \"",node$name,"\": value(s) of attribute \"",n$attrib[2],"\" out of range: ",
            paste(a[ind.not.na,2],collapse=","),sep="")
  }
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.endnode.intpol2d <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.endnode.intpol2d <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  cat("attributes:      ",paste(node$attrib,collapse=" , "),"\n")
  cat("attribute ranges:",node$range[[1]][1],"-",node$range[[1]][2],
      ",",node$range[[2]][1],"-",node$range[[2]][2],"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:   ",funtype,"\n")
  cat("required:        ",node$required,"\n")
  cat("isolines:","\n")
  for ( i in 1:length(node$u) )
  {
    name.u <- ""
    if ( !is.na(node$names.u[i]) ) 
    {
      name.u <- paste(":",node$names.u[i])
    }
    cat("u:",node$u[i],"  ",name.u,"\n")
    names.x <- rep("",length(node$isolines[[i]]$x))
    if ( length(node$isolines[[i]]$names.x) > 0 )
    {    
      names.x <- ifelse(is.na(node$isolines[[i]]$names.x),
                        "",node$isolines[[i]]$names.x)
    }
    names.y <- rep("",length(node$isolines[[i]]$y))
    if ( length(node$isolines[[i]]$names.y) > 0 )
    {    
      names.y <- ifelse(is.na(node$isolines[[i]]$names.y),
                        "",node$isolines[[i]]$names.y)
    }
    print(data.frame(names.x=names.x,
                     x=node$isolines[[i]]$x,
                     y=node$isolines[[i]]$y,
                     names.y=names.y))
  }
}


# plot:
# -----

plot.utility.endnode.intpol2d <- 
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
    utility.endnode.plot2d(node      = n,
                           col       = col,
                           gridlines = gridlines,
                           main      = main,
                           cex.main  = cex.main,
                           ...)
    ind <- order(n$u)
    levels <- n$u[ind]
    isolines <- n$isolines[ind]
    for ( i in 1:length(levels) )
    {
      lines(isolines[[i]],...)
      if ( i > 1 )
      {
        lines(c(isolines[[i-1]]$x[1],isolines[[i]]$x[1]),
              c(isolines[[i-1]]$y[1],isolines[[i]]$y[1]),
              ...)
        lines(c(isolines[[i-1]]$x[length(isolines[[i-1]]$x)],
                isolines[[i]]$x[length(isolines[[i]]$x)]),
              c(isolines[[i-1]]$y[length(isolines[[i-1]]$y)],
                isolines[[i]]$y[length(isolines[[i]]$x)]),
              ...)
      }
    }
  }

