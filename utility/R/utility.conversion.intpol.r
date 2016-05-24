################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# conversion node from values to utilities with interpolation: 
# class "utility.conversion.intpol"
# ==============================================================================


# constructor:
# ------------

utility.conversion.intpol.create <- function(name.node,    # character(1)
                                             node,         # character(1)
                                             x,            # numeric(n)
                                             u,            # numeric(n)
                                             names.x      = rep(NA,length(x)),
                                             names.u      = rep(NA,length(u)),
                                             required     = FALSE,
                                             col          = "black",
                                             shift.levels = 0)
{
  # consistency checks:
  
  check.ok <- T   
  if ( length(x) != length(u) )
  {
    cat("*** Warning: x and u of different length:",
        length(x),length(u),"\n")
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
  if ( ! utility.check.name(name.node,node) )
  {
    cat("*** Warning: Node with same name \"",name.node,"\" exists already ",
        "as sub-node","\n")
    check.ok <- F
  }
  if ( ! check.ok )
  {
    cat("*** Warning: Node \"",name.node,"\" could not be constructed","\n",
        sep="")
    return(NA)
  }
  
  # construct class:
  
  n <- list()
  n$name         <- name.node
  n$description  <- "utility/value interpolation conversion node"
  n$type         <- "conversionnode"
  n$nodes        <- list(node)
  n$x            <- x
  n$u            <- u
  n$names.x      <- names.x
  n$names.u      <- names.u
  n$required     <- required
  n$num.required <- 1
  n$utility      <- TRUE
  n$col          <- col
  n$shift.levels <- shift.levels
  class(n)       <- "utility.conversion.intpol" 
  
  # print and return class
  
  #cat(n$description," \"",name.node,"\" constructed","\n",sep="")   
  return(n)
}


# update parameter values:
# ------------------------

updatepar.utility.conversion.intpol <- function(x,par=NA,...)
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
        warning("Node \"",node$name,"\": multiple occurrences of parameter",
                names(par)[ind[1]],sep="")
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
                names(par)[ind[1]],sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$u[i] <- par[ind]
      }
    } 
  }
  n$nodes[[1]] <- updatepar(n$nodes[[1]],par)
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.cond.utility.conversion.intpol <- function(x,v,...)
{
  node <- x
  u <- approx(x=node$x,y=node$u,xout=v)$y
  return(u)
}


evaluate.utility.conversion.intpol <- function(x,
                                               attrib,   # data.frame, numeric
                                               par = NA,
                                               ...)
{
  node <- x
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # evaluate results:
  
  v <- evaluate(n$nodes[[1]],attrib)
  if ( ! is.data.frame(v) )
  {
    v <- as.data.frame(v)
  }
  u <- evaluate.cond(n,v[,1])
  ind <- !is.na(u) & (u<0 | u>1)
  if ( sum(ind) > 0 )
  {
    warning("Node \"",node$name,"\": node \"",n$name,"\" produced values outside of [0,1]: ",
            paste(u[ind],collapse=","),sep="")
  }
  u <- as.data.frame(u)
  names(u) <- node$name
  
  # return results:
  
  u <- cbind(u,v)
  rownames(u) <- rownames(attrib)
  
  return(u)
}


# print:
# -----

print.utility.conversion.intpol <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.conversion.intpol <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ","utility","\n")
  cat("required:       ",node$required,"\n")
  cat("data pairs:","\n")
  names.x <- ifelse(is.na(node$names.x),"",node$names.x)
  names.u <- ifelse(is.na(node$names.u),"",node$names.u)
  print(data.frame(names.x=names.x,x=node$x,u=node$u,names.u=names.u))
  for ( i in 1:length(node$nodes) ) 
  {
    cat("***","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.conversion.intpol <- 
  function(x,
           u           = NA,
           uref        = NA,
           par         = NA,
           type        = c("hierarchy","table","node","nodes"),
           nodes       = NA,
           col         = utility.calc.colors(),
           gridlines   = c(0.2,0.4,0.6,0.8),
           main        = "",
           cex.main    = 1,
           cex.nodes   = 1,
           cex.attrib  = 1,
           f.reaches   = 0.2,
           f.nodes     = 0.2,
           with.attrib = TRUE,
           levels      = NA,
           plot.val    = TRUE,
           print.val   = TRUE,
           ...)
{
    node <- x
    n <- updatepar(node,par)
    utility.plot(node        = n,
                 u           = u,
                 uref        = uref,
                 type        = type,
                 nodes       = nodes,
                 col         = col,
                 gridlines   = gridlines,
                 main        = main,
                 cex.main    = cex.main,
                 cex.nodes   = cex.nodes,
                 cex.attrib  = cex.attrib,
                 f.reaches   = f.reaches,
                 f.nodes     = f.nodes,
                 with.attrib = with.attrib,
                 levels      = levels,
                 plot.val    = plot.val,
                 print.val   = print.val,
                 ...)
  }

