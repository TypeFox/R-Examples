################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# conversion node from values to utilities with parametric function: 
# class "utility.conversion.parfun"
# ==============================================================================


# constructor:
# ------------

utility.conversion.parfun.create <- function(name.node,    # character(1)
                                             node,         # node
                                             name.fun,     # name of f(a,par)
                                             par,          # numeric(n)
                                             names.par    = rep(NA,length(par)),
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
  if ( ! utility.check.name(name.node,list(node)) )
  {
    cat("*** Warning: node with same name \"",name.node,"\" exists already ",
        "as sub-node","\n")
    check.ok <- F
  }
  if ( ! check.ok )
  {
    cat("*** Warning: node \"",name.node,"\" could not be constructed","\n",
        sep="")
    return(NA)
  }
  
  # construct class:
  
  n <- list()
  n$name         <- name.node
  n$description  <- "utility/value parametric function conversion node"
  n$type         <- "utility.conversion.parfun"
  n$nodes        <- list(node)
  n$name.fun     <- name.fun
  n$par          <- par
  n$names.par    <- names.par
  n$required     <- required
  n$num.required <- 1
  n$utility      <- TRUE
  n$col          <- col
  n$shift.levels <- shift.levels
  class(n)       <- "utility.conversion.parfun" 
  
  # print and return class
  
  #cat(n$description," \"",name.node,"\" constructed","\n",sep="")   
  return(n)
}


# update parameter values:
# ------------------------

updatepar.utility.conversion.parfun <- function(x,par=NA,...)
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
        warning("Node \"",node$name,"\": multiple occurrences of parameter",
                names(par)[ind[1]],sep="")
        ind <- ind[1]
      }
      if ( length(ind) == 1 )
      {
        n$par[i] <- par[ind]
      }
    } 
  }
  n$nodes[[1]] <- updatepar(n$nodes[[1]],par)
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.cond.utility.conversion.parfun <- function(x,v,...)
{
  node <- x
  u <- do.call(node$name.fun,list(v,node$par))
  return(u)
}


evaluate.utility.conversion.parfun <- function(x,
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
  u <- as.data.frame(u)
  names(u) <- n$name
  ind <- !is.na(u) & (u<0 | u>1)
  if ( sum(ind) > 0 )
  {
    warning("Node \"",node$name,"\": node \"",n$name,"\" produced values outside of [0,1]: ",
            paste(u[ind],collapse=","),sep="")
  }
  
  # return results:
  
  u <- cbind(u,v)
  rownames(u) <- rownames(attrib)
  
  # return results:
  
  return(u)
}


# print:
# -----

print.utility.conversion.parfun <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.conversion.parfun <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  cat("node     :      ",node$nodes[[1]]$name,"\n")
  cat("function type:  ","utility","\n")
  cat("required:       ",node$required,"\n")
  cat("function:       ",node$name.fun,"\n")
  cat("parameters:","\n")
  names.par <- ifelse(is.na(node$names.par),"",node$names.par)
  print(data.frame(names.par=names.par,par=node$par))
  for ( i in 1:length(node$nodes) ) 
  {
    cat("***","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.conversion.parfun <- 
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

