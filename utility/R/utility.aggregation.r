################################################################################
#                                                                              #
# utility and value function package                                           #
# ==================================                                           #
#                                                                              #
# version 1.3                                        Peter Reichert 05.10.2014 #
#                                                                              #
################################################################################


# ==============================================================================
# utility node for (potentially) aggregating utility and/or end nodes: 
# class "utility.aggregation"
# ==============================================================================


# constructor:
# ------------

utility.aggregation.create <- 
  function(name.node,          # character(1)
           nodes,              # list of nodes
           name.fun,           # name of aggreg. fun f(u,par)
           par,                # numeric(n)
           names.par    = rep(NA,length(par)),
           required     = FALSE,
           num.required = 1,
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
    utility <- nodes[[1]]$utility
    if ( length(nodes) > 1 )
    {
      for ( i in 2:length(nodes) )
      {
        if ( nodes[[i]]$utility != utility )
        {
          cat("*** Warning: Mixted value and utility nodes",
              "cannot be aggregated","\n")
          check.ok <- F
        }
      }
    }
    if ( ! utility.check.name(name.node,nodes) )
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
    
    node <- list()
    node$name         <- name.node
    node$description  <- "utility/value aggregation node"
    node$type         <- "aggregationnode"
    node$nodes        <- nodes
    node$name.fun     <- name.fun
    node$par          <- par
    node$names.par    <- names.par
    node$required     <- required
    node$num.required <- num.required
    node$utility      <- utility
    node$col          <- col
    node$shift.levels <- shift.levels
    class(node)       <- "utility.aggregation" 
    
    # return class
    
    #cat(node$description," \"",name.node,"\" constructed","\n",sep="")   
    return(node)
  }


# update parameter values:
# ------------------------

updatepar.utility.aggregation <- function(x,par=NA,...)
{
  node <- x
  
  # check availabiliy of named parameter vector:
  
  if ( length(names(par)) == 0 ) return(node)
  
  # update conditional nodes:
  
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
  for ( i in 1:length(n$nodes) )
  {
    n$nodes[[i]] <- updatepar(n$nodes[[i]],par)
  }
  
  # return updated node:
  
  return(n)      
}


# evaluate values or utilities:
# -----------------------------

evaluate.utility.aggregation <- function(x,
                                         attrib,   # data.frame
                                         par=NA,
                                         ...)
{
  node <- x
  
  # check input:
  
  if ( ! is.data.frame(attrib) )
  {
    warning("Node \"",node$name,"\": attrib must be a data frame",sep="")
    return(NA)
  }
  
  # update parameters:
  
  n <- updatepar(node,par)
  
  # evaluate nodes:
  
  u <- evaluate(n$nodes[[1]],attrib)
  ind <- !is.na(u) & (u<0 | u>1)
  if ( sum(ind) > 0 )
  {
    warning("Node \"",node$name,"\": node \"",n$nodes[[1]]$name,"\" produced values outside [0,1]: ",
            paste(u[ind],collapse=","),sep="")
  }
  if ( ! is.data.frame(u) )
  {
    u <- as.data.frame(u)
    names(u) <- n$nodes[[1]]$name
  }
  required <- n$nodes[[1]]$required
  nodenames <- n$nodes[[1]]$name
  if ( length(n$nodes) > 1 )
  {
    for ( i in 2:length(n$nodes) )
    {
      u.i <- evaluate(n$nodes[[i]],attrib)
      ind <- !is.na(u) & (u<0 | u>1)
      if ( sum(ind) > 0 )
      {
        warning("Node \"",node$name,"\": node \"",n$nodes[[i]]$name,"\" produced values outside [0,1]: ",
                paste(u.i[ind],collapse=","),sep="")
      }
      if ( ! is.data.frame(u.i) )
      {
        u.i <- as.data.frame(u.i)
        names(u.i) <- n$nodes[[i]]$name
      }
      u <- cbind(u,u.i)
      nodenames[i] <- n$nodes[[i]]$name
      required[i]  <- n$nodes[[i]]$required 
    }
  }
  if ( length(unique(nodenames)) != length(nodenames) )
  {
    warning("Node \"",node$name,"\": node names are not unique:",
            paste(nodenames,collapse=","))
    u.agg <- as.data.frame(rep(NA,nrow(attrib)))
    names(u.agg) <- n$name
    u <- cbind(u.agg,u)
    rownames(u) <- rownames(attrib)
    return(u)
  }   
  
  # return results:
  
  u.agg.input <- as.matrix(u[,nodenames])
  u.agg <- apply(u.agg.input,1,n$name.fun,n$par)
  res.ok <- apply(u.agg.input,1,utility.check.required,
                  required,n$num.required)
  u.agg <- ifelse(res.ok,u.agg,NA)
  u.agg <- as.data.frame(u.agg)
  names(u.agg) <- n$name
  ind <- !is.na(u.agg) & (u.agg<0 | u.agg>1)
  if ( sum(ind)  > 0 )
  {
    warning("Node \"",node$name,"\": aggregation technique \"",n$name.fun,"\" produced values outside of [0,1]: ",
            paste(u.agg[ind],collapse=","),sep="")
  }
  u <- cbind(u.agg,u)
  rownames(u) <- rownames(attrib)
  
  return(u)
}


# print:
# -----

print.utility.aggregation <- function(x,...)
{
  cat(paste(rep("-",50),collapse=""),"\n")
  summary(x,...)
  cat(paste(rep("-",50),collapse=""),"\n")
}


# summary:
# --------

summary.utility.aggregation <- function(object,...)
{
  node <- object
  cat(node$name,"\n")
  cat(paste(rep("-",nchar(node$name)),collapse=""),"\n")
  cat(node$description,"\n")
  for ( i in 1:length(node$nodes) )
  {
    string1 <- "nodes:          "
    if ( i > 1 ) string1 <- "                "
    string2 <- node$nodes[[i]]$name
    if ( node$nodes[[i]]$type == "endnode" ) 
    {
      num.space <- max(1,15-nchar(node$nodes[[i]]$name))
      string2 <- paste(string2,
                       paste(rep(" ",num.space),collapse=""),
                       "(end node)",sep="") 
    }     
    cat(string1,string2,"\n")
  }
  cat("function:       ",node$name.fun,"\n")
  names.par <- ifelse(is.na(node$names.par),"",node$names.par)
  cat("parameters:","\n")
  print(data.frame(names.par=names.par,par=node$par))
  funtype <- "utility"; if ( !node$utility ) funtype <- "value"
  cat("function type:  ",funtype,"\n")
  cat("required:       ",node$required,"\n")
  cat("required nodes: ",node$num.required,"\n")
  for ( i in 1:length(node$nodes) ) 
  {
    cat("***","\n")
    summary(node$nodes[[i]])
  }
}


# plot:
# -----

plot.utility.aggregation <- 
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


# ==============================================================================
# utility aggregation functions
# ==============================================================================


utility.aggregate.add <- function(u,par)  # par[i]: weight of u[i]
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of additive aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  u.agg <- sum(par[ind]*u[ind])/s
  
  return(as.numeric(u.agg))
}


utility.aggregate.min <- function(u,par=NA)
{
  # check input:
  
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  
  # calculate aggregated value
  
  u.agg <- min(u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.max <- function(u,par=NA)
{
  # check input:
  
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  
  # calculate aggregated value
  
  u.agg <- max(u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.mult <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( length(ind) == 1 )
  {
    return(as.numeric(u[ind]))
  }
  if ( sum( par < 0 | par > 1 ) > 0 )
  {
    warning("Parameter of multiplicative aggregation",
            "smaller than zero or larger than unity")
    return(NA)
  }
  
  # function used in uniroot to determine the scaling constant k:
  
  utility.aggregate.mult.root <- function(k,ki)
  {
    res <- 1
    for ( i in 1:length(ki) )
    {
      res <- res * ( 1 + k * ki[i] )
    }
    res <- 1 + k - res
    return(res)
  }
  
  # define numerical parameter:
  
  eps <- 1e-3   # maximum deviation of sum(par) from unity to use additive fcn 
  
  # rescale weights:
  
  s <- sum(par)   
  fact <- s/sum(par[ind])
  ki <- fact*par[ind]
  
  # calculate additive utility function if sum close to unity:
  
  if ( s > 1-eps & s < 1+eps )
  {
    return(utility.aggregate.add(u,par))
  }
  
  # calculate multiplicative utility function if sum not close to unity:
  
  # calculate k: 
  # (Keeney and Raiffa, Decisions with multiple objectives, 1976,
  # pp. 307, 347-348)
  
  if ( s < 1 )
  {
    lower <- 1
    i <- 0
    while ( utility.aggregate.mult.root(lower,ki) < 0 )
    {
      lower <- 0.1*lower
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    upper <- 1
    i <- 0
    while ( utility.aggregate.mult.root(upper,ki) > 0 )
    {
      upper <- 10*upper
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    k <- uniroot(utility.aggregate.mult.root,ki=ki,
                 lower=lower,upper=upper)$root
  }
  else  # s > 1
  {
    upper <- -0.1
    i <- 0
    while ( utility.aggregate.mult.root(upper,ki) < 0 )
    {
      upper <- 0.1*upper
      i <- i+1
      if ( i > 20 )
      {
        warning("Problem solving equation for scaling constant")
        return(NA)
      }
    }
    k <- uniroot(utility.aggregate.mult.root,ki=ki,
                 lower=-1,upper=upper)$root 
  }
  
  # evaluate multiplicative utility function:
  
  u.agg <- 1  
  for ( i in 1:length(ki) )
  {
    if ( !is.na(u[ind][i]) ) u.agg <- u.agg * (k*ki[i]*u[ind][i]+1) 
  }
  u.agg <- (u.agg - 1)/k
  
  # eliminate values out of range due to numerical inaccuracies:
  
  u.agg <- ifelse(u.agg < 0, 0, u.agg)
  u.agg <- ifelse(u.agg > 1, 1, u.agg)
  
  return(as.numeric(u.agg))
}


utility.aggregate.geo <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of geometric aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  u.agg <- 1
  for ( i in 1:length(ind) )
  {
    if ( par[ind][i]>0 ) u.agg <- u.agg*u[ind][i]^(par[ind][i]/s)
  }
  
  return(as.numeric(u.agg))
}


utility.aggregate.revgeo <- function(u,par) 
{
  return(1-utility.aggregate.geo(1-u,par))
}


utility.aggregate.geooff <- function(u,par)
{
  n <- length(u)
  
  # check input:
  
  if ( length(par) != n + 1)
  {
    warning("Length of parameter vector should be length of utilities/values (for weights) plus one (for offset): ",
            length(par)," ",n)
    return(NA)
  }
  u <- utility.aggregate.geo(u+par[n+1],par[1:n])-par[n+1]
  # correct for numerical errors due to differences of "large" numbers
  u <- ifelse(u>0,u,0)
  u <- ifelse(u<1,u,1)
  return(u) 
}


utility.aggregate.revgeooff <- function(u,par) 
{
  return(1-utility.aggregate.geooff(1-u,par))
}


utility.aggregate.cobbdouglas <- function(u,par) 
{
  return(utility.aggregate.geo(u,par))
}


utility.aggregate.harmo <- function(u,par) 
{
  # check input:
  
  if ( length(u) != length(par) )
  {
    warning("Length of utilities/values and weights not equal: ",
            length(u)," ",length(par))
    return(NA)
  }
  ind <- which(!is.na(u))
  if ( length(ind) == 0 ) return(NA)
  if ( sum( par < 0 ) > 0 )
  {
    warning("Parameter of harmonic aggregation smaller than zero")
    return(NA)
  }
  
  # calculate aggregated value
  
  s <- sum(par[ind])
  if ( s <= 0 ) return(NA)
  if ( sum(u==0) > 0 ) return(0)
  
  u.agg <- s / sum(par[ind]/u[ind])
  
  return(as.numeric(u.agg))
}


utility.aggregate.revharmo <- function(u,par) 
{
  return(1-utility.aggregate.harmo(1-u,par))
}


utility.aggregate.harmooff <- function(u,par)
{
  n <- length(u)
  
  # check input:
  
  if ( length(par) != n + 1)
  {
    warning("Length of parameter vector should be length of utilities/values (for weights) plus one (for offset): ",
            length(par)," ",n)
    return(NA)
  }
  return(utility.aggregate.harmo(u+par[n+1],par[1:n])-par[n+1])
}


utility.aggregate.revharmooff <- function(u,par) 
{
  return(1-utility.aggregate.harmooff(1-u,par))
}


utility.aggregate.mix <- function(u,par)  # par[i]: weight of u[i]
{                                         # par[n+j]: weight of technique j
  # check input:                         # (j = add, min, geo)
  
  n <- length(u)
  if ( n+3 != length(par) )
  {
    warning("Length of parameter vector must be equal to",
            "length of utilities/values plus three:",
            length(par),length(u))
    return(NA)
  }
  s <- sum(par[n+(1:3)])
  if ( s <= 0 | sum(par[n+(1:3)]<0) > 0 )
  {
    warning("Weights of aggregation techniques to average",
            "cannot be negative or not all of them equal to zero")
    return(NA)
  }
  
  u.add <- 0; if ( par[n+1] != 0 ) u.add <- utility.aggregate.add(u,par[1:n])
  u.min <- 0; if ( par[n+2] != 0 ) u.min <- utility.aggregate.min(u)
  u.geo <- 0; if ( par[n+3] != 0 ) u.geo <- utility.aggregate.geo(u,par[1:n])
  
  if ( is.na(u.add) | is.na(u.min) | is.na(u.geo) ) return(NA)
  u.agg <- (par[n+1]*u.add + par[n+2]*u.min + par[n+3]*u.geo)/s
  
  return(u.agg)
}


utility.aggregate.addmin <- function(u,par)
{
  n <- length(u)
  if ( length(par) != n+1 )
  {
    warning("Length of parameter vector should be length of utilities/values ",
            "(for weights) plus one (for weight between methods): ", 
            length(par), " ", n)
    return(NA)
  }
  return(   par[n+1]  * utility.aggregate.add(u,par[1:n]) +
              (1-par[n+1]) * utility.aggregate.min(u,NA))
}
