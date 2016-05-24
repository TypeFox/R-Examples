#' Convert a bayou parameter list into a simmap formatted phylogeny
#' 
#' \code{pars2simmap} takes a list of parameters and converts it to simmap format
#' 
#' @param pars A list that contains \code{sb} (a vector of branches with shifts), \code{loc} (a vector of shift locations),
#' \code{t2} (a vector of theta indices indicating which theta is present after the shift).
#' @param tree A tree of class 'phylo'
#' 
#' @description This function converts a bayou formatted parameter list specifying regime locations into a simmap formatted tree that can
#' be plotted using \code{plotSimmap} from phytools or the \code{plotRegimes} function from bayou.
#' 
#' @return A list with elements: \code{tree} A simmap formatted tree, \code{pars} bayou formatted parameter list, and \code{cols} A named vector of colors.
#' 
#' @examples
#' tree <- reorder(sim.bdtree(n=100), "postorder")
#' 
#' pars <- list(k=5, sb=c(195, 196, 184, 138, 153), loc=rep(0, 5), t2=2:6)
#' tr <- pars2simmap(pars, tree)
#' plotRegimes(tr$tree, col=tr$col)
#' @export
pars2simmap <- function(pars,tree){
  tree <- reorder(tree, "postorder")
  sb <- pars$sb
  loc <- pars$loc
  t2 <- pars$t2
  if(!all(pars$sb %in% 1:nrow(tree$edge))) stop("Invalid parameter list. Specified branches not found in the tree")
  if(!all(pars$loc < tree$edge.length[pars$sb])) stop("Invalid parameter list. Some shift locations specified beyond the length of the branch")
  Th <- NULL
  nbranch <- length(tree$edge.length)
  maps <- lapply(tree$edge.length,function(x){y <- x; names(y) <- 1; y})
  dup <- which(duplicated(sb))
  if(length(dup)>0){
    maps[sb[-dup]] <- lapply(1:length(sb[-dup]),.addshift2map,maps=maps,sb=sb[-dup],loc=loc[-dup],t2=t2[-dup])
  } else {
    maps[sb] <- lapply(1:length(sb),.addshift2map,maps=maps,sb=sb,loc=loc,t2=t2)
  }
  for(i in dup){
    maps[[sb[i]]] <-.addshift2map(i,maps=maps,sb=sb,loc=loc,t2=t2)
  }
  nopt <- rep(1,nbranch)
  for(i in nbranch:1){
    if(i %in% sb){
      opt <- as.integer(names(maps[[i]])[length(maps[[i]])])
      nopt[tree$edge[i,2]] <- opt
      names(maps[[i]])[1] <- nopt[tree$edge[i,1]]
    } else {
      names(maps[[i]])[1] <- nopt[tree$edge[i,1]] 
      nopt[tree$edge[i,2]] <- nopt[tree$edge[i,1]]
    }
  }
  shiftdown <- nopt[tree$edge[,1]]
  new.maps <- lapply(1:nbranch,function(x){names(maps[[x]])[1] <- shiftdown[x]; maps[[x]]})
  new.maps <- maps
  for(j in 1:nbranch){
    names(new.maps[[j]])[1] <-shiftdown[j]
  }
  anc.theta <- unlist(lapply(new.maps[sb],function(x) as.integer(names(x)[length(x)-1])),F,F)
  o <- rev(order(sb,loc*-1))
  shifted.maps <- new.maps[sb[o]]
  t1 <- rep(NA,length(t2))
  for(i in 1:length(t2)){
    nm <- as.integer(names(maps[[sb[o][i]]]))
    t1[nm[2:length(nm)]-1] <- nm[1:(length(nm)-1)]
    Th[t2[o[i]]] <- Th[t1[o[i]]]
  }
  new.tree <- tree
  new.tree$maps <- new.maps
  new.pars <- pars
  col <- c(1,rainbow(pars$k))
  names(col) <- 1:(pars$k+1)
  return(list(tree=new.tree,pars=new.pars,col=col))
}

.pars2map <- function(pars, cache){
  nbranch <- length(cache$edge.length)
  nshifts <- table(pars$sb)
  shifts <- rep(0,nbranch)
  shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
  irow <- rep(1:nbranch,shifts+1)
  segs <- c(cache$edge.length, pars$loc)
  tmp.o <- c(1:nbranch, pars$sb)
  names(segs) <- tmp.o
  add.o <- order(tmp.o,segs)
  segs <- segs[add.o]
  ind <- names(segs)
  t2index <- add.o[which(add.o > nbranch)]
  t2b <- c(rep(1,length(segs)))
  t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
  loc.o <- order(pars$loc,decreasing=TRUE)
  sandwiches <- loc.o[duplicated(pars$sb[loc.o])]
  if(length(sandwiches)>0){
    sb.down <- pars$sb[-sandwiches]
    t2.down <- pars$t2[-sandwiches]
  } else {sb.down <- pars$sb; t2.down <- pars$t2}
  sb.o <- order(sb.down)
  sb.down <- sb.down[sb.o]
  t2.down <- t2.down[sb.o]
  sb.desc <- cache$bdesc[sb.down]
  desc.length <- unlist(lapply(sb.desc, length),F,F)
  sb.desc <- sb.desc[desc.length>0]
  names(t2b) <- names(segs)
  sb.desc2 <- unlist(sb.desc,F,F)
  sb.dup <- duplicated(sb.desc2)
  sb.desc3 <- sb.desc2[!sb.dup]
  t2.names <- rep(t2.down[desc.length>0], unlist(lapply(sb.desc,length),F,F))
  t2.names <- t2.names[!sb.dup]
  t2b[as.character(unlist(sb.desc3,F,F))] <- t2.names
  base <- duplicated(names(segs))*c(0,segs[1:(length(segs)-1)])
  segs <- segs-base
  #maps <- lapply(1:nbranch, function(x) segs[ind==x])
  #maps <- lapply(maps, function(x) if(length(x) >1) {c(x[1],diff(x[1:length(x)]))} else x)
  return(list(segs=segs,theta=t2b))
}

#' Calculates the alpha parameter from a QG model
#' 
#' @param pars A bayou formatted parameter list with parameters h2 (heritability), P (phenotypic variance) and w2 (width of adaptive landscape)
#' 
#' @return An alpha value according to the equation \code{alpha = h2*P/(P+w2+P)}. 
QG.alpha <- function(pars){
  pars$h2*pars$P/(pars$P+pars$w2*pars$P)
}

#' Calculates the sigma^2 parameter from a QG model
#' 
#' @param pars A bayou formatted parameter list with parameters h2 (heritability), P (phenotypic variance) and Ne (Effective population size)
#' 
#' @return An sig2 value according to the equation \code{alpha = h2*P/(Ne)}. 
QG.sig2 <- function(pars){
  (pars$h2*pars$P)/pars$Ne
}


#' Calculates the alpha and sigma^2 from a parameter list with supplied phylogenetic half-life and stationary variance
#' 
#' @param pars A bayou formatted parameter list with parameters halflife (phylogenetic halflife) and Vy (stationary variance)
#' 
#' @return A list with values for alpha and sig2.
OU.repar <- function(pars){
  alpha <- log(2)/pars$halflife
  sig2 <- (2*log(2)/(pars$halflife))*pars$Vy
  return(list(alpha=alpha,sig2=sig2))
}

.toSimmap <- function(map, cache){
  maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[names(map$segs)==x]; names(y) <- map$theta[names(map$theta)==x]; y })  
  tree <- cache$phy
  tree$maps <- maps
  return(tree)
}

#' Converts OUwie data into bayou format
#' 
#' \code{OUwie2bayou} Converts OUwie formatted data into a bayou formatted parameter list
#' 
#' @param tree A phylogenetic tree with states at internal nodes as node labels
#' @param trait A data frame in OUwie format
#' 
#' @return A bayou formatted parameter list
#' @export
OUwie2bayou <- function(tree, trait){
  tree <- reorder(tree, 'postorder')
  tip.states <- trait[,2]
  names(tip.states) <- trait[,1]
  states <- c(tip.states[tree$tip.label], tree$node.label)
  states <- unname(states)
  e1 <- states[tree$edge[,1]]
  e2 <- states[tree$edge[,2]]
  sb <- which(e1 != e2)
  loc <- 0.5*tree$edge.length[sb]
  t2 <- as.numeric(factor(e2[sb]))+1
  k <- length(sb)
  ntheta <- length(unique(t2))+1
  pars <- list(k=k, ntheta=ntheta, sb=sb, loc=loc, t2=t2)
  class(pars) <- c("bayoupars","list")
  return(pars)
}

#' Converts bayou data into OUwie format
#' 
#' \code{bayou2OUwie} Converts a bayou formatted parameter list into OUwie formatted tree and data table that can be analyzed in OUwie
#' 
#' @param pars A list with parameter values specifying \code{sb} = the branches with shifts,
#' \code{loc} = the location on branches where a shift occurs and \code{t2} = the optima to which
#' descendants of that shift inherit
#' @param tree A phylogenetic tree
#' @param dat A vector of tip states
#' 
#' @return A list with an OUwie formatted tree with mapped regimes and an OUwie formatted data table
#' @export
bayou2OUwie <- function(pars, tree, dat){
  if(is.null(names(dat))){
    warning("No labels on trait data, assuming the same order as the tip labels")
  } else {dat <- dat[tree$tip.label]}
  ntips <- length(tree$tip.label)
  cache <- .prepare.ou.univariate(tree, dat)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  tips <- which(tr$edge[,2] <= ntips)
  node.states <- sapply(tr$maps, function(x) names(x)[1])
  names(node.states) <- tr$edge[,1]
  node.states <- rev(node.states[unique(names(node.states))])
  tr$node.label <- as.numeric(node.states)
  tip.states <- sapply(tr$maps[tips], function(x) names(x)[length(x)])
  names(tip.states) <- tr$tip.label[tr$edge[tips,2]]
  tip.states <- as.numeric(tip.states[tr$tip.label])
  OUwie.dat <- data.frame("Genus_species"=tr$tip.label, "Reg"= tip.states, "X"= dat)
  rownames(OUwie.dat) <- NULL
  return(list(tree=tr, dat=OUwie.dat))
}
   
  
