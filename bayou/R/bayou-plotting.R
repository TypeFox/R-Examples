#' Utility for getting the starting and ending ages for each regime
.optima.ages <- function(pars,tree){
  nH <- nodeHeights(tree)
  reg <- sapply(tree$maps,function(x) names(x)[length(x)])
  adj <- sapply(tree$maps,function(x) ifelse(length(x)>1,x[1],0))
  abs.age <- nH[,1]+adj
  start <- tapply(abs.age,reg,min)
  end <- rep(max(nH),pars$ntheta)+1
  o <- as.character(1:pars$ntheta)
  names(end) <- names(start)
  return(cbind(start[o],end[o]))
}

#' Make a color transparent (Taken from an answer on StackOverflow by Nick Sabbe)
#' 
#' @param someColor A color, either a number, string or hexidecimal code
#' @param alpha The alpha transparency. The maxColorValue is set to 255.
#' 
#' @export
makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#' Plot a phylogenetic tree with posterior probabilities from a bayouMCMC chain (function adapted from phytools' plotSimmap)
#' 
#' @param chain A bayouMCMC chain
#' @param burnin The proportion of runs to be discarded, if NULL, then the value stored in the bayouMCMC chain's attributes is used
#' @param lwd The width of the edges
#' @param edge.type Either "theta" (branches will be colored according to their median value of theta), "regimes" (clades will be assigned to distinct regimes if the posterior probability of a shift
#' on that branch is > pp.cutoff), or "pp" (branches will be colored according to the probability of a shift on that branch). If "none" then edge.color will be assigned to all branches.
#' @param pal A color palette function used to paint the branches (unless edge.type="none")
#' @param pp.cutoff If edge.type=="regimes", the posterior probability above which a shift should be reconstructed on the tree.
#' @param circles a logical value indicating whether or not a circle should be plotted at the base of the node with values that correspond to the posterior probability of having a shift.
#' @param circle.cex.max The cex value of a circle with a posterior probability of 1
#' @param circle.col The color used to fill the circles
#' @param circle.pch the type of symbol used to plot at the node to indicate posterior probability
#' @param circle.lwd the line width of the points plotted at the nodes
#' @param circle.alpha a value between 0 and 255 that indicates the transparency of the circles (255 is completely opaque).
#' @param pp.labels a logical indicating whether the posterior probability for each branch should be printed above the branch
#' @param pp.col The color used for the posterior probability labels
#' @param pp.alpha a logical or numeric value indicating transparency of posterior probability labels. If TRUE, then transparency is ramped from invisible (pp=0), to black (pp=1). If numeric, all labels are given the same transparency. If NULL, then no transparency is given. 
#' @param pp.cex the size of the posterior probability labels 
#' @param edge.color The color of edges if edge.type="none"
#' @param parameter.sample When edge.type=="theta", the number of samples used to estimate the median "theta" value from each branch. Since this is 
#' computationally intensive, this enables you to downsample the chain.
#' @param ... Additional arguments passed to ape's plot.phylo
#' 
#' @export

plotSimmap.mcmc <- function(chain, burnin=NULL, lwd=1, edge.type = c("theta", "none", "regimes", "pp"), 
                            pal=rainbow, pp.cutoff=0.3, circles=TRUE, circle.cex.max=3, circle.col="red",
                            circle.pch=21, circle.lwd=0.75, circle.alpha=100, pp.labels=FALSE, pp.col=1, 
                            pp.alpha=255, pp.cex=0.75, edge.color = 1, parameter.sample=1000, ...){
  tree <- attributes(chain)$tree
  edge.type <- match.arg(edge.type, c("theta", "none", "regimes", "pp"))
  cache <- .prepare.ou.univariate(tree, attributes(chain)$dat)
  tree <- cache$phy
  if(is.null(burnin)) burnin = attributes(chain)$burnin
  if(is.null(burnin)) burnin = 0
  if(burnin==0) postburn <- 1:length(chain$gen)  else {
    postburn <- round(burnin*length(chain$gen),0):length(chain$gen)
  }
  L <- Lposterior(chain, tree, burnin=burnin)
  if(!is.null(pp.cutoff)){
    pp <- L$pp
    pars <- list()
    pars$sb <- which(pp > pp.cutoff)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb)+1
    pars$loc <- L$rel.location[pars$sb]*tree$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb)+1)
    if(length(pars$sb)>0){
      tr <- pars2simmap(pars, tree)$tree
      colors <- NULL
    } else {
      tr <- tree
      colors <- NULL
    }
  } else {
    tr <- tree
    tr$maps <- lapply(tr$edge.length, function(x) setNames(x, 1))
    colors <- setNames(1, 1)
  }
  .colorRamp <- function(trait, .pal, nn){
    strait <- (trait-min(trait))/(max(trait-min(trait)))
    itrait <- floor(strait*(nn-1))+1
    if(!is.null(.pal)){
    return(.pal(nn+1)[itrait])
    } else {
      return(itrait)
    }
  }
  if(edge.type=="none"){
    plot(tr, edge.color=edge.color, lwd=lwd, ...)
  }
  if(edge.type == "regimes"){
    plotRegimes(tr, col=colors, lwd=lwd, pal=pal, ...)
  }
  if(edge.type == "theta"){
    .ancestorBranches <- function(branch, cache){
      ancbranches <- which(sapply(cache$bdesc, function(x) branch %in% x))
      sort(ancbranches, decreasing=FALSE)
    }
    .branchRegime <- function(branch, abranches, chain, parameter, seqx, summary=FALSE){
      ancs <- c(branch, abranches[[branch]])
      ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
      ancshifts <- sapply(ancshifts, function(x) ifelse(length(x)==0, 1, x))
      ests <- sapply(1:length(ancshifts), function(x) chain[[parameter]][[seqx[x]]][ancshifts[x]])
      res <- cbind(ests)
      if(summary){
        return(apply(res, 2, median))
      } else {
        return(res)
      }
    }
    if(length(postburn) < parameter.sample){
      warning("Length of post-burnin sample less than the requested parameter sample, using entire post-burnin chain instead")
      seq1 <- postburn
    } else {
      seq1 <- sample(postburn, parameter.sample, replace=FALSE)
    }
    abranches <- lapply(1:nrow(tree$edge), .ancestorBranches, cache=cache)
    allbranches <- suppressWarnings(sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, "theta", seq1, summary=TRUE)))
    plot(tree, edge.color=.colorRamp(allbranches, pal, 100), ...)
  }
  if(edge.type == "pp"){
   plot(tree, edge.color=.colorRamp(L$pp, pal, 100), ...)
  }
  if(circles){
    #theta2 <- L$magnitude.of.theta2
    #root.median <- median(sapply(chain$theta[postburn], function(x) x[1]))
    #theta2[is.na(theta2)] <- root.median
    #theta2 <- theta2 - root.median
    #circle.cols <- sapply(colorRamp(theta2, circle.pal, 100), function(x) makeTransparent(x, circle.alpha))
    circle.cexs <- seq(0, circle.cex.max, length.out=100)[.colorRamp(L$pp, NULL, 100)]
    edgelabels(pch=circle.pch, lwd=circle.lwd, bg=makeTransparent(circle.col, circle.alpha), cex=circle.cexs)
  }
  if(pp.labels){
    edgelabels(round(L$pp,2), col=makeTransparent(pp.col, pp.alpha), cex=pp.cex, frame = "none")
  }
    
  
}

#' Adds visualization of regimes to a plot
#' 
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param cols A vector of colors to give to regimes, in the same order as pars$sb
#' @param type Either "rect", "density" or "lines". "rect" plots a rectangle for the 95\% CI for the stationary
#' distribution of a regime. "density" varies the transparency of the rectangles according to the probability density
#' from the stationary distribution. "lines" plots lines for the mean and 95\% CI's without filling them. 
#' @param transparency The alpha transparency value for the maximum density, max value is 255.
regime.plot <- function(pars,tree,cols,type='rect',transparency=100){
  OA <- .optima.ages(pars,tree)
  CIU95 <- pars$theta+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$theta-2*sqrt(pars$sig2/(2*pars$alpha))
  if(type=="lines"){
    for(i in 1:pars$ntheta){
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=makeTransparent(cols[i],transparency),lwd=3)
      lines(c(OA[i,1],OA[i,2]),rep(CIU95[i],2),col=makeTransparent(cols[i],transparency),lwd=1.25,lty=2)
      lines(c(OA[i,1],OA[i,2]),rep(CIL95[i],2),col=makeTransparent(cols[i],transparency),lwd=1.25,lty=2)
    }
  }
  if(type=="rect"){
    for(i in 1:pars$ntheta){
      rect(OA[i,1],CIL95[i],OA[i,2],CIU95[i],col=makeTransparent(cols[i],transparency),border=NA)
    }
  }
  if(type=="density"){
    ylim <- par('usr')[3:4]
    for(i in 1:pars$ntheta){
      x <- seq(OA[i,1],OA[i,2],length=10)
      y <- seq(ylim[1],ylim[2],length=100)
      Z <- matrix(nrow=length(x),ncol=length(y))
      for(j in 1:length(x)){
        Z[j,] <- dnorm(y,pars$theta[i],sqrt(pars$sig2/(2*pars$alpha)))
      }
      if(sum(Z)!=0){
        densregion(x,y,Z,colmax=makeTransparent(cols[i],transparency),colmin="transparent")
      }
      lines(c(OA[i,1],OA[i,2]),rep(pars$theta[i],2),col=makeTransparent(cols[i],min(255, 50+(transparency))),lwd=2)
    }
  }
}

#' Plot a pheongram with the posterior density for optima values
#' 
#' Plots a phenogram and the posterior density for optima values
#' 
#' @param tree A phylogeny of class 'phylo'
#' @param dat A named vector of tip data
#' @param burnin The initial proportion of the MCMC to be discarded
#' @param chain A bayouMCMC object that contains the results of an MCMC chain
#' @param colors An optional named vector of colors to assign to regimes, \code{NULL} results in no regimes being plotted.
#' @param pp.cutoff The posterior probability cutoff value. Branches with posterior probabilities of having a shift above this value
#' will have the average location of the regime shift painted onto the branches. 
#' @param K A list with the values of K to be plotted. If \code{NULL} all values of K are combined and a total posterior produced. This 
#' allows separate lines to be plotted for different numbers of shifts so that the location of optima can be compared, for example, between 
#' all samples that have 1 vs. 2 shifts in the posterior.
#' @param ... Additional parameters passed to \code{phenogram(...)}
#' 
#' @export
phenogram.density <- function(tree, dat, burnin=0, chain ,colors=NULL, pp.cutoff=NULL, K=NULL, ...){
  tree <- reorder(tree,"postorder")
  dat <- dat[tree$tip.label]
  postburn <- round(length(chain$gen)*burnin,0):length(chain$gen)
  chain2 <- lapply(chain,function(x) x[postburn])
  theta <- chain2$theta
  no.theta <- lapply(theta,length)
  min.theta <- min(unlist(theta))
  max.theta <- max(unlist(theta))
  if(is.null(K)){
    K <- list(unique(unlist(no.theta)))
  }
  if(!is.null(pp.cutoff)){
    L <- Lposterior(chain2, tree)
    pp <- L$pp
    pars <- list()
    pars$sb <- which(pp > pp.cutoff)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb)+1
    pars$loc <- L$rel.location[pars$sb]*tree$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb)+1)
    if(length(pars$sb)>0){
      tr <- pars2simmap(pars, tree)
      tree <- tr$tree
      colors <- tr$col
      names(colors) <- 1:length(colors)
    } else {
      tr <- tree
      colors <- 1; names(colors) <-1
    }
  }
  if(is.null(colors)){
      ntheta <- length(unique(names(unlist(tree$maps))))
      colors <- rainbow(ntheta)
      names(colors) <- 1:ntheta
    }
  nH <- max(nodeHeights(tree))
  plot(c(0,nH+0.3*nH),c(min(dat)-0.25,max(dat)+0.25),type='n',xlab="Time",ylab="Phenotype")
  phenogram(tree, dat, add=TRUE, colors=colors, spread.labels=FALSE, ...)
  dens.theta <- lapply(1:length(K), function(x) density(unlist(theta[no.theta %in% K[[x]]])))
  tmp <- sapply(1:length(dens.theta),function(Q){lines(nH+dens.theta[[Q]]$y*(0.3*nH)/max(dens.theta[[Q]]$y),dens.theta[[Q]]$x,col=Q+1)})
}

#' S3 method for plotting bayouMCMC objects
#' 
#' @param x A mcmc chain of class 'bayouMCMC' produced by the function bayou.mcmc and loaded into the environment using load.bayou
#' @param ... Additional arguments passed to \code{plot.mcmc} from the \code{coda} package
#' 
#' @export
#' @method plot bayouMCMC
plot.bayouMCMC <- function(x, ...){
  if(is.null(attributes(x)$burnin)){
    start <- 1
  } else {
    start <- round(attributes(x)$burnin*length(x$gen),0)
  }
  postburn <- start:length(x$gen)
  chain2 <- lapply(x,function(x) x[postburn])
  chain.length <- length(chain2$gen)
  univariates <- chain2[sapply(chain2,function(x) length(unlist(x)))==length(chain2$gen)]
  univariates$root <- sapply(chain2$theta, function(x) x[1])
  uni.df <- as.data.frame(univariates)
  rownames(uni.df) <- uni.df[,1]
  uni.df <- uni.df[,-1]
  plot(mcmc(uni.df), ...)
}


#' Plot parameter list as a simmap tree
#' 
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param ... Additional arguments passed to plotRegimes
#' 
#' @export
plotBayoupars <- function(pars, tree,...){
  mar <- par()$mar
  tree <- reorder(tree, 'postorder')
  X <- rep(0, length(tree$tip.label))
  names(X) <- tree$tip.label
  cache <- .prepare.ou.univariate(tree, X)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  plotRegimes(tr,...)
  par(mar=mar)
}

#' Experimental function for ancestral state reconstruction for a given OU model
.OU.asr <- function(tree, dat, pars, start=NULL, SE=0){
  phy <- reorder(tree, "postorder")
  dat <- dat[phy$tip.label]
  if(length(SE)>1){
    SE[phy$tip.label]
  }
  if(length(phy$tip.label) > 100) cat("This may take a while for large trees")
  EV <- .vcv.asrOU(phy, dat, pars, SE=SE)
  ntips <- length(phy$tip.label)
  ExpV <- EV$ExpV
  VCV <- EV$VCV
  diag(VCV) <- diag(VCV)+SE^2
  lik.fn <- function(anc){
    -1*dmnorm(as.vector(c(dat, anc)), mean=ExpV[,1], varcov=VCV, log=TRUE)
  }
  if(is.null(start)){
    start = ExpV[(length(dat)+1):(length(phy$edge.length)+1)]
  } 
  result <- optim(start, lik.fn, method="L-BFGS-B")
  x <- c(dat, result$par)
  names(x)[(ntips+1):length(x)] <- (ntips+1):length(x)
  return(x)
}  


.vcv.asrOU <- function(phy, dat, pars, SE, internal=TRUE){
  cache <- .prepare.ou.univariate(phy, dat, SE=SE)
  phy <- cache$phy
  new.pars <- pars
  ntips <- length(phy$tip.label)
  sig2 <- new.pars$sig2
  alpha <- new.pars$alpha
  D <- dist.nodes(phy)
  Cii <- D[ntips+1,]
  C <- D; C[,] <- 0
  ##Covariance[y_i, y_j]= s2/(2*alpha) * Exp[-alpha*t_ij] *  [1 - Exp(-2*alpha*t_a1)]
  for(i in 1:nrow(D)) for(j in 1:ncol(D))
    C[i,j]<- sig2/(2*alpha)*exp(-alpha*D[i,j])*(1-exp(-2*alpha*(Cii[i]+Cii[j]-D[i,j])))
  ##Calculate expectations
  diag(C) <- diag(C) + 1e-10
  mu <- rep(0, nrow(D))
  mu[1:ntips] <- dat
  W <- .allnodes.W(cache, new.pars)
  ExpV <- W %*% new.pars$theta
  return(list(ExpV=ExpV, VCV=C))
}

.allnodes.W <- function(tree, pars){
  a <- pars$alpha
  s2 <- pars$sig2
  nbranch <- length(tree$edge.length)
  if(class(tree)=="phylo"){
    X <- rep(NA,length(tree$tip.label))
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars$ntheta)){
    pars$ntheta <- length(pars$theta)
  }
  plook <- function(x){mapply(paste,x[2:length(x)],x[1:(length(x)-1)],sep=",")}
  tB <- cache$desc$anc[1:(cache$n.node+cache$ntips)]
  tB <- mapply(c,1:(cache$n.node+cache$ntips),tB, SIMPLIFY=FALSE)
  lookup <- lapply(tB,plook)
  edge.names <- mapply(paste,cache$edge[,1],cache$edge[,2],sep=",")
  cache$branchtrace <- t(sapply(lookup,function(x) as.numeric(edge.names %in% x)))
  smtree <- pars2simmap(pars, cache$phy)
  maps <- smtree$tree$maps
  allnodes <- cache$n.node+cache$ntips
  W <- matrix(0, ncol=pars$ntheta, allnodes)
  for(i in 1:allnodes){
    m <- maps[as.logical(cache$branchtrace[i,])]
    m <- c(0, rev(unlist(lapply(m, rev))))
    names(m)[1] <- 1
    TH <- sum(m)
    csm <- cumsum(m)
    eT <- exp(-a*TH)*(exp(a*csm[2:length(csm)])-exp(a*csm[1:(length(csm)-1)]))
    w <- tapply(eT, names(csm)[2:length(csm)], sum)
    W[i, as.numeric(names(w))] <- w
    W[i, 1] <- W[i,1] + exp(-a*TH)
  }
  return(W)
}

#' Experimental phenogram plotting function for set of model of model parameters
#' 
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param dat A named vector of tip data
#' @param regime.col A named vector of colors equal in length to the number of regimes
#' @param SE Standard error of the tip states
#' @param ... Optional arguments passed to \code{phenogram()}
#' 
#' @details This is an experimental plotting utility that can plot a phenogram with a given regime painting from
#' a parameter list. Note that it uses optimization of internal node states using matrix inversion, which is very 
#' slow for large trees. However, what is returned is the maximum likelihood estimate of the internal node states 
#' given the model, data and the parameter values.
#' 
#' @examples
#' \dontrun{
#' tree <- sim.bdtree(n=50)
#' tree$edge.length <- tree$edge.length/max(branching.times(tree))
#' prior <- make.prior(tree, dists=list(dk="cdpois", dsig2="dnorm", 
#'            dtheta="dnorm"), param=list(dk=list(lambda=5, kmax=10), 
#'              dsig2=list(mean=1, sd=0.01), dtheta=list(mean=0, sd=3)), 
#'                plot.prior=FALSE)
#' pars <- priorSim(prior, tree, plot=FALSE, nsim=1)$pars[[1]]
#' pars$alpha <- 4
#' dat <- dataSim(pars, model="OU", phenogram=FALSE, tree)$dat
#' OUphenogram(pars, tree, dat, ftype="off")
#' }
#' @export
OUphenogram <- function(pars, tree, dat, SE=0, regime.col=NULL, ...){
  datanc <- .OU.asr(tree, dat, pars, SE=SE)
  tr <- pars2simmap(pars, reorder(tree,"postorder"))
  OA <- .optima.ages(pars, tr$tree)
  CIU95 <- pars$theta+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$theta-2*sqrt(pars$sig2/(2*pars$alpha))
  if(is.null(regime.col)){
    regime.cols <- tr$col
  } else {regime.cols <- regime.col}
  ylim <- c(min(c(dat, pars$theta-2*sqrt(pars$sig2/(pars$alpha*2)))), max(c(dat, pars$theta-2*sqrt(pars$sig2/(pars$alpha*2)))))
  phenogram(tr$tree, datanc, colors=regime.cols, ylim=ylim, spread.labels=FALSE, ...)
  for(i in 1:pars$ntheta){
    x <- seq(OA[i,1],OA[i,2],length=10)
    y <- seq(ylim[1],ylim[2],length=100)
    Z <- matrix(nrow=length(x),ncol=length(y))
    for(j in 1:length(x)){
      Z[j,] <- dnorm(y,pars$theta[i],sqrt(pars$sig2/(2*pars$alpha)))
    }
    if(sum(Z)!=0){
      densregion(x,y,Z,colmax=makeTransparent(regime.cols[i]),colmin="transparent")
    }
    lines(c(OA[i,1],OA[i,2]),rep(pars$theta[i],2),col=regime.cols[i],lwd=2)
  }
  phenogram(tr$tree, datanc, , colors=regime.cols, add=TRUE, spread.labels=FALSE,  ...)
}

#' Function to plot the regimes from a simmap tree
#' 
#' @param tree A simmap tree of class phylo or simmap with a tree$maps list
#' @param col A named vector of colors to assign to character states, if NULL, then colors are generated from pal
#' @param lwd A numeric value indicating the width of the edges
#' @param pal A color palette function to generate colors if col=NULL
#' @param ... Optional arguments that are passed to plot.phylo
#' 
#' @details This function uses plot.phylo to generate coordinates and plot the tree, but plots the 
#' 'maps' element of phytools' simmap format. This provides much of the functionality of plot.phylo from
#' the ape package. Currently, only types 'phylogram', 'unrooted', 'radial', and 'cladogram' are allowed. Phylogenies must
#' have branch lengths.
#' 
#' @export
plotRegimes <- function(tree, col=NULL, lwd=1, pal=rainbow, ...){
  if(is.null(col)){
    regNames <- unique(names(unlist(tree$maps)))
    nreg <- length(regNames)
    col <- setNames(pal(nreg), regNames)
  }
  #nodecols <- col[sapply(tree$maps, function(x) names(x)[1])]
  tmp <- plot(tree, edge.color="#FFFFFF00", use.edge.length=TRUE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  #if(lastPP$type != "phylogram") stop("Currently only able to plot phylograms")
  nbranch <- nrow(tree$edge)
  .getBranchCoords <- function(i){
    xx <- lastPP$xx[tree$edge[i,]]
    yy <- lastPP$yy[tree$edge[i,]]
    xdist <- diff(xx)
    ydist <- diff(yy)
    map <- tree$maps[[i]]
    cs <- cumsum(c(0, map))/sum(map)
    colmap <- col[names(map)]
    return(list(xx=xx, yy=yy, xdist=xdist, ydist=ydist, cs=cs, colmap=colmap, nsegs=length(cs)-1, segreg = names(colmap)))
  }
  coords <- lapply(1:nbranch, .getBranchCoords)
  .phylogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    if(lastPP$direction %in% c("upwards", "downwards")){
      xcoord <- rbind(xx, matrix(xx[2], nrow=nsegs, ncol=2))
      ycoord <- rbind(rep(yy[1],2), cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
    if(lastPP$direction %in% c("leftwards", "rightwards")){
      ycoord <- rbind(yy, matrix(yy[2], nrow=nsegs, ncol=2))
      xcoord <- rbind(rep(xx[1],2), cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
  }
  .cladogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    xcoord <- cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1])
    ycoord <- cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1])
    rownames(xcoord) <- rownames(ycoord) <- segreg
    cols <- colmap
    dum <- lapply(1:nsegs, function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
  }
  .fanLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    circular.plot(lastPP$edge, lastPP$Ntip, lastPP$Nnode, lastPP$xx, lastPP$yy, )
  }
  if(lastPP$type=="fan") warning("type='fan' not currently supported, plotting a radial cladogram")
  plotfn <- switch(lastPP$type, phylogram=.phylogramLines, cladogram=.cladogramLines, unrooted=.cladogramLines, radial=.cladogramLines, fan=.cladogramLines)
  dum <- lapply(coords, plotfn)
}


