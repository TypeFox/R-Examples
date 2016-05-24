#' Null models for functional-phylogenetic diversity
#'
#' Simulate expectations (under a null model) of mean pairwise distance for 
#' a set of communities with different species richness.
#'
#' @details If \code{plot == TRUE}, then a surface is drawn giving the
#' null distribution.  Lighter shades of gray give larger intervals
#' with categories: 0.005-0.995 = 99\%, 0.025-0.975 = 95\%, 0.05-0.95
#' = 90\%, 0.25-0.75 = 50\%.
#' @param object a \code{\link{comparative.comm}} object, with
#' presence-absence community data.
#' @param type character string giving the type of distance matrix on
#' which the mean pairwise distance is based. Either "trait" or "phy"
#' to a phylogenetic or trait-based distance matrix, or an actual
#' matrix to use (e.g., one from \code{\link{funct.phylo.dist}})
#' @param n.sim The number of permutations of the presence vector used to 
#'  make the estimations.
#' @param plot TRUE or FALSE to make the plot of the expected average 
#'  mean pairwise distance, and the 5-95\% confidence interval.
#' @param disp99 Display the 99\% interval?
#' @note No serious checking of user-provided matrices is performed;
#' this is both useful and dangerous!
#' @seealso \code{\link{sim.phy}} \code{\link{scape}}
#' @return \code{matrix} with quantiles of mean pairwise distances for
#' all quantiles of of mean pairwise distance, with one row for the
#' range of species richnesses in the data (see column SpRich).
#' @author Steve Walker, wrappers by Will Pearse
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' #Must have all species present in at least one community!
#' #...and must be presence-absence data
#' data <- data[,colSums(data$comm) > 0]
#' data$comm[data$comm>1] <- 1
#' sims <- ConDivSim(data)
#' #...without traits...
#' sims.phy <- ConDivSim(data, type="phy")
#' @importFrom ape cophenetic.phylo
#' @importFrom stats quantile sd
#' @importFrom graphics polygon lines points text
#' @importFrom grDevices grey
#' @export
ConDivSim<-function(object, type="traits", n.sim=100, plot = TRUE, disp99 = FALSE){
    #Assertions and argument handling
    if(!inherits(object, "comparative.comm"))
        stop("'data' must be a comparative community ecology object")
    if(is.matrix(type)){
        if(nrow(type)==length(object$phy$tip.label) & nrow(type)==ncol(type))
            dist <- type else stop("'Provided distance matrix is of incorrect dimension'")
    } else {
        if(inherits(type, "dist")){
            dist <- as.matrix(type)
            if(nrow(dist)!=length(object$phy$tip.label) | nrow(dist)!=ncol(dist))
                stop("'Provided distance matrix is of incorrect dimension'")
        } else {
            switch(type,
                   "traits" = {
                       if(is.null(object$data)) 
                           stop("'object' must contain trait data if using a trait matrix")
                       dist <- as.matrix(dist(object$data))
                   },
                   "phy" = dist <- cophenetic.phylo(object$phy),
                   stop("Unknown 'type' of analysis"))
        }
    }
    object$comm <- round(object$comm)
    if(!all(object$comm %in% c(0L, 1L))) stop("Community data must be presence-absence")
    
    ## Calculate the observed species richness
    SpeRich.obs <- rowSums(object$comm)
    Occur.obs <- colSums(object$comm)
	
    ## Check that numbers are high enough for randomisation
    ## (e.g. no communities with one species only)
    if(any(SpeRich.obs < 2)) 
        stop('all communities must have more than one species present')
    if(any(Occur.obs < 1)) 
        stop('all species must be present at more than zero communities')
    
    ## Calculate the range of species richness in the communities
    if(min(SpeRich.obs)>2){
        SpRich <- (min(SpeRich.obs)-1):(max(SpeRich.obs)+1)
    }
    else {
        SpRich <- min(SpeRich.obs):(max(SpeRich.obs)+1)
    }

    ## Calculate the regional species richness(gamma)
    NbSp <- ncol(object$comm)
    
    ## Create the matrix in which to store the results
    Randomiz <- matrix(0, length(SpRich), 11)
    colnames(Randomiz)<-c("SpRich","ExpMeanMPD",
                          "ExpQ0.005MPD","ExpQ0.025MPD",
                          "ExpQ0.05MPD","ExpQ0.25MPD",
                          "ExpQ0.75MPD","ExpQ0.95MPD",
                          "ExpQ0.975MPD","ExpQ0.995MPD",
                          "ExpSdMPD")
    row.names(Randomiz) <- seq_along(SpRich)
    
    ## calculate the observed mean pairwise distance for the community
    ## (with the unweighted method from mpd in picante, it means only 
    ## the lower triangle!)
    MPD.obs <- .mpd(object, dist=as.dist(dist), abundance.weighted = FALSE)
    
    ## Loop on the species richness range to calculate the random 
    ## expectations of mean pairwise distance
    for (k in seq_along(SpRich)){
        
    	# Vector of local richness within the possible range SpRich
    	locrich <- SpRich[k]
        
    	# Create the basic vector of presence
        # browser()
        Vec.Pres <- c(rep(1, locrich), rep(0, NbSp-locrich))
        
        # Estimate the random expectations of mean pairwise distances 
        # for the random community with locrich as the species richness
        Sim.Div <- NULL
        
        for (i in 1:n.sim) {
            Vec.Pres.Sim <- sample(Vec.Pres)
            sample.dis <- dist[as.logical(Vec.Pres.Sim), as.logical(Vec.Pres.Sim)]
            
            # line in mpd from picante with abundance.weighted = FALSE
            Sim.Div[i] <- mean(sample.dis[lower.tri(sample.dis)])
      	}
        
        # Store the results in Randomiz
        Randomiz[k,1]<-locrich
        Randomiz[k,2]<-mean(Sim.Div)
        Randomiz[k, 3:10] <- quantile(Sim.Div, c(0.005, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975, 0.995))
    	Randomiz[k,11]<-sd(Sim.Div)
    }
	
    ## Make the graph if required
    if(plot == TRUE){
        
        plot(Randomiz[,1:2],
             ylim = c(min(Randomiz[, 3], MPD.obs), max(Randomiz[, 10], MPD.obs)),
             type="n", xlab="Species richness",ylab="Mean pairwise distance")
        
        ## Surface for each of the quantile pair:
        if(disp99){
            polygon( ## 0.005-0.995 = 99%
                    c(Randomiz[,1], Randomiz[length(SpRich):1, 1]), 
                    c(Randomiz[,3], Randomiz[length(SpRich):1, 10]),
                    col=grey(0.9), border=NA)
    	}
    	
        polygon( ## 0.025-0.975 = 95%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,4], Randomiz[length(SpRich):1,9]),
                col=grey(0.8),border=NA)
    	
        polygon( ## 0.05-0.95 = 90%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,5], Randomiz[length(SpRich):1,8]),
                col=grey(0.7),border=NA)
        
        polygon( ## 0.25-0.75 = 50%
                c(Randomiz[,1], Randomiz[length(SpRich):1,1]),
                c(Randomiz[,6], Randomiz[length(SpRich):1,7]),
                col=grey(0.6),border=NA)
    	
        ## Average of the expected mean pairwise distance
        lines(Randomiz[,1:2], col = "black", lwd = 2)
        
        ## Labelling the communities
        if(is.null(row.names(object$comm)))
            points(SpeRich.obs, MPD.obs) else
        text(SpeRich.obs, MPD.obs, labels = row.names(object$comm))
        
        
    }
    
    return(Randomiz)
}

#' Produces simulated communities based on species attributes
#' 
#' \code{trait.asm} calculates phylogenetic biodiversity metrics
#' 
#' @param a species attributes (e.g., traits like body size)
#' @param m number of communities to be simulated
#' @param meanSR target mean species richness across simulated communities
#' @param interval adjust to obtain \code{meanSR}
#' @param exponential use the exponential distribution when simulating communities
#' @param Pscale adjust this value when not using the exponential distribution in order to scale the species richnesses in the simulated communities
#' @details Simulates a set of communties based on the supplied attribute (trait) where larger values make it more likely for species to be in the communities.
#
#' @return \code{Y} presence/absence matrix
#' @return \code{P} probabilities
#' @return \code{a} the supplied trait
#' @return \code{exponential} if the exponential distribution was used
#' @return \code{meanSR} supplied \code{meanSR} value
#' @return \code{std} estimated sd
#' @author M.R. Helmus
#' @references Helmus M., Mercado-Silva N. & Vander Zanden M.J. (2013). Subsidies to predators, apparent competition and the phylogenetic structure of prey communities. Oecologia, 173, 997-1007.
#' @examples
#' \dontrun{
#'  data(laja)
#'  trait.asm(laja$fish.pref)
#' }
#' @importFrom stats rbinom optimize
#' @export
trait.asm<-function(a, m=1000,meanSR=NULL,interval=c(.001,10),exponential=TRUE,Pscale=1)
{
    a<-a[!is.na(a)]
	  n<-length(a)
    e<-as.matrix(a)
    
  if(!is.null(meanSR)){
    zz<-function(std,e,m,n,meanSR)
    {
      P <- exp(array(e,c(n,m))/std)/exp(max(e)/std)
      p <- array(P,c(n*m,1))
      Y <- rbinom(n*m, 1, p)
      Y <- t(array(Y, c(n,m)))
      abs(meanSR-mean(rowSums(Y)))
    }
    std<-unlist(optimize(zz,interval=interval,n=n,m=m,e=e,meanSR=meanSR))[1]    
  } else {std<-1}

  #exponential distribution
  if(exponential)
  {
    P <- (exp(array(e,c(n,m))/std)/exp(max(e)/std))
  } else {
    #do not use exponential distribution
    P <- Pscale*array(e,c(n,m))
  }
  p <- array(P,c(n*m,1))
	# convert distribution to presence/absence
	Y <- rbinom(n*m, 1, p)
	Y <- t(array(Y, c(n,m)))
  colnames(Y)<-names(a)
  return(list(Y=Y,P=P,a=a,exponential=exponential,meanSR=meanSR,std=std))
}

##' Traitgram for comparative community object
##'
##' \code{traitgram.cc} A wrapper for the
##' \code{\link[picante:traitgram]{traitgram}} function in the
##' \code{picante} package.
##'
##' @param object A \code{\link{comparative.comm}} object.
##' @param trait Which trait to plot.  If \code{\link{missing}}, use
##' the first trait.  If a positive \code{\link{numeric}} vector of
##' \code{\link{length}} one, use the \code{as.integer(trait)}th
##' trait.  If a \code{\link{numeric}} vector, use it instead of the
##' trait data in \code{object}.  If a \code{\link{character}} vector
##' of \code{\link{length}} one, use the trait with that name.  If a
##' \code{\link{function}} pass the trait data frame through that
##' function and use the result (\code{\link{princompOne}} is a
##' wrapper).  If an \code{\link{expression}}, evaluate that
##' expression in the environment of the trait data and use the
##' result.  If a \code{\link{character}} vector, then convert to an
##' expression and evaluate in the environment of the trait data and
##' use the result.
##' @param moreArgs List of more arguments to pass on to \code{trait}
##' (if its a \code{\link{function}}).
##' @param ... Additional arguments to be passed on to
##' \code{\link[picante:traitgram]{traitgram}} or \code{prcomp} for
##' \code{traitgram.cc} and \code{princomOne} respectively.
##' @return \code{traitgram.cc}: see
##' \code{\link[picante:traitgram]{traitgram}}
##' @importFrom picante traitgram
##' @importFrom ape multi2di
##' @export
traitgram.cc <- function(object, trait, moreArgs = NULL, ...) {
    if(is.null(object$data)) stop("must supply trait information")
    if(is.null(object$phy)) stop("must supply phylogeny")
    if(missing(trait)) {
        if(is.null(dim(object$data))) {
            tt <- object$data
        } else {
            tt <- object$data[, 1]
        }
    } else if(is.numeric(trait)) {
        if(length(trait) == 1) {
            if(trait < 1) stop("trait can't be a negative number")
            tt <- object$data[, as.integer(trait)]
        } else {
            tt <- trait
        }
    } else if(is.function(trait)) {
        tt <- do.call(trait, c(list(object$data), moreArgs))
    } else if(is.language(trait)) {
        tt <- with(object$data, eval(trait))
    } else if(is.character(trait)) {
        trait <- parse(text = paste(trait, collapse = "; "))
        tt <- with(object$data, eval(trait))
    }

    # need names for traitgram
    if(is.null(names(tt)))
        names(tt) <- rownames(object$data)

    # resolve possible polytomies
    pp <- multi2di(object$phy)

    # plot
    traitgram(tt, pp, ...)
}

##' First axis of a principal components analysis
##'
##' \code{princompOne} A very soft wrapper for \code{\link{princomp}}
##' 
##' @param x A matrix-like object
##' @return \code{princompOne}: the first axis of a PCA
##' @rdname traitgram.cc
##' @importFrom stats princomp
##' @export
princompOne <- function(x, ...) princomp(x, ...)$scores[,1]
