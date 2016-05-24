#' Calculate any metric(s) (and compare with null distributions)
#'
#' Allow the calculation of any metric within \code{pez}.
#'
#' \code{generic.null} Calculate metrics and compare with null
#' distributions. Very light wrapper around the utility functions
#' \code{\link{generic.null}} and \code{\link{generic.metrics}} (which
#' is, itself, a very simple function!).
#'
#' @note \code{comp.fun} can be \emph{anything}; much ink has been
#' written about the use of standard effect sizes in eco-phylogenetic
#' analyses (\emph{e.g.}, Kembel 2009). That this function makes it
#' easy for you to compute Standard Effect Sizes does not necessarily
#' mean that you should (see Pearse et al. 2013).
#'
#' Calculating null permutations on a dispersion metric makes little
#' sense, since (by definition; see Pearse et al. 2014) a dispersion
#' metric \emph{require} the use of a null distribution to be
#' calculated. There is nothing to stop you doing so, however! The
#' code makes no attempt to stop you calculating null dissimilarity
#' metrics, but I am not certain that doing so is a good idea using
#' this code as I don't know what to do with such null models!
#'
#' The \code{\link{pez.shape}}, \code{\link{pez.evenness}},
#' \code{\link{pez.dispersion}}, and \code{\link{pez.dissimilarity}}
#' wrapper functions go to some trouble to stop you calculating
#' metrics using inappropriate data (see their notes). These functions
#' give you access to the underlying code within \code{pez}; there is
#' nothing I can do to stop you calculating a metric that, in my
#' opinion, doesn't make any sense. You have been warned :D
#' @param data data \code{\link{comparative.comm}} object
#' @param permute number of null permutations to perform (default
#' 1000)
#' @param metrics vector of functions to be calculated on \code{data};
#' see \link{pez.metrics} for a list of them.
#' @param null.model one of "taxa.labels", "richness", "frequency",
#' "sample.pool", "phylogeny.pool", "independentswap", or
#' "independentswap". These correspond to the null models available in
#' \code{\link{picante}}
#' @param comp.fun comparison function to compare observed values with
#' null values. Default is \code{\link{.ses}}; this is a Standard
#' Effect Size (obs - mean)/SEmean. You may supply your own function;
#' it should take the observed site-metric matrix as its first
#' argument, and a site-metric-permutation array as its second. See
#' the internals of \code{\link{generic.null}} for an example of its
#' use.
#' @param trait if using \code{trait.asm} \code{null.model}, which
#' trait to use (as in \code{\link{trait.asm}})
#' @return \code{generic.null} Output from \code{comp.fun}, by default
#' an array (site-metric-type), where type is the observed value, the
#' mean of the null permutations, the Standard Error of that mean, the
#' Standard Effect Size of the metric (obs-null.mean)/SE, and then the
#' rank of the observed value in the permutations. The rank can be
#' considered a bootstrapped p-value of significance, but remember
#' that this is a rank: at the 95% (two-tailed) level, a value of 0.99
#' would be significant.
#' @param ... additional arguments (e.g, \code{dist},
#' \code{abundance.weighted}) to be passed to any metric functions
#' (see \code{\link{generic.metrics}} for possible arguments)
#' @author Will Pearse
#' @references Kembel S.W. (2009) Disentangling niche and neutral
#' influences on community assembly: assessing the performance of
#' community phylogenetic structure tests. Ecology letters, 12(9),
#' 949-960.
#' @references Pearse W.D., Jones F.A. & Purvis A. (2013) Barro
#' Colorado Island's phylogenetic assemblage structure across fine
#' spatial scales and among clades of different ages. Ecology, 94(12),
#' 2861-2872.
#' @export
#' @rdname generic.metrics
#' @name generic.metrics
#' @examples
#' #Setup data
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' #Calculate some metrics
#' generic.metrics(data, c(.mpd, .pse))
#' #Compare with a trait-based null model (trait.asm)
#' generic.null(data, c(.mpd, .pse), "trait.asm", permute=10, trait="fish.pref")
#' #...be patient when running large (e.g., 1000) sets of null simulations
#' #You can also do this in pieces, giving even more flexibility
#' observed <- generic.metrics(data, c(.mpd, .pse))
#' #null <- .metric.null(data, c(.mpd, .pse))
#' #ses <- .ses(observed, null)
#' #...this is how everything works within generic.null
#' #...and, as with everything in pez, all internal functions start with a "."
generic.null <- function(data, metrics, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap", "trait.asm"), permute=1000, comp.fun=.ses, ...){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(permute < 0) stop("Can't have negative null permutations!")
    null.model <- match.arg(null.model)
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")
    
    #Calculate real values
    observed <- generic.metrics(data, metrics, ...)
    
    #Calculate null distributions and structure
    null <- .metric.null(data, metrics, null.model, permute, ...)

    #Perform comparison and return
    return(comp.fun(observed, null))
}

#' \code{.ses} Calculate Standard Effect Sizes of metrics
#' @param observed observed metric values in site-metric matrix
#' (\emph{e.g.}, from \code{\link{generic.metrics}})
#' @param null null distributions (\emph{e.g.}, from
#' \code{\link{.metric.null}}) in a site-metric-permutation array
#' @return \code{.ses} Vector of standard effect sizes
#' @export
#' @rdname generic.metrics
#' @name generic.metrics
.ses <- function(observed, null){
    means <- apply(null, 1:2, mean)
    ses <- apply(null, 1:2, function(x) sd(x)/sqrt(length(x)))
    rank <- matrix(nrow=nrow(observed), ncol=ncol(observed))
    for(i in seq_len(nrow(observed)))
        for(j in seq_len(ncol(observed)))
            rank[i,j] <- rank(c(observed[i,j], null[i,j,]), ties.method="average")[1]
    output <- array(dim=c(dim(observed),5), dimnames=list(NULL, NULL, c("observed", "null.mean", "SE", "SES", "rank")))
    output[,,1] <- observed
    output[,,2] <- means; output[,,3] <- ses
    output[,,4] <- (observed - means)/ses; output[,,5] <- rank/dim(null)[3]
    return(output)
}

#' \code{.metric.null} Produce null randomisations and compute metrics
#' across them
#' @return \code{.metric.null} site-metric-permutation array
#' @export
#' @rdname generic.metrics
#' @name generic.metrics
.metric.null <- function(data, metrics, null.model=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap", "trait.asm"), permute=1000, trait=-1, ...){
    #Assertions and argument handling
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(permute < 0) stop("Can't have negative null permutations!")
    null.model <- match.arg(null.model)
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")
    if(null.model == "trait.asm"){
        if(is.null(data$data))
            stop("Need traits to calculate trait-based null models")
        trait <- data$data[,trait]
        if(length(trait) != ncol(data$comm))
            stop("'trait' is not a column in trait data")
    }
    
    #Calculate null distributions and structure
    if(null.model=="trait.asm"){
        null <- trait.asm(trait, permute)[[1]]
        rownames(null) <- seq_len(permute)
        colnames(null) <- colnames(data$comm)
        null <- comparative.comm(data$phy, null, data$traits, data$env)
        null <- array(rep(sapply(metrics, function(x) x(null)), permute), c(nrow(data$comm), length(metrics), permute))
    } else {
        null <- apply(replicate(permute, .eco.null(data$comm, method=null.model)), 3, comparative.comm, phy=data$phy)
        null <- lapply(null, function(x) do.call(cbind, lapply(metrics, function(y) y(x))))
        null <- array(unlist(null), c(dim(null[[1]]), length(null)))
    }

    return(null)
}

#' \code{generic.metrics} Calculate arbitrary metrics across data
#' @return \code{generic.metrics} site-metric matrix
#' @export
#' @rdname generic.metrics
#' @name generic.metrics
generic.metrics <- function(data, metrics, ...){
    if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative.comm object")
    if(!is.list(metrics) | !all(sapply(metrics, is.function)))
        stop("'metrics' must be a list of functions")
    
    return(do.call(cbind, lapply(metrics, function(x) x(data, ...))))
}
