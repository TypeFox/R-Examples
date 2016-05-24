#' Calculate similarity score for neuron morphologies
#'
#' Uses the NBLAST algorithm that compares the morphology of two neurons. For
#' more control over the parameters of the algorithm, see the arguments of
#' \code{\link{NeuriteBlast}}.
#'
#' @details when \code{smat=NULL} options("nat.nblast.defaultsmat") will be
#'   checked and if NULL, then \code{smat.fcwb} or \code{smat_alpha.fcwb} will
#'   be used depending on the value of \code{UseAlpha}.
#'
#'   When \code{OmitFailures} is not \code{NA}, individual nblast calls will be
#'   wrapped in \code{try} to ensure that failure for any single neuron does not
#'   abort the whole nblast call. When \code{OmitFailures=FALSE}, missing values
#'   will be left as \code{NA}. \code{OmitFailures=TRUE} is not (yet)
#'   implemented. If you want to drop scores for neurons that failed you will
#'   need to set \code{OmitFailures=FALSE} and then use \code{\link{na.omit}} or
#'   similar to post-process the scores.
#'
#'   Note that when \code{OmitFailures=FALSE} error messages will not be printed
#'   because the call is wrapped as \code{try(expr, silent=TRUE)}.
#'
#'   Internally, the \code{\link{plyr}} package is used to provide options for
#'   parallelising NBLASTs and displaying progress. To display a progress bar as
#'   the scores are computed, add \code{.progress="text"} to the arguments
#'   (non-text progress bars are available -- see
#'   \code{\link[plyr]{create_progress_bar}}). To parallelise, add
#'   \code{.parallel=TRUE} to the arguments. In order to make use of parallel
#'   calculation, you must register a parallel backend that will distribute the
#'   computations. There are several possible backends, the simplest of which is
#'   the multicore option made available by \code{doMC}, which spreads the load
#'   across cores of the same machine. Before using this, the backend must be
#'   registered using \code{registerDoMC} (see example below).
#'
#' @param query the query neuron.
#' @param target a \code{\link[nat]{neuronlist}} to compare neuron against.
#'   Defaults to \code{options("nat.default.neuronlist")}. See
#'   \code{\link{nat-package}}.
#' @param smat the scoring matrix to use (see details)
#' @param sd Standard deviation to use in distance dependence of nblast v1
#'   algorithm. Ignored when \code{version=2}.
#' @param version the version of the algorithm to use (the default, 2, is the
#'   latest).
#' @param normalised whether to divide scores by the self-match score of the
#'   query
#' @param UseAlpha whether to consider local directions in the similarity
#'   calculation (default: FALSE).
#' @param OmitFailures Whether to omit neurons for which \code{FUN} gives an
#'   error. The default value (\code{NA}) will result in nblast stopping with an
#'   error message the moment there is an eror. For other values, see details.
#' @param ... Additional arguments passed to NeuriteBlast or the function used
#'   to compute scores from distances/dot products. (expert use only).
#' @return Named list of similarity scores.
#' @section NBLAST Versions: The \code{nblast} version argument presently
#'   exposes two versions of the algorithm; both use the same core procedure of
#'   aligning two vector clouds, segment by segment, and then computing the
#'   distance and absolute dot product between the nearest segment in the target
#'   neuron for every segment in the query neuron. However they differ
#'   significantly in the procedure used to calculate a score using this set of
#'   distances and absolute dot products.
#'
#'   \bold{Version 1} of the algorithm uses a standard deviation (argument
#'   \code{sd}) as a user-supplied parameter for a negative exponential
#'   weighting function that determines the relationship between score and the
#'   distance between segments. This corresponds to the parameter \eqn{\sigma}
#'   in the weighting function:
#'
#'   \eqn{f=\sqrt{|\vec{u_{i}}\cdot\vec{v_{i}}|\exp\left(-d_{i}^{2}/2\sigma^{2}\right)}}
#'
#'
#'   This is the same approach described in Kohl et al 2013 and the similarity
#'   scores in the interval (0,1) described in that paper can exactly
#'   recapitulated by setting \code{version=1} and \code{normalised=TRUE}.
#'
#'   \bold{Version 2} of the algorithm is described in Costa et al 2014. This
#'   uses a more sophisticated and principled scoring approach based on a
#'   log-odds ratio defined by the distribution of matches and non-matches in
#'   sample data. This information is passed to the nblast function in the form
#'   of a \emph{scoring matrix} (which can be computed by
#'   \code{\link{create_scoringmatrix}}); a default scoring matrix
#'   \code{\link{smat.fcwb}} has been constructed for \emph{Drosophila} neurons.
#'
#'   \bold{Which version should I use?} You should use version 2 if you are
#'   working with \emph{Drosophila} neurons or you have sufficient training data
#'   (in the form of validated matching and random neuron pairs to construct a
#'   scoring matrix). If this is not the case, you can always fall back to
#'   version 1, setting the free parameter (sd or \eqn{\sigma}) to a value that
#'   encapsulates your understanding of the location precision of neurons in
#'   your species/brain region of interest. In the fly brain we have used
#'   \eqn{\sigma=3} microns, since previous estimates of the localisation of
#'   identifiable features of neurons (Jefferis, Potter et al 2007) are of this
#'   order.
#'
#' @references Kohl, J. Ostrovsky, A.D., Frechter, S., and Jefferis, G.S.X.E
#' (2013). A bidirectional circuit switch reroutes pheromone signals in male and
#' female brains. Cell 155 (7), 1610--23
#' \href{http://dx.doi.org/10.1016/j.cell.2013.11.025}{doi:
#' 10.1016/j.cell.2013.11.025}.
#'
#' Costa, M., Ostrovsky, A.D., Manton, J.D., Prohaska, S., and Jefferis,
#' G.S.X.E. (2014). NBLAST: Rapid, sensitive comparison of neuronal structure
#' and construction of neuron family databases. Biorxiv preprint.
#' \href{http://dx.doi.org/10.1101/006346}{doi: 10.1101/006346}.
#'
#' Jefferis G.S.X.E., Potter C.J., Chan A.M., Marin E.C., Rohlfing T., Maurer
#' C.R.J., and Luo L. (2007). Comprehensive maps of Drosophila higher olfactory
#' centers: spatially segregated fruit and pheromone representation. Cell 128
#' (6), 1187--1203.
#' \href{http://dx.doi.org/10.1016/j.cell.2007.01.040}{doi:10.1016/j.cell.2007.01.040}
#'
#' @seealso \code{\link{nat-package}}, \code{\link{nblast_allbyall}},
#'   \code{\link{create_scoringmatrix}}, \code{\link{smat.fcwb}}
#' @export
#' @importFrom nat neuronlist
#' @examples
#' # load sample Kenyon cell data from nat package
#' data(kcs20, package='nat')
#' # search one neuron against all neurons
#' scores=nblast(kcs20[['GadMARCM-F000142_seg002']], kcs20)
#' # scores from best to worst, top hit is of course same neuron
#' sort(scores, decreasing = TRUE)
#' hist(scores, breaks=25, col='grey')
#' abline(v=1500, col='red')
#'
#' # plot query neuron
#' open3d()
#' # plot top 3 hits (including self match with thicker lines)
#' plot3d(kcs20[which(sort(scores, decreasing = TRUE)>1500)], lwd=c(3,1,1))
#' rest=names(which(scores<1500))
#' plot3d(rest, db=kcs20, col='grey', lwd=0.5)
#'
#' # normalised scores (i.e. self match = 1) of all neurons vs each other
#' # note use of progress bar
#' scores.norm=nblast(kcs20, kcs20, normalised = TRUE, .progress="text")
#' hist(scores.norm, breaks=25, col='grey')
#' # produce a heatmap from normalised scores
#' jet.colors <- colorRampPalette( c("blue", "green", "yellow", "red") )
#' heatmap(scores.norm, labCol = with(kcs20,type), col=jet.colors(20), symm = TRUE)
#'
#' \dontrun{
#' # Parallelise NBLASTing across 4 cores using doMC package
#' library(doMC)
#' registerDoMC(4)
#' scores.norm2=nblast(kcs20, kcs20, normalised=TRUE, .parallel=TRUE)
#' stopifnot(all.equal(scores.norm2, scores.norm))
#' }
nblast <- function(query, target=getOption("nat.default.neuronlist"),
                   smat=NULL, sd=3, version=c(2, 1), normalised=FALSE,
                   UseAlpha=FALSE, OmitFailures=NA, ...) {
  version <- as.character(version)
  version <- match.arg(version, c('2', '1'))

  if(is.character(target)) {
    target=get(target, envir = parent.frame())
  }
  if(is.null(target)) stop("No target neurons provided. Please set them directly",
                           " or use option('nat.default.neuronlist')")

  # Convert target to neuronlist if passed a single object
  if("dotprops" %in% class(target)) target <- neuronlist(target)

  if(version == '2') {
    if(is.null(smat)) {
      smat=getOption("nat.nblast.defaultsmat")
      if(is.null(smat)) {
        if(UseAlpha) smat="smat_alpha.fcwb"
        else smat="smat.fcwb"
      }
    }
    if(is.character(smat)) smat=get(smat)
    NeuriteBlast(query=query, target=target, NNDistFun=lodsby2dhist, smat=smat,
                 UseAlpha=UseAlpha, normalised=normalised, OmitFailures=OmitFailures, ...)
  } else if(version == '1') {
    NeuriteBlast(query=query, target=target, NNDistFun=WeightedNNBasedLinesetDistFun,
                 UseAlpha=UseAlpha, sd=sd, normalised=normalised,
                 OmitFailures=OmitFailures, ...)
  } else {
    stop("Only NBLAST versions 1 and 2 are currently implemented. For more advanced control, see NeuriteBlast.")
  }
}

#' Wrapper function to compute all by all NBLAST scores for a set of neurons
#'
#' @description Calls \code{nblast} to compute the actual scores. Can accept
#'   either a neuronlist or neuron names as a character vector. This is a thin
#'   wrapper around nblast and its main advantage is the option of "mean"
#'   normalisation for forward and reverse scores, which is the most sensible
#'   input to give to a clustering algorithm as well as the choice of returning
#'   a distance matrix.
#'
#' @details Note that \code{nat} already provides a function
#'   \code{\link{nhclust}} for clustering, which is a wrapper for R's
#'   \code{hclust} function. \code{nhclust} actually expects \bold{raw} scores
#'   as input.
#'
#' @section TODO: It would be a good idea in the future to implement a parallel
#'   version of this function.
#' @param x Input neurons (neuronlist or character vector)
#' @param smat the scoring matrix to use (see details of \code{\link{nblast}}
#'   for meaning of default \code{NULL} value)
#' @param ... Additional arguments for methods or \code{nblast}
#' @export
#' @seealso \code{\link{nblast}, \link{sub_score_mat}, \link{nhclust}}
#' @examples
#' library(nat)
#' kcs20.scoremat=nblast_allbyall(kcs20)
#' kcs20.hclust=nhclust(scoremat=kcs20.scoremat)
#' plot(kcs20.hclust)
nblast_allbyall<-function(x, ...) UseMethod("nblast_allbyall")

#' @rdname nblast_allbyall
#' @param db A neuronlist or a character vector naming one. Defaults to value of
#'   options("nat.default.neuronlist")
#' @export
nblast_allbyall.character<-function(x, smat=NULL, db=getOption("nat.default.neuronlist"), ...){
  if(is.character(db)) {
    db=get(db, envir = parent.frame())
  }
  nblast_allbyall(db[x], smat=smat, ...)
}

#' @rdname nblast_allbyall
#' @inheritParams sub_score_mat
#' @export
nblast_allbyall.neuronlist<-function(x, smat=NULL, distance=FALSE,
                                     normalisation=c("raw","normalised","mean"),
                                     ...){
  normalisation=match.arg(normalisation)
  scores=nblast(x, x, smat=smat, normalised=FALSE, ...)
  if(length(x)==1) {
    # turn this into a matrix
    scores=matrix(scores, dimnames=list(names(x), names(x)))
  }
  if(normalisation !='raw')
    sub_score_mat(scoremat=scores, distance = distance, normalisation=normalisation)
  else scores
}

#' Produce similarity score for neuron morphologies
#'
#' A low-level entry point to the NBLAST algorithm that compares the morphology
#' of a neuron with those of a list of other neurons. For most use cases, one
#' would probably wish to use \code{\link{nblast}} instead.
#'
#' @details For detailed description of the \code{OmitFailures} argument, see
#'   the details section of \code{\link{nblast}}.
#' @param query either a single query neuron or a \code{\link[nat]{neuronlist}}
#' @param target a \code{neuronlist} to compare neuron against.
#' @param targetBinds numeric indices or names with which to subset
#'   \code{target}.
#' @param normalised whether to divide scores by the self-match score of the
#'   query
#' @param simplify whether to simplify the scores from a list to a vector.
#'   \code{TRUE} by default. The only time you might want to set this false is
#'   if you are collecting something other than simple scores from the search
#'   function. See \code{\link{simplify2array}} for further details.
#' @param ... extra arguments to pass to the distance function.
#' @inheritParams nblast
#' @return Named list of similarity scores.
#' @importFrom nat is.neuronlist
#' @export
#' @seealso \code{\link{WeightedNNBasedLinesetMatching}}
NeuriteBlast <- function(query, target, targetBinds=NULL, normalised=FALSE,
                         OmitFailures=NA, simplify=TRUE, ...){
  if(isTRUE(OmitFailures))
    stop("OmitFailures=TRUE is not yet implemented!\n",
         "Use OmitFailures=FALSE and handle NAs to taste.")

  if("neuron" %in% class(target)) target <- neuronlist(target)

  if(nat::is.neuronlist(query)) {
    res=plyr::llply(query, NeuriteBlast, target=target, targetBinds=targetBinds,
                  normalised=normalised, OmitFailures=OmitFailures, simplify=simplify, ...=...)
    if (!identical(simplify, FALSE) && length(res))
      res=simplify2array(res, higher = (simplify == "array"))
    return(res)
  } else {
    if(is.null(targetBinds))
      targetBinds=seq_along(target)
    else if(is.character(targetBinds))
      targetBinds=charmatch(targetBinds,names(target))
    scores=rep(NA_real_, length(targetBinds))
    FUN = if(is.na(OmitFailures)) {
      # will error out if something goes wrong
      WeightedNNBasedLinesetMatching
    } else {
      function(...) {
        score=try(WeightedNNBasedLinesetMatching(...), silent = TRUE)
        if(inherits(score,'try-error')) NA_real_ else score
      }
    }
    scores=plyr::llply(targetBinds, function(i, ...) FUN(target[[targetBinds[i]]], query, ...), ... )
    names(scores)=names(target)[targetBinds]
    # unlist if these are simple numbers
    if (!identical(simplify, FALSE) && length(scores))
      scores=simplify2array(scores, higher = (simplify == "array"))
  }

  if(normalised){
    if(is.list(scores)) stop("Cannot normalise results when they are not a single number")
    scores=scores/NeuriteBlast(query, neuronlist(query), normalised=FALSE, ...)
  }
  scores
}


#' @importFrom stats dnorm
WeightedNNBasedLinesetDistFun<-function(nndists,dotproducts,sd,...){
  summaryfun=function(x) 1-mean(sqrt(x),na.rm=T)
  sapply(sd,function(sd) summaryfun(dnorm(nndists,sd=sd)*dotproducts/dnorm(0,sd=sd)))
}


#' Compute point & tangent vector similarity score between two linesets
#'
#' @description WeightedNNBasedLinesetMatching will work with 2 objects of class
#'   \code{dotprops} or \code{neuron}. The code to calculate scores directly for
#'   \code{neuron} objects gives broadly comparable scores to that for
#'   \code{dotprops} objects, but has been lightly tested. Furthermore only
#'   objects in \code{dotprops} form were used in the constriction of the
#'   scoring matrices distributed in this package. It is therefore recommended
#'   to convert neuron objects to dotprops objects using the
#'   \code{\link{dotprops}} function.
#' @rdname WeightedNNBasedLinesetMatching
#' @param target,query dotprops or neuron objects to compare (must be of the
#'   same class)
#' @return Value of NNDistFun passd to WeightedNNBasedLinesetMatching
#' @importFrom nabor knn
#' @export
WeightedNNBasedLinesetMatching <- function(target, query, ...) {
  UseMethod("WeightedNNBasedLinesetMatching")
}


#' @details \code{UseAlpha} determines whether the alpha values
#'   (eig1-eig2)/sum(eig1:3) are passed on to WeightedNNBasedLinesetMatching.
#'   These will be used to scale the dot products of the direction vectors for
#'   nearest neighbour pairs.
#' @param UseAlpha Whether to scale dot product of tangent vectors (default=F)
#' @param ... extra arguments to pass to the distance function.
#' @export
#' @seealso \code{\link[nat]{dotprops}}
#' @rdname WeightedNNBasedLinesetMatching
#' @importFrom nat dotprops
WeightedNNBasedLinesetMatching.dotprops<-function(target, query, UseAlpha=FALSE, ...) {
  if(!"dotprops" %in% class(query)) query <- dotprops(query)
  if(UseAlpha)
    WeightedNNBasedLinesetMatching(target$points,query$points,dvs1=target$vect,dvs2=query$vect,
                                   alphas1=target$alpha,alphas2=query$alpha,...)
  else
    WeightedNNBasedLinesetMatching(target$points,query$points,dvs1=target$vect,dvs2=query$vect,...)
}


#' @export
#' @param OnlyClosestPoints Whether to restrict searches to the closest points
#'   in the target (default FALSE, only implemented for dotprops).
#' @rdname WeightedNNBasedLinesetMatching
#' @importFrom nat dotprops
WeightedNNBasedLinesetMatching.neuron<-function(target, query, UseAlpha=FALSE,
                                                OnlyClosestPoints=FALSE,...) {
  if(!"neuron" %in% class(query)) {
    target <- dotprops(target)
    return(WeightedNNBasedLinesetMatching(target=target, query=query, UseAlpha=UseAlpha, OnlyClosestPoints=OnlyClosestPoints, ...))
  }
  if(UseAlpha)
    stop("UseAlpha is not yet implemented for neurons!")
  if(OnlyClosestPoints)
    stop("OnlyClosestPoints is not yet implemented for neurons!")
  target=data.matrix(target$d[,c("X","Y","Z","Parent")])
  query=data.matrix(query$d[,c("X","Y","Z","Parent")])
  WeightedNNBasedLinesetMatching(target, query, ...)
}


WeightedNNBasedLinesetMatching.default<-function(target,query,dvs1=NULL,dvs2=NULL,alphas1=NULL,
                                        alphas2=NULL,NNDistFun=WeightedNNBasedLinesetDistFun,Verbose=FALSE,
                                        BothDirections=FALSE,BothDirectionsFun=list,OnlyClosestPoints=FALSE,...){
  # return a score based on the similarity of nearest neighbour location
  # and the dot product of the direction vectors

  NNDistFun=match.fun(NNDistFun)
  if(BothDirections){
    BothDirectionsFun=match.fun(BothDirectionsFun)
    f=WeightedNNBasedLinesetMatching(target,query,dvs1,dvs2,alphas1,alphas2,
                                     NNDistFun=NNDistFun,Verbose=Verbose,BothDirections=FALSE,
                                     OnlyClosestPoints=OnlyClosestPoints,...)
    b=WeightedNNBasedLinesetMatching(query,target,dvs1,dvs2,alphas1,alphas2,
                                     NNDistFun=NNDistFun,Verbose=Verbose,BothDirections=FALSE,
                                     OnlyClosestPoints=OnlyClosestPoints,...)
    if(length(f)==1 && length(b)==1) return (BothDirectionsFun(f,b))
    if(length(dim(f))==1 && length(f)==length(b)) return (cbind(f,b))
    return(BothDirectionsFun(f,b))
  }

  a=target[,c("X","Y","Z")]
  b=query[,c("X","Y","Z")]

  nntarget=nabor::knn(a,b,k=1)

  idxArray=cbind(nntarget$nn.idx,seq(length(nntarget$nn.idx)))

  # Need to supply a set of pairs of points.
  # will use the parent of each chosen point.
  # if parent undefined, then ignore that point

  if(is.null(dvs1) || is.null(dvs2)){
    if(OnlyClosestPoints==TRUE)
      stop("OnlyClosestPoints is not yet implemented for neurons")
    # Calculate the direction vectors
    dvs=findDirectionVectorsFromParents(target,query,idxArray,ReturnAllIndices=TRUE,Verbose=Verbose)

    # Calculate segment lengths
    l1.seglengths=normbyrow(dvs[,1:3])
    l2.seglengths=normbyrow(dvs[,4:6])
    # normalise the direction vectors
    dvs[,1:3]=dvs[,1:3]/l1.seglengths
    dvs[,4:6]=dvs[,4:6]/l2.seglengths
    # Calculate absolute dot products
    # nb absolute, because we don't really care about directionality here
    dps=abs(dotprod(dvs[,1:3],dvs[,4:6]))
  } else {
    # OnlyClosestPoints prunes the list of query-target pairs so that no
    # points in the target are duplicated (points in query are already unique)
    if(OnlyClosestPoints){
      # sort by increasing distance between pairs
      # remove duplicates in target
      targetdupes=duplicated(nntarget$nn.idx[order(nntarget$nn.dist)])
      idxArray=idxArray[!targetdupes,,drop=FALSE]
      nntarget$nn.dists=nntarget$nn.dists[!targetdupes]
    }
    dps=abs(dotprod(dvs1[idxArray[,1],],dvs2[idxArray[,2],]))
    if(!is.null(alphas1) && !is.null(alphas2)){
      # for perfectly aligned points, alpha = 1, at worst alpha = 0
      # sqrt seems reasonable since if alpha1=alpha2=0.5 then
      # the scalefac will be 0.5
      # zapsmall to make sure there are no tiny negative numbers etc
      scalefac=sqrt(zapsmall(alphas1[idxArray[,1]]*alphas2[idxArray[,2]]))
      dps=dps*scalefac
    }
  }

  NNDistFun(as.vector(nntarget$nn.dists),dps,...)
}


dotprod=function(a,b){
  # expects 2 matrices with n cols each
  c=a*b
  if(length(dim(c))>1) 	rowSums(c)
  else sum(c)
}

# Originally from AnalysisSuite PotentialSynapases.R
normbyrow=function(a){
  # returns euclidean norm (by row if reqd)
  c=a*a
  if(length(dim(c))>1) 	sqrt(rowSums(c))
  else sqrt(sum(c))
}

findDirectionVectorsFromParents<-function(d1,d2,idxArray,ReturnAllIndices=FALSE,Verbose=FALSE){
  # rather than dropping root points, just use the vector from them rather than to them
  if(Verbose) cat(".")
  p1=.CleanupParentArray(d1[,"Parent"])
  p2=.CleanupParentArray(d2[,"Parent"])
  parentPointsArray=cbind(p1[idxArray[,1]],p2[idxArray[,2]])
  # find any points with bad parents and instead use their offspring
  if(any(parentPointsArray[,1]<1 | parentPointsArray[,2]<1)){
    stop ("Some points do not have a parent: therefore impossible to calculate direction vector")
  }
  dvs=cbind(d1[idxArray[,1],c("X","Y","Z"),drop=FALSE]-d1[parentPointsArray[,1],c("X","Y","Z"),drop=FALSE],
            d2[idxArray[,2],c("X","Y","Z"),drop=FALSE]-d2[parentPointsArray[,2],c("X","Y","Z"),drop=FALSE])

  if(ReturnAllIndices){
    attr(dvs,"idxArray")=idxArray
    attr(dvs,"parentPointsArray")=parentPointsArray
  }
  dvs
}

.CleanupParentArray<-function(pa){
  # takes a list of parents for points and replaces any <1
  # with the first offspring of that point
  #if(length(dim(pa))>1) apply(pa,2,.CleanupParentArray)
  pointsNeedingWork<-which(pa<1)
  if(length(pointsNeedingWork)<1) return( pa )
  for(p in pointsNeedingWork){
    wp=which(pa==p)
    if(length(wp)>1){
      warning(cat("more than 1 point in .CleanupParentArray, choosing first from:",wp))
      pa[p]=wp[1]
    } else if(length(wp)<1){
      warning("no points to choose in .CleanupParentArray using original value")
    } else pa[p]=wp
  }
  pa
}


lodsby2dhist <- function(nndists, dotprods, smat, Return=c('sum', 'weightedlodtable', 'elements')) {
  Return <- match.arg(Return)

  if(missing(dotprods)) {
    if(!is.list(nndists))
      stop("must provide nndists and dotprods or a list with both")
    dotprods <- nndists[[2]]
    nndists <- nndists[[1]]
  }
  if(length(nndists)!= length(dotprods))
    stop("nndists and dotprods must have the same length.")

  c1 <- findInterval(nndists, attr(smat,"distbreaks"), all.inside=TRUE)
  # NB use of all.inside fixes NAs that would otherwise result
  # when dot product falls outside (0,1) due to floating point errors
  c2 <- findInterval(dotprods, attr(smat,"dotprodbreaks"), all.inside=TRUE)

  if(Return == 'elements') return(smat[cbind(c1,c2)])

  nlc1 <- length(attr(smat, "distbreaks")) - 1
  nlc2 <- length(attr(smat, "dotprodbreaks")) - 1

  this.countstable <- fast2dintegertable(c1, c2, nlc1, nlc2)

  weightedlodtable <- smat*this.countstable
  if(Return == 'sum') return(sum(weightedlodtable))
  weightedlodtable
}


fast2dintegertable <- function(a, b, nlevelsa=max(a), nlevelsb=max(b)) {
  # Fast 2D cross-tabulation (joint histogram) for two integer inputs

  nlevelsab <- nlevelsa*nlevelsb
  if(nlevelsab > .Machine$integer.max) stop("Number of levels exceeds integer type.")
  ab <- (a - 1) * nlevelsb + b
  matrix(tabulate(ab, nlevelsab), ncol=nlevelsb, nrow=nlevelsa, byrow=T)
}
