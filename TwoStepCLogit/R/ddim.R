#' Data Dimension Statistics
#' 
#' Function that computes dimension statistics for a data set with clusters and strata
#' and its \code{print} method.
#' 
#' @param formula A formula object, with the response on the left of a \code{~} operator 
#'                and, on the right hand side, a \code{strata} and a \code{cluster} term 
#'                (ex. \code{formula = Y ~ strata(var_strata) + cluster(var_cluster)}).
#'                The \code{strata} and \code{cluster} functions (from the package survival) are
#'                used to identify the stratification and the cluster variables, respectively.
#' @param data A data frame (or object coercible by as.data.frame to a data frame) containing the 
#'             variables in the model.
#' 
#' @return \item{Sc}{ The number of strata in each cluster.}
#' @return \item{Ystat}{ A data.frame with \code{n}, the numbers of observations per stratum 
#'                       (\eqn{n^c_s}{ncs}), and \code{m}, the sum of the responses per stratum 
#'                       (\eqn{m^c_s}{mcs}).}
#' 
#' @seealso \code{\link{Ts.estim}}
#' @export
#' @importFrom stats aggregate
#' @examples
#' dimstat <- ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison)
#' dimstat
ddim <- function(formula, data){
  
  call <- match.call()
  
  # Validation des arguments formula et data
  if (is.null(call$data)) stop("a 'data' argument is required", call. = FALSE)
  if (is.null(call$formula)) stop("a 'formula' argument is required", call. = FALSE)
  data.name <- deparse(call$data)
  out.shared <- shared(formula = formula, data=data, data.name=data.name)
  info.cluster <- mm <- NULL # inutilisees dans ddim(), mais dans la sortie de shared car utile a Ts.estim()
  info.strata <- y <- var.cluster <- NULL
  for(i in 1:length(out.shared)) assign(names(out.shared)[i],out.shared[[i]]) 
  
  # Pour retirer l'information sur les strates
  Ts.strata <- function(x){x}  ## genre de fonction strata seulement pour la fonction Ts.estim car la fonction strata 
  ## de survival n'est pas adequate ici, a l'image de la fonction cluster de survival
  var.strata <- eval(parse(text=paste("Ts",info.strata$vars,sep=".")), envir=data)
  if (length(dim(var.strata)) > 1) stop("in 'formula', the stratum identifier must be a single variable")
  Sc <- tapply(var.strata, var.cluster, function(x){length(unique(x))})
  
  # Pour retirer de l'info sur la dimension des donnees
  ncs <- aggregate(y, list(cluster=var.cluster, strata=var.strata), length)
  mcs <- aggregate(y, list(cluster=var.cluster, strata=var.strata), sum)
  Ystat <- cbind(ncs[, 1:2], "n"=ncs[, 3], "m"=mcs[, 3])
  Ystat <- Ystat[order(Ystat$cluster, Ystat$strata),]
  rownames(Ystat) <- NULL
  
  # Sortie des resultats
  out <- list(Sc=Sc, Ystat=Ystat)
  class(out) <- "ddim"
  return(out)
}


#' @rdname ddim
#' @method print ddim
#' @param x An object, produced by the \code{\link{ddim}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
#' @export
#' @importFrom stats median
"print.ddim" <- function(x, ...) { 
  cat("Number of clusters =",length(x$Sc))
  cat("\n\nNumbers of strata per cluster")
  if (length(unique(x$Sc))==1) {
    cat(" =",x$Sc[1],"\n")  # if the number of strata doesn't vary among clusters
  } else {
    cat(":\n")
    stat <- c(Min=min(x$Sc), Median=median(x$Sc), Mean=mean(x$Sc), 
        Max=max(x$Sc), Sum=sum(x$Sc))
    print.default(round(stat,0), print.gap = 2, quote = FALSE, right=TRUE, ...)
  }
  cat("\nNumbers of observations per strata (n) and\nsums of the response variable values per strata (m) :\n")
  nm <- unique(x$Ystat[,c("n","m")])
  if (nrow(nm)==1) { # if n and m don't vary among strata
    nmf <- as.matrix(x$Ystat[1,c("n","m"),drop=FALSE])
    rownames(nmf) <- " "
    print.default(nmf, print.gap = 2, quote = FALSE, right=TRUE, ...)  
  } else if(nrow(nm) < 20) {
    tab <- table(x$Ystat[,c("n","m")])
    freq <- vector(length=nrow(nm))
    for (i in 1:nrow(nm)) {
      freq[i] <- tab[paste(nm[i,1]),paste(nm[i,2])]  
    }
    nmf <- as.matrix(cbind(nm, "NumberOfStrata"=freq))
    rownames(nmf) <- rep(" ", nrow(nmf))
    print.default(nmf, print.gap = 2, quote = FALSE, right=TRUE, ...)
  } else {
    print(summary(x$Ystat[,c("n","m")], digits=1))
  }
  cat("\nTotal number of observations =",sum(x$Ystat[,"n"]),"\n")
  
  invisible(x)
}


