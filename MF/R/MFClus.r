#' Estimates mitigated fraction from clustered or stratified data.
#' 
#' Averages the U statistic over the clusters and computes MF from it. Clusters are excluded if they do not include both treatments.
#' 
#' @title Clustered mitigated fraction
# @usage MFClus(formula, data, compare = c("con", "vac"), trace.it = FALSE)
#' @param formula Formula of the form \code{y ~ x + cluster(w)}, where y is a continuous response, x is a factor with two levels of treatment, and w is a factor indicating the clusters.
#' @param data Data frame.  See \code{Note} for handling of input data with more than two levels.
#' @param compare Text vector stating the factor levels - \code{compare[1]} is the control or reference group to which \code{compare[2]} is compared
#' @param trace.it Verbose tracking of the cycles? Default FALSE.
#' @return a \code{\link{mfcluster-class}} data object
#' @note
#' If input data contains more than two levels of treatment, rows associated with unused treatment levels will be removed. \cr
#' Factor levels for treatments not present in the input data will be ignored. \cr
#' Clusters with missing treatments will be excluded. See \code{\link{mfbootcluster-class}} or use \code{trace.it} to identify excluded clusters.
#' @export
#' @references Siev D. (2005). An estimator of intervention effect on disease severity. \emph{Journal of Modern Applied Statistical Methods.} \bold{4:500--508}
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @seealso \code{\link{mfcluster-class}}
#' @examples
#' \dontrun{
#' MFClus(lesion ~ group + cluster(litter), piglung)
#' 
#' #  Comparing vac to con 
#' #  
#' #  MF = 0.3533835 
#' #  
#' #  By Cluster
#' #     w  u         r n1 n2         mf
#' #  U 25 10 0.4000000  5  5 -0.2000000
#' #  K 12  2 0.2500000  4  2 -0.5000000
#' #  Z 16 10 0.8333333  3  4  0.6666667
#' #  D  3  2 1.0000000  1  2  1.0000000
#' #  N  1  0 0.0000000  1  3 -1.0000000
#' #  T  8  5 0.8333333  2  3  0.6666667
#' #  P  4  1 0.5000000  2  1  0.0000000
#' #  L  3  2 0.6666667  1  3  0.3333333
#' #  G 15  9 0.7500000  3  4  0.5000000
#' #  J 15  9 1.0000000  3  3  1.0000000
#' #  W  6  3 0.7500000  2  2  0.5000000
#' #  A  9  3 0.3333333  3  3 -0.3333333
#' #  X 12  6 1.0000000  3  2  1.0000000
#' #  F 13  7 0.7777778  3  3  0.5555556
#' #  S 21 11 0.9166667  4  3  0.8333333
#' #  H 14  8 0.8888889  3  3  0.7777778
#' #  Y  2  1 1.0000000  1  1  1.0000000
#' #  E  2  1 1.0000000  1  1  1.0000000
#' #  
#' #  All
#' #        w  u         r n1 n2        mf
#' #  All 181 90 0.6766917 50 52 0.3533835
#' #  
#' #  Excluded Clusters
#' #  [1] M, Q, R, B, O, V, I, C
#' }
##
##--------------------------------------------------------------------
## Clustered or Stratified MF
##--------------------------------------------------------------------
##
MFClus <- function(formula, data, compare = c("con", "vac"), trace.it = FALSE){
    # formula of the form response ~ treatment + cluster(clustername)
    # based on prob{F(y)<F(x)}
    # within-cluster ranking only 
    # 3/19/01 initial coding
    # revised 10/3/06 to eliminate clusters without both treatments represented
	# revised 8/27/13 - remove group levels if no observations from that level are present in original data 
	# revised 9/03/13 - subset initial data by comparison group levels
	# revised 9/03/13 - move data reshaping shared by MFClusBoot and MFClus to external function 
	dat <- NULL
	group <- NULL
	clusters <- NULL
	strat <- NULL
	reshapeCluster(data = data, formula = formula, compare = compare, envir = environment()) 
	id <- compare
    out <- matrix(NA, length(strat), 6, dimnames = list(strat, c("w", "u", "r", 
		"n1", "n2", "mf")))
    excluded.clusters <- NULL
    for(stratum in strat) {     
        x <- dat[group == id[1] & as.character(clusters) == stratum]
        y <- dat[group == id[2] & as.character(clusters) == stratum]
        n.x <- length(x)
        n.y <- length(y)
        if(n.x > 0 & n.y > 0) {
            N <- n.x + n.y
            x.y <- c(x, y)
            w <- sum(rank(x.y)[1:n.x])
            u <- w - (n.x * (n.x + 1))/2
            r <- u/(n.x * n.y)
            mf <- 2 * r - 1
            out[stratum,  ] <- c(w, u, r, n.x, n.y, mf)
        }
        else {
            if(trace.it) cat('Cluster',stratum,'missing a treatment\n')
            out[stratum,  ] <- c(NA, NA, NA, n.x, n.y, NA)
            excluded.clusters <- c(excluded.clusters, stratum)
        }
    }
    All <- apply(out, 2, sum, na.rm = T)
    R <- All["u"]/sum(out[, "n1"] * out[, "n2"])
    MF <- 2 * R - 1
    All["r"] <- R
    All["mf"] <- MF
    All <- data.frame(t(All))
    dimnames(All)[[1]] <- 'All'
	return(mfcluster$new(All = All, byCluster = out, excludedClusters = excluded.clusters,
		call = match.call(), compare = compare))
}
