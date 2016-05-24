
#'Distribution of the number of responses assuming marginal compatibility.
#'
#'The \code{mc.est} function estimates the distribution of the number of
#'responses in a cluster under the assumption of marginal compatibility:
#'information from all cluster sizes is pooled. The estimation is performed
#'independently for each treatment group.
#'
#'The EM algorithm given by Stefanescu and Turnbull (2003) is used for the binary data.
#'
#'@useDynLib CorrBin
#'@export
#'@param object a \code{\link{CBData}} or \code{\link{CMData}} object
#'@param \dots other potential arguments; not currently used
#'@return For \code{CBData}: A data frame giving the estimated pdf for each treatment and
#'clustersize.  The probabilities add up to 1
#'for each \code{Trt}/\code{ClusterSize} combination. It has the following columns: 
#'@return \item{Prob}{numeric, the probability of \code{NResp} responses in a
#'cluster of size \code{ClusterSize} in group \code{Trt}}
#'@return \item{Trt}{factor, the treatment group}
#'@return \item{ClusterSize}{numeric, the cluster size}
#'@return \item{NResp}{numeric, the number of responses}
#'@author Aniko Szabo
#'@references Stefanescu, C. & Turnbull, B. W. (2003) Likelihood inference for
#'exchangeable binary data with varying cluster sizes.  \emph{Biometrics}, 59,
#'18-24
#'@keywords nonparametric models
#'@examples
#'
#'data(shelltox)
#'sh.mc <- mc.est(shelltox)
#'
#'library(lattice)
#'xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, subset=ClusterSize>0, 
#'    type="l", as.table=TRUE, auto.key=list(columns=4, lines=TRUE, points=FALSE),
#'    xlab="Number of responses", ylab="Probability P(R=r|N=n)")
#'
#'@name mc.est

mc.est <- function(object,...) UseMethod("mc.est")
 

#'Test the assumption of marginal compatibility
#'
#'\code{mc.test.chisq} tests whether the assumption of marginal compatibility is
#'violated in the data.
#'
#'The assumption of marginal compatibility (AKA reproducibility or interpretability) implies that
#'the marginal probability of response does not depend on clustersize.
#'Stefanescu and Turnbull (2003), and Pang and Kuk (2007) developed a
#'Cochran-Armitage type test for trend in the marginal probability of success
#'as a function of the clustersize. \code{mc.test.chisq} implements a
#'generalization of that test extending it to multiple treatment groups.
#'
#'@export
#'@param object a \code{\link{CBData}} or \code{\link{CMData}} object
#'@param \dots other potential arguments; not currently used
#'@return A list with the following components:
#'@return \item{overall.chi}{the test statistic; sum of the statistics for each
#'group}
#'@return \item{overall.p}{p-value of the test}
#'@return \item{individual}{a list of the results of the test applied to each
#'group separately: \itemize{ \item chi.sq the test statistic for the group
#'\item p p-value for the group}}
#'@author Aniko Szabo
#'@seealso \code{\link{mc.est}} for estimating the distribution under marginal
#'compatibility.
#'@references Stefanescu, C. & Turnbull, B. W. (2003) Likelihood inference for
#'exchangeable binary data with varying cluster sizes. \emph{Biometrics}, 59,
#'18-24
#'
#'Pang, Z. & Kuk, A. (2007) Test of marginal compatibility and smoothing
#'methods for exchangeable binary data with unequal cluster sizes.
#'\emph{Biometrics}, 63, 218-227
#'@keywords htest
#'@examples
#'
#'data(shelltox)
#'mc.test.chisq(shelltox)
#' 

mc.test.chisq <- function(object,...) UseMethod("mc.test.chisq")
 
