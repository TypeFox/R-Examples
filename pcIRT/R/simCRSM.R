#' simulate data according to CRSM
#'
#' With this function data sets according to the Continous Rating Scale Model
#' are simulated
#'
#' The midpoint and the length of the response scale define the interval of the
#' data set generated. The default of the function generates data according to
#' a response scale between 0 and 1 - that is midpoint 0.5 and length 1.
#'
#' @param itempar a numerical vector with item parameters
#' @param disp a number setting the dispersion parameter for the item set
#' @param perspar a numerical vector with the person parameters
#' @param mid the midpoint of the response scale (on which the data set is
#' generated)
#' @param len the length of the response scale (on which the data set is
#' generated)
#' @param seed a seed for the random number generated can optionally be set
#' @return \item{datmat}{simulated data set} \item{true_itempar}{the fixed item
#' parameters according to the input} \item{true_disppar}{the fixed
#' dispersion parameter according to the input} \item{true_perspar}{the fixed
#' person parameters according to the input}
#' @author Christine Hohensinn
#' @seealso \code{\link{simMPRM}}
#' @references Mueller, H. (1987). A Rasch model for continuous ratings.
#' Psychometrika, 52, 165-181.
#' @keywords continuous rating scale model simulation
#' @examples
#'
#' #set item parameters
#' item_p <- c(-1.5,-0.5,0.5,1)
#'
#' #set dispersion parameter for items
#' dis_p <- 5
#'
#' #generate person parameters by a standard normal dispersion
#' pp <- rnorm(50, 0,1)
#'
#' #simulate data set
#' #this is only an illustrating example for simulating data!
#' #In practice, a sample size of n=50 will be too small for most application
#' #demands
#' simdatC <- simCRSM(item_p, dis_p, pp)
#'
#'
#' @export simCRSM
simCRSM <-
function(itempar, disp, perspar, mid=0.5, len=1, seed=NULL){

mat.intercept <- outer(perspar, itempar, "-")

forint <- as.vector(mat.intercept)


funk.n <- function(st, int, mid, disp){exp(st*int+ st*(2*mid-st)*disp)}

nenn  <- sapply(forint, function(i) {
  integrate(funk.n, int=i, mid=mid, disp=disp, lower=(mid-0.5*len), upper=(mid+0.5*len))$value
})

interv <- seq(mid-0.5*len,mid+0.5*len,length.out=500)
zahl.int <- lapply(forint, function(o) {
    sapply(1:(length(interv)-1), function(j) {
      integrate(funk.n, int=o, mid=mid, disp=disp, lower=interv[j], upper=interv[j+1])$value
  })
})

pval <- lapply(1:length(forint), function(y) zahl.int[[y]]/nenn[y])

if (!is.null(seed)) {set.seed(seed)}

p.response.sampled <- sapply(1:length(forint), function(s) sample(1:(length(interv)-1), size=1, prob=pval[[s]]))
p.response <- interv[p.response.sampled]

datmat <- matrix(p.response, ncol=length(itempar))
list(datmat=datmat, true_itempar=itempar, true_disppar=disp, true_perspar=perspar)
}
