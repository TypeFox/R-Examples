#' Run a simple species migration model
#'
#' \command{kissmig} runs a simple, raster-based, stochastic migration model to simulate species migration and 
#' range shifts. It uses a geographic area of origin along with suitability maps to iteratively run a simple 
#' 3x3 cell algorithm. Specifically, it allows for generating accessibility maps for easy integration of limited 
#' migration in species distribution models (\href{http://onlinelibrary.wiley.com/doi/10.1111/ecog.00930/abstract}{Nobis and Normand 2014}). 
#' @usage
#' kissmig(O, S=NULL, it, type='FOC', signed=FALSE, pext=1.0, pcor=0.2, seed=NULL)
#' @param O a single RasterLayer of the geographic origin
#' @param S a Raster* object of suitability, i.e., a RasterLayer, RasterStack, or RasterBrick
#' @param it number of iteration steps
#' @param type type of result: final distribution ('DIS'), iteration step of first occurrence 
#' ('FOC'), iteration step of last occurrence ('LOC'), or number of iteration steps with occurrence ('NOC')
#' @param signed if TRUE, the sign indicates whether the cells was colonized (positive) or uncolonized (negative)
#' after the last iteration step
#' @param pext propability [0,1] a colonized cell becomes uncolonized between iteration steps, i.e., the species gets locally extinct 
#' @param pcor propability [0,1] corner cells are considered in the 3x3 cell neighborhood 
#' @param seed integer used to set the seed of the random number generator
#' @details
#' Starting from origin "O" \command{kissmig} simulates migration for "it" iteration steps in a heterogeneous environment 
#' characterised by the suitability layer(s) "S". The colonized cells of the origin "O" have value 1, uncolonized
#' cells value 0. In case "S" consists of several suitability layers to cover environmental change, "it" is applied to each
#' layer. Suitability ranges between 0 (unsuitable) and 1 (suitability maximum). \command{kissmig} uses a 3x3 algorithm
#' for species spread/migration. All cells get exstinct before an iteration step with probability "pext", and for
#' a recolonization or new colonization event corner cells within the 3x3 neighborhood are considers
#' with probability "pcor" ("pcor"=0.2 produces more realistic circular spread patterns - see Nobis & Normand 2014).
#' For runtime optimization, signed results are generate for "signed"=TRUE, i.e, in addtion to
#' the result type 'FOC, 'LCO', or 'NOC', the sign indicates the final distribution ('DIS') with positive values
#' beeing colonized and negative values beeing previously colonized but uncolonized after the last iteration step.
#' To get reproducible results, the seed of the R random number generator can be set using the "seed" parameter.
#' @examples
#' \donttest{
#' library(kissmig)
#' 
#' s <- kissmigDummyS(mean=12,sd=3, download=TRUE) # a suitability map
#' o <- kissmigOrigin(s, 8, 44.5, 0.5)             # geographic origin
#' l <- s>=0                                       # land mask
#' plot(s+o, asp=1.0, main='suitability + origin')
#' 
#' # run kissmig with different type of output
#' 
#' k <- kissmig(o, s, it=150, type='DIS')
#' plot(k*l, asp=1.0, main='final distribution')
#' 
#' k <- kissmig(o, s, it=150, type='FOC')
#' plot(k*l, asp=1.0, main='first iteration step of occurrence')
#' 
#' a <- kissmigAccess(k)
#' plot(a*l, asp=1.0, main='accessibility based on "FOC", absolute values')
#' 
#' a <- kissmigAccess(k, rel=TRUE)
#' plot(a*l, asp=1.0, main='accessibility based on "FOC", relative values')
#' 
#' k <- kissmig(o, s, it=150, type='LOC')
#' plot(k*l, asp=1.0, main='last iteration step of occurrence')
#' 
#' k <- kissmig(o, s, it=150, type='NOC')
#' plot(k*l, asp=1.0, main='number of iteration steps with occurrences')
#' }
#'
#' @references \itemize{
#' \item{\href{http://onlinelibrary.wiley.com/doi/10.1111/ecog.00930/abstract}{Nobis MP and Normand S (2014)} KISSMig - a simple model for R to account for
#' limited migration in analyses of species distributions. \cite{Ecography} 37: 1282-1287.}
#' \item{KISSMig homepage \url{http://purl.oclc.org/wsl/kissmig}}}
#' @import raster
#' @useDynLib kissmig
#' @seealso \code{\link{kissmigAccess}, \link{kissmigOrigin}}
#' @export kissmig
#'

kissmig <- function(O, S=NULL, it, type='FOC', signed=FALSE, pext=1.0, pcor=0.2, seed=NULL) {

  # check class of origin 'O' and read data
  if (class(O) != 'RasterLayer') stop("origin 'O' must be a RasterLayer")
  ans <- O
  ov  <- values(O)
  
  # check 'type'
  type <- toupper(type)
  ifelse(type %in% c('DIS','FOC','LOC','NOC'),
    ty <- which(c('DIS','FOC','LOC','NOC')==toupper(type)),
    stop("'type' must be 'DIS', 'FOC', 'LOC', or 'NOC'", call. = FALSE)
  )
  
  # check 'S', prepare suitability vector and dimensions
  ifelse(is.null(S), 
    { sv <- rep(1.0, ncell(O))
      dh <- dim(O)
      warning('no suitability data found - globally set to 1.0', call. = FALSE)
    },{
      if (!(class(S) %in% c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
        stop("suitability 'S' must be a RasterLayer, RasterStack, or RasterBrick")
      }  
      compareRaster(O,S)           # stops if not the same
      sv <- as.vector(values(S))   # converts matrix in calse of brick or stack
      dh <- dim(S)
    }) 

  # R RNG (seed isn't a parameter of kissmig_C)
  if (!is.null(seed)) set.seed(seed)
  
  # signed?
  si <- ifelse(signed, 1, 0)
  
  v <- .Call('kissmig_c', 
             as.double(ov), 
             as.double(sv), 
             as.integer(dh), 
             as.integer(it),
             as.double(pext),
             as.double(pcor),
             as.integer(ty),
             as.integer(si),
			 PACKAGE='kissmig'
  )
  
  values(ans) <- v
  return(ans)
  
}
