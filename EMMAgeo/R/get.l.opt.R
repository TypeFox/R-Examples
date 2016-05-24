#' Identify optimum weight transformation value
#' 
#' This function returns for a series of input vaules that weight 
#' transformation value, which yielded the highest measure of model quality.#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param l Numeric vector of weight transformation values to test.
#' @param quality Charcter scalar, qualtiy measure for against which to test
#' the influence of \code{l}. One out of \code{"mRm"}, \code{"mRn"}, 
#' \code{"mRt"}, \code{"mEm"}, \code{"mEn"} and \code{"mEt"}. Default is 
#' \code{"mRt"}.
#' @param Vqn Numeric matrix specifying optional unscaled user-defined
#' end-member loadings.
#' @param rotation Character scalar, rotation type, default is "Varimax" (cf.
#' Dietze et al., 2012). One out of the rotations provided in GPArotation is
#' possible (cf. \code{\link{rotations}}).
#' @param plot Logical scalar, optional graphical output of the result.
#' @param \dots Further arguments passed to the function.
#' @return Numeric scalar, weight tranformation value with optimal EMMA result.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## TO BE DONE!!!
#' 
#' @export get.l.opt
get.l.opt <- function(
  X, 
  l,
  quality = "mRt",
  Vqn,
  rotation,
  plot = TRUE,
  ...
){
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## test input data
  if(missing(l) == TRUE) {
    l <- 0
  }
  
  if(missing(Vqn) == TRUE) {
    stop("Vqn missing!")
  }
  
  if(missing(rotation) == TRUE) {
    rotation <- "Varimax"
  }
  
  ## create result vector
  quality.out <- matrix(nrow = length(l),
                        ncol = 6)
  
  ## test quality for all l-values
  for(i in 1:length(l)) {
    
    ## run EMMA
    E <- EMMAgeo::EMMA(X = X, 
                       q = nrow(Vqn),
                       l = l[i], 
                       Vqn = Vqn, 
                       rotation = rotation)
    
    ## assign quality measure
    quality.out[i, c(1, 2, 4, 5)] <- c(mean(E$Rm, na.rm = TRUE),
                                       mean(E$Rn, na.rm = TRUE),
                                       mean(E$Em, na.rm = TRUE),
                                       mean(E$En, na.rm = TRUE))
  }
  
  ## complete quality matrix
  quality.out[,3] <- apply(X = quality.out[,1:2], 
                           MARGIN = 1, 
                           FUN = mean,
                           na.rm = TRUE)
  
  quality.out[,6] <- apply(X = quality.out[,4:5], 
                           MARGIN = 1, 
                           FUN = mean,
                           na.rm = TRUE)
  
  ## invert error measures for computation convenience
  quality.out[,4:6] <- 1 / quality.out[,4:6]
 
  ## isolate quality measure of interest
  if(quality == "mRm") {
    quality.out.single <- quality.out[,1]
  } else if(quality == "mRn") {
    quality.out.single <- quality.out[,2]
  } else if(quality == "mRt") {
    quality.out.single <- quality.out[,3]
  } else if(quality == "mEm") {
    quality.out.single <- quality.out[,4]
  } else if(quality == "mEn") {
    quality.out.single <- quality.out[,5]
  } else if(quality == "mEt") {
    quality.out.single <- quality.out[,6]
  }
  
  ## identify l with optimal quality
  l.opt <- l[quality.out.single == max(quality.out.single, na.rm = TRUE)]
  
  ## optionally plot result
  if(plot == TRUE) {
    plot(x = l,
         y = quality.out.single,
         main = "Model quality due to weight transformation values",
         xlab = "l",
         ylab = quality)
    points(x = l.opt,
           y = quality.out.single[quality.out.single == 
                                    max(quality.out.single, na.rm = TRUE)],
           col = 2)
  }
  
  ## return result
  return(l.opt)
}
