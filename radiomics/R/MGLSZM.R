#' @include GLSZM.R
NULL

#' Multiple gray level size zone matrix.
#'
#' \code{mglszm} returns a matrix of class "mglszm", the multiple gray level size zone
#'  matrix for a given matrix.
#'
#' The function creates
#' a GLSZM using grey levels: 2, 4, 8, 16, 32, 64, 128, and 256. The values of these
#' GLSZM's are then weighted and combined using a gaussian distribution with mean 
#' of 0 and sd of 1.
#' 
#' 
#' Can be visualized using \code{image(mglszm(data))}. For visualization info
#' see \code{?image.radiomics}
#' 
#' @param data A 2D image matrix.
#' @param truncate Logical, removes any sizes or gray levels that have no entries.
#' @param ... Can be given verbose=FALSE to suppress output from the n_grey conversion.       
#' @return a matrix of dimension n_grey by region size, the MGLSZM. The column 
#'   names represent the region size, row names represent grey level, and 
#'   the entries represent the count of how many times a given size of given grey level
#'   occur.
#'   
#' @references \url{http://thibault.biz/Research/ThibaultMatrices/MGLSZM/MGLSZM.html}   
#' @importFrom stats dnorm
#' @examples
#' \dontrun{
#' image(psf)
#' mglszm(psf)
#' 
#' image(discretizeImage(psf, n_grey=5, verbose=F))
#' mglszm(psf, n_grey=5, verbose=F) 
#' }
#' @export
mglszm <- setClass("mglszm",
                   contains="matrix"
)


setMethod("initialize", 
          signature = "mglszm", 
          definition = function(.Object, data, truncate, ...){
            #TODO: Make weights a function argument
            #TODO: Make number of bits a function argument
            unique_vals <- unique(c(data))
            
            #pass through discretizeImage to access warnings and errors
            data <- discretizeImage(data, n_grey=length(unique_vals), ...)
            
            #If composed entirely of NA, return <0,0> matrix
            if(sum(is.na(data)) == dim(data)[1]*dim(data)[2]){
              .Object@.Data <- matrix()[-1,-1]
              return(.Object)
            }
            
            # Check if image consists of single pixel 
            # Multiple behaviours are potentially desirable -
            # Here I have chosen to return nothing, as the glszm will 
            # contain all relevant information already.
            #Another option would be to return the glszm. Unclear which is the better choice.
            if(length(sort(unique_vals)) == 1){
              .Object@.Data <- matrix()[-1,-1]
              return(.Object)
            }
            #create weights
            #  - -3.5 to 3.5 makes the sum of the weights ~1
            #  - 8 is the number of bits we will use (2^k)
            #  - 2, 4, 8, 16, 32, 64, 128, 256
            #  - Therefore 16 and 32 are the highest weighted
            weights <- dnorm(seq(-3.5,3.5,length.out=8), mean=0, sd=1)
            
            #initialize matrix at max grey level
            #such that it begins as the maximum number of possible rows
            all_grey <- as.character(1:prod(dim(data)))
            MGLSZM <- matrix(0, nrow=256, ncol=prod(dim(data)), dimnames=list(1:256, 1:prod(dim(data))))
            
            
            #Loop over bit values, gray levels = 2^k
            for(k in 1:8){
              n_grey = 2^k
              KGLSZM <- glszm(data, n_grey=n_grey, verbose=FALSE, ...)
              #Scale by weighting factor:
              KGLSZM <- weights[k] * KGLSZM 
              
              #Add on a new column for each connected size not represented in the glszm
              # TODO: This can probably be largely replaced using match()
              # See glrlm definition
              buffer_matrix <- matrix(0, nrow = nrow(KGLSZM),
                                      ncol = sum( ! all_grey %in% colnames(KGLSZM)),
                                      dimnames = list(
                                        row.names(KGLSZM),
                                        all_grey[which(! all_grey %in% colnames(KGLSZM))] ))
              
              KGLSZM <- cbind(KGLSZM, buffer_matrix)  
              KGLSZM <- KGLSZM[,order(as.numeric(colnames(KGLSZM)))]
              
              #KGLSZM must be expanded to be added to the MGLSZM
              rows <- sort(rep_len(1:nrow(KGLSZM), length.out=256))
              KGLSZM <- KGLSZM[rows, ]
              rownames(KGLSZM) <- 1:256
              
              MGLSZM <- MGLSZM + KGLSZM
            }
            
            if(truncate){
              truncated <- MGLSZM[which(rowSums(MGLSZM) > 0),which(colSums(MGLSZM) > 0)]
              if(is.matrix(truncated)){
                MGLSZM <- truncated
              }
            }

            .Object@.Data <- MGLSZM
            
            .Object
            
          }
)

#' @export          
mglszm <- function(data, truncate = TRUE, ...){
  return(new("mglszm", data, truncate, ...))
}
