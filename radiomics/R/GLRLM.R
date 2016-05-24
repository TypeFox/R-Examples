#' Gray level run length matrix.
#'
#' \code{glrlm} returns a gray level run length matrix for a given matrix.
#'
#'
#' Can be visualized using \code{image(glrlm(data))}. For visualization info
#' see \code{?image.radiomics}
#' 
#' @param data A numeric 2D matrix.
#' @param angle One of 0, 45, 90 or 135, the direction the run is calculated.
#' @param n_grey an integer value, the number of grey levels the image should
#'   be quantized into.
#' @param max_run_length An integer value, the default is the maximum possible
#'   run length. Setting it to a smaller value truncates the output. Desirable
#'   in cases where the matrix is extremely sparse, for example when
#'   there are few long runs.
#' @param truncate Logical Remove run lengths which have no entries 
#' @param ... Can be given verbose=FALSE to suppress output from the n_grey conversion.       
#' @return a matrix of class "glrlm" of dimension n_grey by run length. The column 
#'   names represent the length of the run, and row names represent 
#'   grey values in the image.
#'   
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0146664X75800086}
#' @examples
#' \dontrun{
#' hallbey
#' glrlm(hallbey)
#' glrlm(hallbey, angle="90") 
#' }
#' @export
glrlm <- setClass("glrlm",
                  contains="matrix"
)


setMethod("initialize", 
          signature = "glrlm", 
          definition = function(.Object, data, angle, n_grey, max_run_length, truncate, ...){
            
            #send to discretizeImage for error checking
            data <- discretizeImage(data, n_grey=n_grey, ...)

            #initialize rlm
            unique_vals <- sort(unique(c(data)))
            rlm <- matrix(0, nrow=length(unique_vals), ncol=max_run_length)
            rownames(rlm) <- unique_vals
            colnames(rlm) <- c( 1:max_run_length )
            
            if(identical( angle, 45) | identical(angle, 135) ){
              if(identical(angle, 45)) data <- data[ nrow(data):1, ]
              
              # check rle of each diagonal of the matrix, add to rlm
              # start at bottom left, work to top right
              bottom <- nrow(data) - 1
              top <- -ncol(data) + 1
              for(i in bottom:top){
                
                runs <- t(table(rle(data[row(data)==(col(data) + i)])))
                rlm <- add_to_rlm(runs, rlm, max_run_length)
                
              }
              
            } else if (identical(angle, 0) | identical(angle, 90)){
              if(identical(angle, 90)) data <- t(data)
              
              for(i in 1:nrow(data)){
                
                runs <- t(table(rle(data[i,])))
                rlm <- add_to_rlm(runs, rlm, max_run_length)
                
              }
              
            } else {
              stop("Angle must be one of '0', '45', '90', '135'.")
            }
            
            if(truncate){
              truncated_rlm <- rlm[,which(colSums(rlm) > 0)]
              if(is.matrix(truncated_rlm)){
                rlm <- truncated_rlm
              }
            }
            .Object@.Data <- rlm
            
            .Object
          }
          
          
)

#' @export          
glrlm <- function(data, angle = 0, n_grey = 32, max_run_length = min(dim(data)), truncate=TRUE, ...){
  return(new("glrlm", data, angle, n_grey, max_run_length, truncate, ...))
}

add_to_rlm <- function(runs, rlm, max_run_length){
  #Intermediate function, not meant to be called directly
  # adds matching rows and columns from rle tables to rlm
  
  mrow <- match(rownames(runs), rownames(rlm))
  mrow <- mrow[which(!is.na(mrow))] 
  mcol <-  match(colnames(runs), colnames(rlm))
  mcol <- mcol[which(!is.na(mcol))] 
  
  rlm[mrow,mcol] <- rlm[mrow,mcol] + runs[,which(as.numeric(colnames(runs)) <= max_run_length)]
  return(rlm)
}
