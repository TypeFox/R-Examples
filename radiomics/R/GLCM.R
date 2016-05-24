#' Gray level co-occurrence matrix.
#'
#' \code{glcm} returns a gray level co-occurrence matrix for a given matrix.
#' 
#' 
#' Can be visualized using \code{image(glcm(data))}. For visualization info
#' see \code{?image.radiomics}
#'
#' 
#' @param data A numeric 2D matrix.
#' @param angle One of "0", "45", "90" or "135", the pixel to which the 
#'   current pixel is compared.
#' @param d an integer value, the distance between the current pixel, and the
#'   pixel to which it is compared.
#' @param n_grey an integer value, the number of grey levels the image should
#'   be quantized into. If greater than the number of unique values in the image,
#'   no action will be taken.
#' @param normalize Logical value, if TRUE (default) the matrix will be normalized such that 
#'   the sum of it's components is 1.
#' @param ... Can be given verbose=FALSE to suppress output from the n_grey conversion.       
#' @return a matrix of dimension n_grey by n_grey, the GLCM. The column and row names represent 
#'   grey values in the image.
#'   
#' @references \url{http://www.fp.ucalgary.ca/mhallbey/tutorial.htm}
#' @examples
#' \dontrun{
#' hallbey
#' glcm(hallbey)
#' glcm(hallbey, angle="90") #vertical GLCM
#' }
#' @export
glcm <- setClass("glcm",
                 contains="matrix"
)

setMethod("initialize", 
          signature = "glcm", 
          definition = function(.Object, data, angle, d, n_grey, normalize, ...){
            
            #Send to discretizeImage for error checking
            #Discretize grey values if required
            #discretize image and initialize GLCM based on discretized image
            data <- discretizeImage(data, n_grey=n_grey, ...)
            
            unique_vals <- sort(unique(c(data)))
            
            #If matrix is composed of a single value glcm is undefined
            #return empty matrix
            if(length(unique_vals) == 1 && length(which(data == unique_vals)) == 1){
              .Object@.Data <- matrix(NA)[-1,-1]
              return(.Object)
            }
            
            #the value of 0 is reserved for NAs in the matrix, 
            #if there are any 0's in the DF, add 1 to all values
            #original values will be replaced after
            if(is.element(0, data)) data <- data + 1
            
            #Convert All NAs to 0
            data[is.na(data)] <- 0
            
            
            if(identical(angle, 0)){
              counts <- glcm0(data, n_grey = max(data), d)
              
            } else if (identical(angle, 45)){
              counts <- glcm45(data, n_grey = max(data), d)
              
            } else if (identical(angle, 90)){
              counts <- glcm90(data, n_grey = max(data), d)
              
            } else if (identical(angle, 135)){
              counts <- glcm135(data, n_grey = max(data), d)
              
            } else {
              stop("angle must be one of '0', '45', '90', '135'.")
            }
                            
            #Row 1 and Col 1 hold NA values, remove them
            counts <- counts[-1, -1]
            
            #Situation where matrix is composed of a single NA
            if(length(counts) == 0){
              .Object@.Data <- counts
              return(.Object)
            }

            #Replace proper values in column and row names
            #Two situations:
            #1. No zeroes were present, thus nothing was added
            #2. One was added to all entries because there were zeros in the matrix
            if(is.matrix(counts)){
              if(dim(counts)[1] == max(unique_vals)){ #ie. 1 wasn't added
                counts <- counts[unique_vals, unique_vals]
                #counts <- counts[which(rownames(counts) %in% unique_vals), which(colnames(counts) %in% unique_vals)]
              } else if (dim(counts)[1] == max(unique_vals)+1) {
                #counts <- counts[which((as.numeric(rownames(counts)) - 1) %in% unique_vals), which((as.numeric(colnames(counts)) - 1) %in% unique_vals)]
                counts <- counts[unique_vals + 1, unique_vals + 1]
              }  
            }
            
            if(!is.matrix(counts)) {
              #Edge case where only a single grey value present - leads to a numeric, rather than a matrix
              #Therefore case to 1x1 matrix
              counts <- matrix(counts)
            }
            
            rownames(counts) <- colnames(counts) <- unique_vals
                      
            #GLCMs should be symmetrical, so the transpose is added
            counts <- counts + t(counts)
            #Normalize
            if(normalize){
              count_sum <- sum(counts)
              if(count_sum > 0){
                counts <- counts/count_sum
              }
            }
              
            
            
            .Object@.Data <- counts
            
            .Object
            
          }   )

#' @export          
glcm <- function(data, angle = 0, d=1, n_grey = 32, normalize=TRUE, ...){
  return(new("glcm", data, angle, d, n_grey, normalize, ...))
}
