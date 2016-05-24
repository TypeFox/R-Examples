#' Gray level size zone matrix.
#'
#' \code{glszm} returns a gray level size zone matrix for a given matrix.
#'
#'
#' Can be visualized using \code{image(glszm(data))}. For visualization info
#' see \code{?image.radiomics}
#' 
#' @param data A numeric 2D matrix.
#' @param n_grey an integer value, the number of grey levels the image should
#'   be quantized into.
#' @param truncate Logical. Remove values for sizes that have no entries
#' @param ... Can be given verbose=FALSE to suppress output from the n_grey conversion. 
#' @return a matrix of dimension n_grey by region size, the GLSZM. The column 
#'   names represent the region size, row names represent grey level, and 
#'   the entries represent the count of how many times a given size of given grey level
#'   occur.
#'
#' @examples
#' \dontrun{
#' image(psf)
#' glszm(psf)
#' 
#' image(discretizeImage(psf, n_grey=5, verbose=F))
#' glszm(psf, n_grey=5, verbose=F) 
#' }
#' @references \url{http://thibault.biz/Research/ThibaultMatrices/GLSZM/GLSZM.html}
#' @importFrom reshape2 acast
#' @importFrom spatstat as.im
#' @importFrom spatstat levelset
#' @importFrom spatstat connected
#' @importFrom methods new 
#' @export
glszm <- setClass("glszm",
                  contains="matrix"
)


setMethod("initialize", 
          signature = "glszm", 
          definition = function(.Object, data, n_grey, truncate, ...){
            #send to discretizeImage for error checking
            data <- discretizeImage(data, n_grey=n_grey, ...)
            if(sum(is.na(data)) == dim(data)[1]*dim(data)[2]){
              .Object@.Data <- matrix()[-1,-1]
              return(.Object)
            } 
            grey_lvls <- unique(c(data))
            grey_lvls <- grey_lvls[!is.na(grey_lvls)]
            #convert to data for use with spatstats functions
            data <- spatstat::as.im(data)
            
            #Initialize dataframe to hold count data
            count_data <- data.frame()
            
            
            for(i in grey_lvls){
              # Threshold the data
              imBinary <- spatstat::levelset(data, i, compare="==")
              connections <- spatstat::connected(imBinary)
              
              # Extract counts of each uniqe value 
              counts <- table(table(as.matrix(connections)))
              count_data <- rbind(count_data, data.frame(i, counts))
            }
            
            #Clean up names 
            colnames(count_data) <- c("greylvl", "size", "counts")
            #cast to matrix
            count_data <- reshape2::acast(count_data, greylvl~size, value.var="counts")
            #sort columns, if there is only a single size a vector is returned, hence the if
            if(length(colnames(count_data)) > 1 && nrow(count_data) > 1){
              count_data <- count_data[,order(as.numeric(as.character(colnames(count_data))))]
            }
            count_data[is.na(count_data)] <- 0
            
            if(truncate){
              truncated <- count_data[,which(colSums(count_data) > 0)]
              if(is.matrix(truncated)){
                count_data <- truncated
              }
            }
            
            .Object@.Data <- count_data
            
            .Object
          }
)

#' @export          
glszm <- function(data, n_grey = 32, truncate=TRUE, ...){
  return(new("glszm", data, n_grey,truncate, ...))
}
