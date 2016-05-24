#' The simple version of CUTOFF
#' @param data a data matrix with missing values
#' @param cutoff the cutoff value for the CUTOFF method
#' @param method CUTOFF method to be used. 
#' @param ID If the reference information needs to be retained during the imputation
#' if TRUE, then reference information can be retained from the returned list by 
#' calling ID. If FALSE, then no reference information will be retained.
#' @param ... other arguments
#' @return If ID = FALSE, then return the imputed data matrix with no 
#' missing values. If ID = TRUE, then return a list of two components:
#' 
#' \item{imputed}{The imputed data matrix with no missing values}
#' \item{ID}{The reference information during the imputation}
#' 
#' @references Lingbing Feng, Gen Nowak, Alan. H. Welsh and Terry. J. O'Neill 
#' (2014): CUTOFF: A Spatio-temporal Imputation Method, 
#'\emph{Journal of Hydrology}. (submitted)
#' 
#' @export
#' @examples
#' data(hqmr.data)
#' #' # check the number of missing values
#' nmissing(hqmr.data[, -79])
#' # impute the data by the CUTOFF method
#' impdata <- Cut(data = hqmr.data)
#' nmissing(impdata)
Cut <- function (data, cutoff = 0.75, method = "pearson", ID = FALSE, ...) {
    chunk <- as.matrix(subset(data, select = -date))
    cor_matrix <- cor(chunk, use = "complete.obs", method = method)
    diag(cor_matrix) <- 0
    dimnames(cor_matrix) <- NULL
    year <- 1900 + as.POSIXlt(data$date)$year
    month <- 1 + as.POSIXlt(data$date)$mon
    rainfall_new <- chunk
    nc <- ncol(rainfall_new)
    nr <- nrow(rainfall_new)
    id <- vector("list", nc)
    for (j in 1:nc) {
        if (any(cor_matrix[, j] >= cutoff)) {
            R_id <- which(cor_matrix[, j] >= cutoff)
        }
        else {
            big <- max(cor_matrix[, j])
            R_id <- which(cor_matrix[, j] == big)
        }
        for (i in 1:nr) {
            if (is.na(rainfall_new[i, j])) {
                C_id <- which(year != year[i] & month == month[i])
                C_bar_data <- chunk[C_id, j]
                C_bar <- mean(C_bar_data, na.rm = TRUE)
                R_bar_data <- chunk[C_id, R_id]
                R_bar <- mean(as.vector(R_bar_data), na.rm = TRUE)
                R_ix <- which(year == year[i] & month == month[i])
                R_data <- chunk[R_ix, R_id]
                if (!all(is.na(R_data))) {
                  R <- mean(as.vector(R_data), na.rm = TRUE)
                }
                else {
                  find.na <- function(x) {
                    for (i in 1:length(x)) while (!is.na(x[i])) return(x[i])
                  }
                  idj <- order(cor_matrix[, j], decreasing = TRUE)
                  rvec <- as.vector(chunk[R_ix, idj])
                  r <- find.na(rvec)
                  R <- r
                }
                rainfall_new[i, j] <- R * C_bar/R_bar
            }
            if (ID) {
                id[[j]] <- R_id
            }
        }
    }
    if (ID) {
        return(list(imputed = rainfall_new, ID = id))
    }
    else {
        return(rainfall_new)
    }
}
