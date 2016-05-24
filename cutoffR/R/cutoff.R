#' The CUTOFF Spatio-temporal Imputation Method
#' 
#' @details This function implements the CUTOFF spatio-temporal imputation 
#' method that is described in Feng et al.(2014)
#' @param data a matrix or data frame with missing values
#' @param N a number indicating the number used for the "CUTOFF by number" method
#' @param cutoff a number indicating the cutoff value used fot he "CUTOFF by correlation" method
#' @param P a number for the "penalty" imputation option for CUTOFF. That is, for those candidate missing
#' station with too many reference stations, we can penalise and fix the number of reference stations to P
#' @param M a number used for the "relaxation" imputation option for CUTOFF. That is, for those cadidate
#' missing station with too few reference stations, we can relax and add its number of reference stations to M
#' @param Adj a number used for the "adjacent" method in CUTOFF. That is, the missing value's adjacent points in time
#' is also used for imputation. The default values is 1. 2 is also avaialbe. Any number bigger than 2 has not been implemented
#' yet. This options is useful when the length of the time series is short so may be more temporal information can be useful 
#' to improve the imputation performance.
#' @param space.weight a logical value. If true, then space weighting strategy is carried out. The default is FALSE.
#' @param method the imputation method to be used. There are three options: "correlation", "number" and "penalty". Details
#' can be found in Feng et al.(2014).
#' @param time.opts options for the temporal dimension; either "average" or "adjacent" can be used. "average" refers to
#' simple averaging, "adjacent" refers to the "aajacent" method.
#' @param kernel logical, if TRUE then kernel smoothing can be used to smooth the averaging. Default is FALSE. If TURE, then kerFUN
#' has to be specified. 
#' @param kerFUN the kernel function to be used for kernel smoothing. There are four kernel functions available in this package: Epank, UnifK
#' GaussK and CosK. User can define their own kernel function to pass to this function.
#' @param lambda a number indicating the bandwidth parameter value for kernel smoothing. 
#' @param corr the type of correlation coefficient to be used for the "CUTOFF by correlation" method. Default is "pearson", 
#' "spearman" and "kendall" are alternatives.
#' @param keep.ID if the reference ID for each missing stations need to be kept. If TRUE, relevant ID information can be retrieved after imputation.
#' Default is FALSE.
#' @param ... other arguments that can passed
#' @return If keep.ID = FALSE, then return the imputed data matrix with no 
#' missing values. If keep.ID = TRUE, then return a list of two components:
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
#' # check the number of missing values
#' nmissing(hqmr.data[, -79])
#' # impute the data by the CUTOFF method
#' impdata <- cutoff(data = hqmr.data)
#' nmissing(impdata)
cutoff <- function (data, N = 4, cutoff = 0.75, P = 5, M = floor(P/2), Adj = 1, 
                    space.weight = FALSE, method = c("correlation", "number", "penalty"), 
                    time.opts = c("average", "adjacent"), 
                    kernel = FALSE, kerFUN = NULL, lambda = NULL, 
                     corr = "pearson", keep.ID = FALSE, ...) {
  method <- match.arg(method)
  time.opts <- match.arg(time.opts)
  if (kernel) {
    if (is.null(kerFUN) || is.null(lambda)){
      stop("Kernel function and window size have to be set beforehand if you want to use kernel methods")
    }
    kerF <- match.fun(kerFUN)
  }
   
  chunk <- as.matrix(subset(data, select = - date))
  cor_matrix <- cor(chunk, use = "complete.obs", method = corr)
  diag(cor_matrix) <- 0
  dimnames(cor_matrix) <- NULL
  year <- 1900 + as.POSIXlt(data$date)$year
  month <- 1 + as.POSIXlt(data$date)$mon
  
  rainfall_new <- chunk
  nc <- ncol(rainfall_new)
  nr <- nrow(rainfall_new)
  if (keep.ID) {
    # keep the reference file information for each station
    id <- vector("list", nc)
  }
    for (j in 1:nc) {
      
      if (method == "correlation") {
        
        if (any(cor_matrix[, j] >= cutoff)) {
          R_id <- which(cor_matrix[, j] >= cutoff)
        }
        else {
            big <- max(cor_matrix[, j])
            R_id <- which(cor_matrix[, j] == big)
        }
        
      } else if (method == "number") {
          big <- order(cor_matrix[, j], decreasing = TRUE)
          R_id <- big[1:N]
      } else if (method == "penalty") { 
        big <- order(cor_matrix[, j], decreasing = TRUE)
        if (any(cor_matrix[, j] >= cutoff)) {
          R_id <- which(cor_matrix[, j] >= cutoff)
          }  
        if (length(R_id) > P) {            
          R_id <- big[1:P]
          } else { 
            R_id <- big[1:M]
          }
      }
      if (space.weight) {
        
        w.Fun <- function(x, weight){
          return(x * weight/sum(weight))
        }
        
        cors <- cor_matrix[R_id, j] # correlation vector for reference stations
        len <- length(cors)
        if (len > 2) {
          wvec <- ((len - 2) * cors^2 / (1 - cors^2))
          } else {
          wvec <- rep(1, len)
        }
      }
    
        for (i in 1:nr) {
            if (is.na(rainfall_new[i, j])) {
              if (time.opts == "average") {
                a_id <- which(year != year[i] & month == month[i])
              }
              if (time.opts == "adjacent") {
                if (is.null(Adj)) {
                  cat("Adjacent value has not been set and 1 will be chosen as default", "\n")
                  Adj  <- 1
                }
                if (Adj == 1) {
                  
                  if (month[i] != 1 & month[i] != 12) {
                    a_id <- which((month %in% c(month[i], month[i] - 1, month[i] + 1)) & year != year[i])
                    } else if (month[i] == 1) {
                      a_id <- which((month %in% c(1, 2, 12)) & year != year[i])
                      } else {
                        a_id <- which((month %in% c(12, 1, 11)) & year != year[i])
                      }
                }
                if (Adj == 2) {
                  # adjacent positions 
                  if (month[i] != 1 & month[i] != 2  & month[i] != 11 & month[i] != 12) {
                    a_id <- which((month %in% c(month[i], month[i] - 1, month[i] - 2,  month[i] + 1, month[i] + 2)) & year != year[i])
                    } else if (month[i] == 1) {
                      a_id <- which((month %in% c(1, 2, 3, 11, 12)) & year != year[i])
                      } else if (month[i] == 2) {
                        a_id <- which((month %in% c(1, 2, 3, 4, 12)) & year != year[i])
                        } else if (month[i] == 11) {
                          a_id <- which((month %in% c(9, 10, 11, 12, 1)) & year != year[i])
                        }
                  else {
                    a_id <- which((month %in% c(10, 11, 12, 1, 2)) & year != year[i])
                  }
                }
              }
              
              C_bar_data <- chunk[a_id, j]     
              C_bar <- mean(C_bar_data, na.rm = TRUE)
              R_bar_data <- chunk[a_id, R_id]
              R_bar <- mean(as.vector(R_bar_data), na.rm = TRUE)
              if(space.weight && !is.null(dim(R_bar_data))) {              
                R_bar <- mean(apply(R_bar_data, 1, w.Fun, weight = wvec), na.rm = TRUE)
              }
              
              if (kernel) {
                yd <- abs(year - year[i])    
                yd_i <- yd[a_id]  
                z <- yd_i / lambda
                ker <- kerF(z)
                weightK <- function(x) sum(ker * x, na.rm = TRUE) / sum(ker)

                if (time.opts == "average") {
                  
                  C_bar_data <- chunk[a_id, j]
                  C_bar <- weightK(C_bar_data)
                  R_bar_data <- chunk[a_id, R_id]
                  if(is.null(dim(R_bar_data))){
                    R_bar <- weightK(R_bar_data)
                  } else {
                    R_bar <- mean(apply(R_bar_data, 2, weightK))
                    if(space.weight) {
                      R_bar <- mean(apply(R_bar_data, 1, w.Fun, weight = wvec), na.rm = TRUE)
                    }
                  }
                }
                if (time.opts == "adjacent") {
                  
                  split_fun <- function(x) {
                    splitdata <- split(x, as.factor(year[a_id]))
                    meandata <- sapply(splitdata, mean, na.rm = TRUE)
                    return(meandata)
                  }
            
                  C_bar_data <- split_fun(chunk[a_id, j])
                  C_bar <- weightK(C_bar_data)
                  R_bar_data <- sapply(as.data.frame(chunk[a_id, R_id]), split_fun)
                  if(is.null(dim(R_bar_data))){
                    R_bar <- weightK(R_bar_data)
                    } else {
                      R_bar <- mean(apply(R_bar_data, 2, weightK))
                      if(space.weight) {
                        R_bar <- mean(apply(R_bar_data, 1, w.Fun, weight = wvec), na.rm = TRUE)
                      }
                    }
                }
                
              }
              R_ix <- which(year == year[i] & month == month[i])
              R_data <- chunk[R_ix, R_id]
              if (space.weight) {
                R_data <- w.Fun(R_data, weight = wvec)
              }
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
            if (keep.ID) {
              id[[j]] <- R_id
            }
        }
    }
  if (keep.ID) {
    return(invisible(list(imputed = rainfall_new, ID = id)))
  }
  else {
    return(invisible(rainfall_new))
  }
}
