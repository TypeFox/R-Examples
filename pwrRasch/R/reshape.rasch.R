#' Reshape data frame in wide format into a long format
#' 
#' This function reshapes a matrix from 'wide' into a 'long' format. This is necessary for
#' the three-way analysis of variance with mixed classification for testing the Rasch model.
#'
#' In order to apply the three-way analysis of variance with mixed classification for
#' testing the Rasch model, data need to be in 'long' format. That is, Rasch model 
#' data design is interpreted as a analysis of variance design (A > \strong{B}) x C, 
#' where items are levels of a fixed factor C and the testees are levels of a random 
#' factor B, nested within a fixed factor A of different subgroups.
#'  
#' @param data    Matrix or data frame in 'wide' format.
#' @param group   Vector which assigns each person to a certain subgroup (external split criterion). 
#'                Note, that this function is restricted to A = 2 subgroups.
#' 
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#' 
#' @export
#' 
#' @seealso 
#' \code{\link{aov.rasch}}
#'
#' @references
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model 
#' calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model.
#' \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' @return 
#' Returns a data frame with following entries:
#' \tabular{ll}{ 
#'   \code{group}    \tab fixed factor A (subgroup) \cr
#'   \code{person}   \tab random factor B (testees) \cr
#'   \code{item}     \tab fixed factor C (items) \cr
#'   \code{response} \tab dependent variable, 0 (item not solved) and 1 (item solved) \cr
#' }
#' 
#' @examples
#' \dontrun{
#' 
#' # simulate Rasch model based data
#' # 100 persons, 20 items,
#' dat <- simul.rasch(100, items = seq(-3, 3, length.out = 20))
#' # reshape simulated data into 'long' format with balanced assignment 
#' # of testees into two subgroups.
#' dat.long <- reshape.rasch(dat, group = rep(0:1, each = nrow(dat) / 2))
#' head(dat.long)
#' 
#' # extract variable names of items
#' vnames <- grep("it", names(aid_st2), value = TRUE)
#' # reshape aid subtest 2 data into 'long' format with split criterium sex
#' aid_long.sex <- reshape.rasch(aid_st2[, vnames], group = aid_st2[, "sex"])
#' }
reshape.rasch <- function(data, group) {

  #--------------------------------------------------------------------------------------------------------#
  # Input Check
  
  # A = 2 subgroups
  if (length(unique(group)) != 2) {
    
    stop("Subgroups A != 2 (external split criterion) specified")
    
  } 
  
  # length of split criterion
  if (length(group) != nrow(data)) {
    
    stop("Length of A (external split criterion) does not match the number of observations in the data")
    
  }  
  
  #--------------------------------------------------------------------------------------------------------#
  
  resp1 <- data[which(group == unique(group)[1]), ]
  resp2 <- data[which(group == unique(group)[2]), ]
  
  b1 <- nrow(resp1)
  b2 <- nrow(resp2)
  
  c <- ncol(data)   
  
  data.long <- data.frame(group = c(rep(1, times = b1*c), rep(2, times = b2*c)),
                          person = c(rep(1:b1, times = c, each = 1), rep((b1 + 1):(b1 + b2), times = c, each = 1)),
                          item = c(rep(1:c, each = b1), rep(1:c, each = b2)),
                          response = unlist(c(resp1, resp2)))                                      
  
  return(data.long)
  
}