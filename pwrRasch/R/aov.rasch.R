#' Three-Way Analysis of Variance with Mixed Classification for Testing the Rasch Model
#' 
#' This function applies the three-way analysis of variance with mixed classification 
#' for testing the Rasch model.
#' 
#' The F-test in a three-way analysis of variance design (A > \strong{B}) x C with mixed classification
#' (fixed factor A = subgroup, random factor B = testees, and fixed factor C = items) is used
#' to test the Rasch model. Rasch model fitting means that there is no interaction A x C.
#' A statistically significant interaction A x C indicates differential item functioning (DIF) 
#' of the items with respect of the two groups of testees Note, if a main effect of A (subgroup)
#' exists, an artificially high type I risk of the A x C interaction F-test results - that is, 
#' the approach works as long as no statistically significant main effect of A occurs.
#' Note that in case of unbalanced groups computation can take a long time.
#'
#' @param data        A data frame in which the variables specified in the model will be found.
#'                    Note that data needs to be in 'long' format. 
#' @param group       Column name of the data frame containing the grouping variable. 
#' @param person      Column name of the data frame containing the person number variable.
#' @param item        Column name of the data frame containing the item number variable. 
#' @param response    Column name of the data frame containing the response variable.
#' @param output      If \code{TRUE}, an output will be shown on the console.
#'
#' @author 
#' Takuya Yanagida \email{takuya.yanagida@@univie.ac.at},
#' Jan Steinfeld \email{jan.steinfeld@@univie.ac.at}
#'
#' @references
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model 
#' calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model.
#' \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' @return
#' Returns an ANOVA table
#' 
#' @seealso 
#' \code{\link{reshape.rasch}}, \code{\link{pwr.rasch}}
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # simulate Rasch model based data
#' # 100 persons, 20 items,
#' dat <- simul.rasch(100, items = seq(-3, 3, length.out = 20))
#' # reshape simulated data into 'long' format with balanced assignment 
#' # of testees into two subgroups
#' dat.long <- reshape.rasch(dat, group = rep(0:1, each = nrow(dat) / 2))
#' # apply three-way analysis of variance with mixed classification for testing the Rasch model
#' aov.rasch(dat.long) 
#' 
#' # extract variable names of items
#' vnames <- grep("it", names(aid_st2), value = TRUE)
#' # reshape aid subtest 2 data into 'long' format with split criterium sex
#' aid_long.sex <- reshape.rasch(aid_st2[, vnames], group = aid_st2[, "sex"])
#' # apply three-way analysis of variance with mixed classification for testing the Rasch model
#' aov.rasch(aid_long.sex)

#' }
aov.rasch <- function(data, group = "group", person = "person", item = "item", response = "response",
                      output = TRUE) {
  
  data.test <- data[!is.na(data$response), ]
  
  if (eval(parse(text = paste0("length(unique(table(data$", group, "))) == 1 &",
                               "all(tapply(data.test$", person, ", data.test$", item, ", function(x) all(table(x) == 1))) & ",
                               "length(unique(tapply(data.test$", item, ", data.test$", person, ", function(x) sum(table(x) == 1)))) == 1")))) {
    
     restab <- aov.rasch.balanced(data, group = group, person = person, item = item, response = response,
                                 output = output)
    
  } else {
    
    restab <- aov.rasch.unbalanced(data, group = group, person = person, item = item, response = response,
                                   output = TRUE)

  }

  class(restab) <- "aovrasch"
  return(invisible(restab))

}