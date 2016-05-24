#
# model settings
# mean: E, D
# between: E, D 
# within: E, D
# cov: 0, D        # 0 - zero diagonal, D - same as variance

getModelStructure <- function(mean = "D", between = "D", within = "D", cov = "D") {
  if (any(!(mean %in% c("D", "E"))) | 
      any(!(between %in% c("D", "E"))) | 
      any(!(within %in% c("D", "E"))) | 
      any(!(cov %in% c("D", "0"))) | 
      length(mean) == 0 |
      length(between) == 0 |
      length(within) == 0 |
      length(cov) == 0) {
         stop("Wrong argument specification (mean, between, within should be %in% c('D','E'), cov %in% c('D', '0'))")
      }
  list(mean = mean, between = between, within = within, cov = cov)
}
