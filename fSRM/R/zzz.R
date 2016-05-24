#' @importFrom grDevices gray
#' @importFrom methods show
#' @importFrom stats addmargins formula pchisq
#' @importFrom utils packageDescription packageVersion read.csv read.table str


# Define skinning variables in a new environment
#------------------------------------
style <- new.env(parent=globalenv())
style$actor <- "A"
style$partner <- "P"
style$relationship <- "R"
style$familyeffect <- "FE"
style$self <- "S"