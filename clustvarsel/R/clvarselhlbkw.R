#############################################################################
##  Headlong backward search
#############################################################################

clvarselhlbkw <- function(X, G = 1:9, 
                          emModels1 = c("E","V"), 
                          emModels2 = mclust.options("emModelNames"),
                          samp = FALSE, sampsize = 2000, 
                          hcModel = "VVV",
                          allow.EEE = TRUE, forcetwo = TRUE, 
                          BIC.upper = 0, BIC.lower = -10,
                          itermax = 100)
{
  warning("Headlong backward search is not currently available! Please use another type of search and/or direction.", immediate. = TRUE)
  invisible(NULL)
}
