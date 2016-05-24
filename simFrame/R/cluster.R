# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

# clusterAssign <- function(cl, x, value) {
#     ans <- clusterCall(cl, 
#         function(x, value) {
#             assign(x, value, envir=.GlobalEnv)
#             NULL
#         }, 
#         x=x, value=value)
#     invisible(ans)
# }
