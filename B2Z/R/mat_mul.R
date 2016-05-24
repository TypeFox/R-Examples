#################################################################
#Suppose in the likelihood (Y(ti)-log(C(theta1,ti))^T = (a b)^T #
#and S = Sigma^-1. This function computes (a b)^T%*%S%*%(a b)   #
#################################################################

mat_mul <- function(a, b, S) a*(a*S[1] + b*S[2]) + b*(a*S[3] + b*S[4])