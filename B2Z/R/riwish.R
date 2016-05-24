###########################################
#This function randomly generates a       # 
#value from an inverse Wishart. Here,     #
#instead of returning a matrix, it returns# 
#a vector of size 3 (tauN, tauF, tauNF).  #
#                                         #
###########################################

riwish <-
function (i, v, S) 
{
#This function is from the package MCMCpack by  
#Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park
mat <- solve(rwish(v, solve(S)))
return( c(mat[1,1], mat[2,2], mat[1,2]) )
}

