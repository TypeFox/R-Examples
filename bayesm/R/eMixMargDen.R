eMixMargDen=
function(grid,probdraw,compdraw) 
{
#
# Revision History:
#   R. McCulloch 11/04
#
# purpose: plot the marginal density of a normal mixture averaged over MCMC draws
#
# arguments:
#    grid -- array of grid points, grid[,i] are ordinates for ith component
#    probdraw -- ith row is ith draw of probabilities of mixture comp
#    compdraw -- list of lists of draws of mixture comp moments (each sublist is from mixgibbs)
#
# output:
#    array of same dim as grid with density values
#
#
den=matrix(0,nrow(grid),ncol(grid))
for(i in 1:length(compdraw)) den=den+mixDen(grid,probdraw[i,],compdraw[[i]])
return(den/length(compdraw))
}
