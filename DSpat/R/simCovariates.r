simCovariates = function(hab.range=30, probs=c(1/3,2/3) , river.loc=50)
################################################################################
# Create a set of covariates in a 100x100 world with a vertical linear feature
# and discrete habitats
#
# Arguments:
#
#  hab.range   - habitat range that controls patchiness
#  probs       - ordered probablities that define habitat cutoffs
#  river.loc   - x coordinate for north-south river location
#
# Value:  dataframe with columns x,y,river and habitat #
#
# Devin Johnson; Jeff Laake
# 10 April 2008
################################################################################
{
 coords = cbind(x=rep(seq(.5,99.5,1),each=100), y=rep(seq(.5,99.5,1), 100))
 #habRaw = GaussRF(coords, model="gauss", grid=FALSE, param=c(0,1,0,))
 habRaw = RFsimulate(model=RMgauss(var=1, scale=hab.range/sqrt(-log(0.05))), x=coords[,"x"], y=coords[,"y"])@data[,1]
 quan = c(-Inf,quantile(habRaw, probs),Inf)
 habitat=as.numeric(cut(habRaw,quan))
 d = abs(coords[,'x']-river.loc)/max(abs(coords[,'x']-river.loc))
 return(as.data.frame(cbind(x=coords[,1], y=coords[,2], river=d, habitat=habitat)))
}
