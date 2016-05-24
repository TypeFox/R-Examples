#----------------------------------------------------------
#
# Calculate RCRV reduced centered random variable
# with bias and dispersion
# Return:
#     bias: Bias 
#     disp: dispersion
#     y : vector of y used to calculate bias and dispersion
#     obsError: observation error used in the calculation 
# Author: Ronald Frenette, Severe Weather Lab, Quebec region
#         Jun 2009
# 
#-----------------------------------------------------------

rcrv<-function(obs,epsMean,epsVariance,obsError)
{
  y<-(obs-epsMean)/(sqrt(epsVariance+(obsError*obsError)))
  bias<-mean(y) 
  disp<-sqrt(var(y))
  return(list(bias=bias,disp=disp,y=y,obsEror=obsError))
}

