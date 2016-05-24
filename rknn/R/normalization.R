################################################################################
# Normalization Functions                                                      #
# File:   Normalization.R                                                      #
# Author: Shengqiao Li                                                         #
# Date:   Octber 8, 2009  (initial)                                            #
# Dependency: None                                                             #
################################################################################
normalize.unit<- function(data)
{
  ####normalize to [0,1]
  minvect<- apply(data, 2, min);
  rangevect<- apply(data, 2, function(x)diff(range(x)))  
  scale(data, center=minvect, scale=rangevect)
}

normalize.sigmoidal<- function(data)
{
  #nonlinearly transform data into [-1,1] using a sigmoid function
  z<- scale(data);
  tanh(z/2)

}

normalize.softmax<- function(data)
{
  #more or less linear in the middle range, and has a nonlinearity at both ends
  z<- scale(data);
  plogis(z)
}
normalize.decscale<- function (data)
{
    #decimal scaling to a matrix or dataframe. Decimal scaling transforms the
    #data into [-1,1] by finding k such that the absolute value of the maximum
    #value of each attribute divided by 10^k is less than or equal to 1.
    
    maxvect <- apply(abs(data), 2, max)
    kvector <- ceiling(log10(maxvect))
    scalefactor <- 10^kvector
    scale(data, center = FALSE, scale = scalefactor)
}

################################################################################