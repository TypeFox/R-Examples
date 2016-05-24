###############################################################################
## Example: Neg Binomial Family
###############################################################################

require(ROptEst)
options("newDevice"=TRUE)

#-------------------------------------------------------------------------------
## Preparations
#-------------------------------------------------------------------------------
## generates Neg.Binomial Family with
## m = 25 and probability of success theta = 0.25
N <- NbinomFamily(size = 25, prob = 0.25) 
N       # show N 

#An object of class "NbinomFamily"
#### name:       Negative Binomial family
#
#### distribution:       Distribution Object of Class: Nbinom
# size: 25
# prob: 0.25
#
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### props:
#[1] ""
#
#plot(N) # plot of Nbinom(size = 25, prob = 0.25) and L_2 derivative
#checkL2deriv(N)
#
##precision of centering:  -3.103725e-15
##precision of Fisher information:
##              prob
##prob -2.842171e-14
##$maximum.deviation
##[1] 2.842171e-14

#-------------------------------------------------------------------------------
## classical optimal IC
#-------------------------------------------------------------------------------
IC0 <- optIC(model = N, risk = asCov())
IC0       # show IC

#An object of class “IC” 
#### name:        Classical optimal influence curve for Negative Binomial family 
#### L2-differentiable parametric family:         Negative Binomial family 
#
#### 'Curve':    An object of class “EuclRandVarList” 
#Domain: Real Space with dimension 1 
#[[1]]
#length of Map:   1 
#Range:  Real Space with dimension 1 
#
#### Infos:
#     method  message                                  
#[1,] "optIC" "optimal IC in sense of Cramer-Rao bound"

plot(IC0) # plot IC
checkIC(IC0)

#precision of centering:  7.796853e-15 
#precision of Fisher consistency:
#              prob
#prob -2.166600e-12
#maximum deviation 
#     2.166600e-12 

Risks(IC0)

#$asCov
#         prob
#prob 0.001875
#
#$trAsCov
#[1] 0.001875

#-------------------------------------------------------------------------------
## lower case radius
#-------------------------------------------------------------------------------
lowerCaseRadius(L2Fam = N, neighbor = ContNeighborhood(), risk = asMSE())

#lower case radius
#         4.153322

lowerCaseRadius(L2Fam = N, neighbor = TotalVarNeighborhood(), risk = asMSE())

#lower case radius
#         1.840705 
         
#-------------------------------------------------------------------------------
## L_2 family + infinitesimal neighborhood
#-------------------------------------------------------------------------------
RobN1 <- InfRobModel(center = N, neighbor = ContNeighborhood(radius = 0.5))
RobN1     # show RobN1

#An object of class “InfRobModel” 
####### center:  An object of class "NbinomFamily"
#### name:       Negative Binomial family
#
#### distribution:       Distribution Object of Class: Nbinom
# size: 25
# prob: 0.25
#
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### props:
#[1] ""
#
####### neighborhood:    An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5 

(RobN2 <- InfRobModel(center = N, neighbor = TotalVarNeighborhood(radius = 0.5)))

#An object of class “InfRobModel” 
####### center:  An object of class "NbinomFamily"
#### name:       Negative Binomial family
#
#### distribution:       Distribution Object of Class: Nbinom
# size: 25
# prob: 0.25
#
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### props:
#[1] ""
#
####### neighborhood:    An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.5 

#-------------------------------------------------------------------------------
## OBRE solution
#-------------------------------------------------------------------------------

system.time(ICA <-  optIC(model=RobN1, risk=asAnscombe(),
            verbose=TRUE,lower=NULL,upper=10))

#-------------------------------------------------------------------------------
## MSE solution
#-------------------------------------------------------------------------------
system.time(IC1 <- optIC(model=RobN1, risk=asMSE()))

#   user  system elapsed
#  10.53    0.02   11.21 

IC1

#An object of class “ContIC” 
#### name:        IC of contamination type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clip:       [1] 0.06143505
#### cent:       [1] 0.003717104
#### stand:
#           prob
#prob 0.00310076
#
#### Infos:
#     method  message                          
#[1,] "optIC" "optimally robust IC for ‘asMSE’"

checkIC(IC1)

#precision of centering:  4.297326e-14 
#precision of Fisher consistency:
#             prob
#prob 2.555733e-13
#maximum deviation 
#     2.555733e-13 

Risks(IC1)

#$asCov
#            prob
#prob 0.002157193
#
#$asBias
#$asBias$value
#[1] 0.06143505
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "ContNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.002157193
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#           prob
#prob 0.00310076
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5 


getRiskIC(IC1, asBias(), ContNeighborhood()) # standardized bias

#$asBias
#$asBias$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) convex contamination neighborhood"
#
#$asBias$value
#[1] 0.06143505


getRiskIC(IC1, asMSE(), ContNeighborhood(radius = 0.5))

#$asMSE
#$asMSE$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) convex contamination neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5
#
#$asMSE$value
#[1] 0.00310076


(Cov1 <- getRiskIC(IC1, asCov()))

#$asCov
#$asCov$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asCov$value
#            prob
#prob 0.002157193



(mse1 <- getRiskIC(IC1, asMSE(), TotalVarNeighborhood(radius = 0.5)))

#$asMSE
#$asMSE$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5
#
#$asMSE$value
#[1] 0.00310076

(bias1 <- getRiskIC(IC1, asBias(), TotalVarNeighborhood()))

#$asBias
#$asBias$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$value
#[1] 0.06143505


## only suboptimal -> ToDo-List
addRisk(IC1) <- list(Cov1, mse1, bias1)
Risks(IC1)

#$asCov
#$asCov[[1]]
#[1] 0.002157193
#
#$asCov$asCov
#$asCov$asCov$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asCov$asCov$value
#            prob
#prob 0.002157193
#
#
#
#$asBias
#$asBias$value
#[1] 0.06143505
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "ContNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#$asBias$asBias
#$asBias$asBias$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asBias$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$asBias$value
#[1] 0.06143505
#
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.002157193
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#           prob
#prob 0.00310076
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5 
#
#$asMSE$asMSE
#$asMSE$asMSE$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asMSE$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$asMSE$radius
#[1] 0.5
#
#$asMSE$asMSE$value
#[1] 0.00310076


plot(IC1)

system.time(IC2 <- optIC(model=RobN2, risk=asMSE()))

#   user  system elapsed
# 75.57    0.22   81.59

IC2

#An object of class “TotalVarIC” 
#### name:        IC of total variation type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clipLo:     [1] -0.06032159
#### clipUp:     [1] 0.052032
#### stand:
#            prob
#prob 0.005585238
#
#### Infos:
#     method  message                          
#[1,] "optIC" "optimally robust IC for ‘asMSE’"

checkIC(IC2)

#precision of centering:  9.638494e-16 
#precision of Fisher consistency:
#             prob
#prob 2.164935e-13
#maximum deviation 
#     2.164935e-13 

Risks(IC2)

#$asCov
#            prob
#prob 0.008255488
#
#$asBias
#$asBias$value
#[1] 0.1123536
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "TotalVarNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.008255488
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#           prob
#prob 0.01141132
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.5 
#

getRiskIC(IC2, asMSE(), TotalVarNeighborhood(radius = 0.5))

#$asMSE
#$asMSE$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5
#
#$asMSE$value
#[1] 0.01141132

getRiskIC(IC2, asBias(), TotalVarNeighborhood())

#$asBias
#$asBias$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$value
#[1] 0.1123536

getRiskIC(IC2, asBias(), ContNeighborhood())

#$asBias
#$asBias$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) convex contamination neighborhood"
#
#$asBias$value
#[1] 0.1123536

Cov2 <- getRiskIC(IC2, asCov())
addRisk(IC2) <- Cov2
Risks(IC2)

#$asCov
#$asCov[[1]]
#[1] 0.008255488
#
#$asCov$asCov
#$asCov$asCov$distribution
#[1] "Nbinom(25, 0.25)"
#
#$asCov$asCov$value
#            prob
#prob 0.008255488
#
#
#
#$asBias
#$asBias$value
#[1] 0.1123536
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "TotalVarNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.008255488
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#           prob
#prob 0.01141132
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.5 


plot(IC2)

#-------------------------------------------------------------------------------
## lower case solutions
#-------------------------------------------------------------------------------
(IC3 <- optIC(model=RobN1, risk=asBias()))


#An object of class “ContIC” 
#### name:        IC of contamination type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clip:       [1] 0.05458835
#### cent:       [1] 1.333333
#### stand:
#     [,1]
#[1,]    1
#### lowerCase:  [1] -0.3268154
#
#### Infos:
#     method  message                                        
#[1,] "optIC" "minimum asymptotic bias (lower case) solution"


checkIC(IC3)

#precision of centering:  7.555495e-16 
#precision of Fisher consistency:
#             prob
#prob 8.881784e-16
#maximum deviation 
#     8.881784e-16 

Risks(IC3)

#$asBias
#$asBias$value
#[1] 0.05458835
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "ContNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$asCov
#[1] 0.002918187
#
#$trAsCov
#$trAsCov$value
#[1] 0.002918187
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#[1] 0.00310443
#
#$asMSE$r
#[1] 0.25
#
#$asMSE$at
#An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5 


plot(IC3)

(IC4 <- optIC(model=RobN2, risk=asBias()))

#An object of class “TotalVarIC” 
#### name:        IC of total variation type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clipLo:     [1] -0.0574604
#### clipUp:     [1] 0.05147243
#### stand:
#     [,1]
#[1,]    1
#
#### Infos:
#     method  message                                        
#[1,] "optIC" "minimum asymptotic bias (lower case) solution"

checkIC(IC4)

#precision of centering:  1.111799e-15                      
#precision of Fisher consistency:                     
#            prob                     
#prob 2.14051e-13                     
#maximum deviation                      
#      2.14051e-13                      

Risks(IC4)

#$asBias
#$asBias$value
#[1] 0.1089328
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "TotalVarNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$asCov
#[1] 0.002889749
#
#$trAsCov
#$trAsCov$value
#[1] 0.002889749
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#[1] 0.003631397
#
#$asMSE$r
#[1] 0.25
#
#$asMSE$at
#An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.5 


plot(IC4)


#-------------------------------------------------------------------------------
## Hampel solution
#-------------------------------------------------------------------------------
(IC5 <- optIC(model=RobN1, risk=asHampel(bound=clip(IC1))))

#minimal bound:   0.05458835 
#An object of class “ContIC” 
#### name:        IC of contamination type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clip:       [1] 0.06143505
#### cent:       [1] 0.003717089
#### stand:
#            prob
#prob 0.003100752
#
#### Infos:
#     method  message                                                
#[1,] "optIC" "optimally robust IC for 'asHampel' with bound = 0.061"


checkIC(IC5)

#precision of centering:  4.289037e-14 
#precision of Fisher consistency:
#              prob
#prob -5.088488e-07
#maximum deviation 
#     5.088488e-07 

Risks(IC5)

#$asCov
#            prob
#prob 0.002157190
#
#$asBias
#$asBias$value
#[1] 0.06143505
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "ContNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.002157190
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#            prob
#prob 0.003100757
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5 



plot(IC5)

(IC6 <- optIC(model=RobN2, risk=asHampel(bound=Risks(IC2)$asBias$value), maxiter = 200))

#minimal bound:   0.1089328 
#maximum iterations reached!
# achieved precision:     5.583403e-07 
#An object of class “TotalVarIC” 
#### name:        IC of total variation type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5 
#
#### clipLo:     [1] -0.06032159
#### clipUp:     [1] 0.052032
#### stand:
#            prob
#prob 0.005585234
#
#### Infos:
#     method  message                                                
#[1,] "optIC" "optimally robust IC for 'asHampel' with bound = 0.112"


checkIC(IC6)

#precision of centering:  9.55774e-16 
#precision of Fisher consistency:
#              prob
#prob -4.638466e-08
#maximum deviation 
#     4.638466e-08 


Risks(IC6)

#$asCov
#            prob
#prob 0.008255479
#
#$asBias
#$asBias$value
#[1] 0.1123536
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "TotalVarNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.008255479
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#           prob
#prob 0.01141131
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.5 



plot(IC6)


#-------------------------------------------------------------------------------
## radius minimax IC
#-------------------------------------------------------------------------------
system.time(IC7 <- radiusMinimaxIC(L2Fam=N, neighbor=ContNeighborhood(),
                        risk=asMSE(), loRad=0.01, upRad=3.9))
#   user  system elapsed
#  33.26    0.02   33.64 

IC7

#An object of class “ContIC” 
#### name:        IC of contamination type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5814856 
#
#### clip:       [1] 0.05991006
#### cent:       [1] 0.004363098
#### stand:
#            prob
#prob 0.003424638
#
#### Infos:
#     method            message                                            
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0.01, 3.9]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.581"                    
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.177"    

checkIC(IC7)

#precision of centering:  5.901277e-12 
#precision of Fisher consistency:
#              prob
#prob -1.295138e-09
#maximum deviation 
#     1.295138e-09 

Risks(IC7)

#$asCov
#            prob
#prob 0.002211033
#
#$asBias
#$asBias$value
#[1] 0.05991006
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "ContNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.002211033
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#            prob
#prob 0.003424638
#
#$asMSE$r
#[1] 0.5814856
#
#$asMSE$at
#An object of class “ContNeighborhood” 
#type:    (uncond.) convex contamination neighborhood 
#radius:  0.5814856 


plot(IC7)

system.time(IC8 <- radiusMinimaxIC(L2Fam=N, neighbor=TotalVarNeighborhood(),
                        risk=asMSE(), loRad=0.01, upRad=1.8))
#   user  system elapsed
# 565.58    0.21  586.05

IC8

#An object of class “TotalVarIC” 
#### name:        IC of total variation type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.2944932 
#
#### clipLo:     [1] -0.06497353
#### clipUp:     [1] 0.05431373
#### stand:
#            prob
#prob 0.003432502
#
#### Infos:
#     method            message                                            
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0.01, 1.8]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.294"                    
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.168"   


checkIC(IC8)

#precision of centering:  6.400051e-12 
#precision of Fisher consistency:
#              prob
#prob -1.404645e-09
#maximum deviation 
#     1.404645e-09 


Risks(IC8)

#$asCov
#            prob
#prob 0.003987582
#
#$asBias
#$asBias$value
#[1] 0.1192873
#
#$asBias$biastype
#An object of class "symmetricBias"
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#$asBias$neighbortype
#[1] "TotalVarNeighborhood"
#attr(,"package")
#[1] "RobAStBase"
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.003987582
#
#$trAsCov$normtype
#An object of class "NormType"
#Slot "name":
#[1] "EuclideanNorm"
#
#Slot "fct":
#function (x) 
#{
#    if (is.vector(x)) 
#        return(abs(x))
#    else return(sqrt(colSums(x^2)))
#}
#<environment: namespace:distrMod>
#
#
#
#$asMSE
#$asMSE$value
#            prob
#prob 0.005221649
#
#$asMSE$r
#[1] 0.2944932
#
#$asMSE$at
#An object of class “TotalVarNeighborhood” 
#type:    (uncond.) total variation neighborhood 
#radius:  0.2944932 

plot(IC8)


#-------------------------------------------------------------------------------
## least favorable radius
#-------------------------------------------------------------------------------
## (may take quite some time!)
system.time(r.rho1 <- leastFavorableRadius(L2Fam=N, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))

#current radius:  0.3820278      inefficiency:    1.040429 
#current radius:  0.6180722      inefficiency:    1.044464 
#current radius:  0.7639556      inefficiency:    1.041844 
#current radius:  0.5931816      inefficiency:    1.044631 
#current radius:  0.5546088      inefficiency:    1.044691 
#current radius:  0.5644109      inefficiency:    1.044698 
#current radius:  0.5640279      inefficiency:    1.0447 
#current radius:  0.5599945      inefficiency:    1.044697 
#current radius:  0.5624873      inefficiency:    1.044701 
#current radius:  0.5631909      inefficiency:    1.044701 
#current radius:  0.5627771      inefficiency:    1.044701 
#current radius:  0.5625595      inefficiency:    1.044701 
#current radius:  0.5626002      inefficiency:    1.044701 
#current radius:  0.5625595      inefficiency:    1.044701 
#   user  system elapsed 
# 361.25    0.12  369.05 

## same as for binomial????
 
r.rho1

#$rho
#[1] 0.5
#
#$leastFavorableRadius
#[1] 0.5625595
#
#$`asMSE-inefficiency`
#[1] 1.044701

system.time(r.rho2 <- leastFavorableRadius(L2Fam=N, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=0.5))

#current radius:  0.3820278      inefficiency:    1.041727 
#current radius:  0.6180722      inefficiency:    1.027733 
#current radius:  0.2361444      inefficiency:    1.043317 
#current radius:  0.2660735      inefficiency:    1.044275 
#current radius:  0.2943633      inefficiency:    1.044409 
#current radius:  0.2852759      inefficiency:    1.04443 
#current radius:  0.2866889      inefficiency:    1.044456 
#current radius:  0.2893884      inefficiency:    1.044426 
#current radius:  0.2872589      inefficiency:    1.044427 
#current radius:  0.2862418      inefficiency:    1.044439 
#current radius:  0.2869066      inefficiency:    1.044442 
#current radius:  0.2865879      inefficiency:    1.044429 
#current radius:  0.2867296      inefficiency:    1.044453 
#current radius:  0.2866482      inefficiency:    1.044425 
#current radius:  0.2866889      inefficiency:    1.044456 
#    user  system elapsed 
# 4891.07    1.90 5063.44 
 
r.rho2

#$rho
#[1] 0.5
#
#$leastFavorableRadius
#[1] 0.2866889
#
#$`asMSE-inefficiency`
#[1] 1.044456


###############################################################################
## k-step (k >= 1) estimation
################################################################################

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05)
x <- rnbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.01)

### MLE:

(estML <- MLEstimator(x=x, NbinomFamily(size=25)))

#Evaluations of Maximum likelihood estimate:
#-------------------------------------------
#An object of class “Estimate”
#generated by call
#  MLEstimator(x = x, ParamFamily = NbinomFamily(size = 25))
#samplesize:   100
#estimate:
#      prob
#  0.135624871
# (0.002521857)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.0006359763
#Criterion:
#negative log-likelihood
#               1868.396

## 2. Kolmogorov(-Smirnov) minimum distance estimator

(est0 <- MDEstimator(x=x, NbinomFamily(size=25)))

#Evaluations of Minimum Kolmogorov distance estimate:
#----------------------------------------------------
#An object of class “Estimate”
#generated by call
#  MDEstimator(x = x, ParamFamily = NbinomFamily(size = 25))
#samplesize:   100
#estimate:
#     prob
#0.2471440
#fixed part of the parameter:
#size
#  25
#Criterion:
#Kolmogorov distance
#         0.05226461

### 3.1. one-step estimation: radius known

### ac) Define infinitesimal robust model
RobN3 <- InfRobModel(center=NbinomFamily(size=25, prob=estimate(est0)),
                      neighbor=ContNeighborhood(radius=0.5))

## bc) Compute optimally robust IC

IC9 <- optIC(model=RobN3, risk=asMSE())
checkIC(IC9)

#precision of centering:  5.93858e-07
#precision of Fisher consistency:
#             prob
#prob 8.014067e-05
#maximum deviation
#     8.014067e-05

(est1c <- oneStepEstimator(x, IC=IC9, start=est0))

#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  oneStepEstimator(x = x, IC = IC9, start = est0)
#samplesize:   100
#estimate:
#      prob
#  0.251879937
# (0.004632576)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002146076
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03064529
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249199096954169
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5
#
#### clip:       [1] 0.06129058
#### cent:       [1] 0.003713275
#### stand:
#           prob
#prob 0.00308521
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 1

est1c1 <- roptest(x, NbinomFamily(size = 25), eps = 0.05, initial.est = est0)
checkIC(pIC(est1c1))

#precision of centering:  3.376922e-18
#precision of Fisher consistency:
#            prob
#prob 7.04511e-05
#maximum deviation
#      7.04511e-05

est1c2 <- roptest(x, NbinomFamily(size = 25), eps = 0.05, distance = KolmogorovDist)
checkIC(pIC(est1c2))

#precision of centering:  3.376922e-18
#precision of Fisher consistency:
#            prob
#prob 7.04511e-05
#maximum deviation
#      7.04511e-05

est1c3 <- roptest(x, NbinomFamily(size = 25), eps = 0.025)
checkIC(pIC(est1c3))

#precision of centering:  8.748485e-11
#precision of Fisher consistency:
#             prob
#prob 8.488734e-05
#maximum deviation
#     8.488734e-05


estimate(est1c)
#     prob
#0.2518799
estimate(est1c1)
#     prob
#0.2518792
estimate(est1c2)
#     prob
#0.2518792
estimate(est1c3)
#     prob
#0.2502157
confint(est1c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2426583 0.2611016
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#oneStepEstimator(x = x, IC = IC9, start = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1c1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2426577 0.2611008
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.05, initial.est = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1c2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2426577 0.2611008
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.05, distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1c3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2413840 0.2590475
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
## av) Define infinitesimal robust model
RobN4 <- InfRobModel(center=NbinomFamily(size=25, prob=estimate(est0)),
                 neighbor=TotalVarNeighborhood(radius=0.25))
## bv) Compute optimally robust IC
IC10 <- optIC(model=RobN4, risk=asMSE())
checkIC(IC10)
#precision of centering:  6.484067e-07
#precision of Fisher consistency:
#             prob
#prob 7.233166e-05
#maximum deviation
#     7.233166e-05
## cv) Determine 1-step estimate
(est1v <- oneStepEstimator(x, IC=IC10, start=est0))
#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  oneStepEstimator(x = x, IC = IC10, start = est0)
#samplesize:   100
#estimate:
#      prob
#  0.251562909
# (0.005826841)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.003395207
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03055969
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249199096954169
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.25
#
#### clipLo:     [1] -0.06692041
#### clipUp:     [1] 0.05531835
#### stand:
#            prob
#prob 0.003063591
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 1

## instead of av)-cv) you can also use function roptest
est1v1 <- roptest(x, NbinomFamily(size = 25), eps = 0.025, initial.est = est0,
                   neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v1))
#precision of centering:  4.984589e-14
#precision of Fisher consistency:
#             prob
#prob 6.263991e-05
#maximum deviation
#     6.263991e-05

## you can also omit step 2
est1v2 <- roptest(x, NbinomFamily(size = 25), eps = 0.025,
                   neighbor = TotalVarNeighborhood(), distance = KolmogorovDist)
checkIC(pIC(est1v2))
#precision of centering:  4.984589e-14
#precision of Fisher consistency:
#             prob
#prob 6.263991e-05
#maximum deviation
#     6.263991e-05

## Using Cramer-von-Mises MD estimator (default)
est1v3 <- roptest(x, NbinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v3))
#precision of centering:  3.647907e-14
#precision of Fisher consistency:
#             prob
#prob 6.370115e-05
#maximum deviation
#     6.370115e-05

## comparison of estimates
estimate(est1v)
#     prob
#0.2515629
estimate(est1v1)
#     prob
#0.2515621
estimate(est1v2)
#     prob
#0.2515621
estimate(est1v3)
#     prob
#0.2516431

## confidence intervals
confint(est1v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2399644 0.2631614
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#oneStepEstimator(x = x, IC = IC10, start = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1v1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2399638 0.2631603
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    initial.est = est0, neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1v2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2399638 0.2631603
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    neighbor = TotalVarNeighborhood(), distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1v3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2400033 0.2632828
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 3.2. k-step estimation: radius known
IC9 <- optIC(model=RobN3, risk=asMSE())
(est2c <- kStepEstimator(x, IC=IC9, start=est0, steps = 3L))
#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC9, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.252059546
# (0.004676936)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002187373
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03093139
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.252050979065772
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5
#
#### clip:       [1] 0.06186279
#### cent:       [1] 0.003740501
#### stand:
#            prob
#prob 0.003144124
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3


est2c1 <- roptest(x, NbinomFamily(size = 25), eps = 0.05, initial.est = est0, steps = 3L)
checkIC(pIC(est2c1))
#precision of centering:  -2.077441e-18
#precision of Fisher consistency:
#             prob
#prob 6.645565e-05
#maximum deviation
#     6.645565e-05

est2c2 <- roptest(x, NbinomFamily(size = 25), eps = 0.05, steps = 3L,
                   distance = KolmogorovDist)
checkIC(pIC(est2c2))
#precision of centering:  -2.077441e-18
#precision of Fisher consistency:
#             prob
#prob 6.645565e-05
#maximum deviation
#     6.645565e-05

# Using Cramer-von-Mises MD estimator
est2c3 <- roptest(x, NbinomFamily(size = 25), eps = 0.05, steps = 3L)
checkIC(pIC(est2c3))
#precision of centering:  5.311002e-18
#precision of Fisher consistency:
#             prob
#prob 6.642177e-05
#maximum deviation
#     6.642177e-05

## comparison of estimates
estimate(est2c)
#     prob
#0.2520595
estimate(est2c1)
#     prob
#0.2520595
estimate(est2c2)
#     prob
#0.2520595
estimate(est2c3)
#     prob
#0.2520597

## confidence intervals
confint(est2c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2427483 0.2613708
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#kStepEstimator(x = x, IC = IC9, start = est0, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2c1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2427483 0.2613708
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.05, initial.est = est0,
#    steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2c2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2427483 0.2613708
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.05, steps = 3L,
#    distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2c3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2427483 0.2613712
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.05, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25

IC10 <- optIC(model=RobN4, risk=asMSE())
(est2v <- kStepEstimator(x, IC=IC10, start=est0, steps = 3L))
#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC10, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.251762581
# (0.005876303)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.003453094
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03081803
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.251747739764054
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.25
#
#### clipLo:     [1] -0.06750276
#### clipUp:     [1] 0.05576935
#### stand:
#           prob
#prob 0.00311586
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3


checkIC(pIC(est2v))
#precision of centering:  5.496728e-13
#precision of Fisher consistency:
#             prob
#prob 6.135889e-05
#maximum deviation
#     6.135889e-05

est2v1 <- roptest(x, NbinomFamily(size = 25), eps = 0.025, initial.est = est0,
                   steps = 3L, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v1))
#precision of centering:  5.49678e-13
#precision of Fisher consistency:
#            prob
#prob 6.13594e-05
#maximum deviation
#      6.13594e-05

est2v2 <- roptest(x, NbinomFamily(size = 25), eps = 0.025, steps = 3L,
                   distance = KolmogorovDist, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v2))
#precision of centering:  5.49678e-13
#precision of Fisher consistency:
#            prob
#prob 6.13594e-05
#maximum deviation
#      6.13594e-05

## Using Cramer-von-Mises MD estimator
est2v3 <- roptest(x, NbinomFamily(size = 25), eps = 0.025, steps = 3L,
                   neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v3))
#precision of centering:  5.493408e-13
#precision of Fisher consistency:
#             prob
#prob 6.130999e-05
#maximum deviation
#     6.130999e-05

## comparison of estimates
estimate(est2v)
#     prob
#0.2517626
estimate(est2v1)
#     prob
#0.2517626
estimate(est2v2)
#     prob
#0.2517626
estimate(est2v3)
#     prob
#0.2517631

## confidence intervals
confint(est2v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %   97.5 %
#prob 0.2400641 0.263461
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#kStepEstimator(x = x, IC = IC10, start = est0, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2v1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %   97.5 %
#prob 0.2400641 0.263461
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    initial.est = est0, neighbor = TotalVarNeighborhood(), steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2v2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %   97.5 %
#prob 0.2400641 0.263461
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    neighbor = TotalVarNeighborhood(), steps = 3L, distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2v3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2400644 0.2634618
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps = 0.025,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 4.1. one-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=NbinomFamily(size=25, prob=estimate(est0)),
                 neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est3c <- oneStepEstimator(x, IC=IC11, start=est0))
#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  oneStepEstimator(x = x, IC = IC11, start = est0)
#samplesize:   100
#estimate:
#      prob
#  0.25268280
# (0.00470825)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002216761
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03607630
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249199096954169
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.6077586
#
#### clip:       [1] 0.05935958
#### cent:       [1] 0.004542211
#### stand:
#            prob
#prob 0.003518260
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, Inf]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.608"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.189"
#steps:
#[1] 1


checkIC(pIC(est3c))
#precision of centering:  5.751481e-07
#precision of Fisher consistency:
#             prob
#prob 7.761578e-05
#maximum deviation
#     7.761578e-05


IC12 <- radiusMinimaxIC(L2Fam=NbinomFamily(size=25, prob=estimate(est0)),
                 neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est3v <- oneStepEstimator(x, IC=IC12, start=est0))
#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  oneStepEstimator(x = x, IC = IC12, start = est0)
#samplesize:   100
#estimate:
#      prob
#  0.252242480
# (0.006465472)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.004180232
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03648773
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249199096954169
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3088797
#
#### clipLo:     [1] -0.06424131
#### clipUp:     [1] 0.05388796
#### stand:
#            prob
#prob 0.003537267
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, Inf]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.309"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.183"
#steps:
#[1] 1
checkIC(pIC(est3v))
#precision of centering:  6.224275e-07
#precision of Fisher consistency:
#             prob
#prob 7.046135e-05
#maximum deviation
#     7.046135e-05

## maximum radius for given sample size n: sqrt(n)*0.5
(est3c1 <- roptest(x, NbinomFamily(size = 25), eps.lower= 0.001, eps.upper = 0.5))

#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate” 
#generated by call
#  roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.lower = 0.001, 
#    eps.upper = 0.5)
#samplesize:   100
#estimate:
#      prob    
#  0.246522069 
# (0.004640427)
#fixed part of the parameter:
#size 
#  25 
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002153356
#Infos:
#     method    message                                                  
#[1,] "roptest" "1-step estimate for Negative Binomial family"           
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03488933
#(partial) influence curve:
#An object of class “ContIC” 
#### name:        IC of contamination type 
#
#### L2-differentiable parametric family:         Negative Binomial family 
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.245546087232864
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5928425 
#
#### clip:       [1] 0.05885092
#### cent:       [1] 0.004367926
#### stand:
#            prob
#prob 0.003370622
#
#### Infos:
#     method            message                                          
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0.01, 5]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.593"                  
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.181"            
#steps:
#[1] 1


checkIC(pIC(est3c1))

#precision of centering:  -2.995982e-18 
#precision of Fisher consistency:
#             prob
#prob 6.959669e-05
#maximum deviation 
#     6.959669e-05 

(est3v1 <- roptest(x, NbinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood()))

#maximum iterations reached!
# achieved precision:     0.02142607
#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood())
#samplesize:   100
#estimate:
#      prob
#  0.25237905
# (0.00646375)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.004178006
#Infos:
#     method    message
#[1,] "roptest" "1-step estimate for Negative Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03641902
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.250253492473826
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3068304
#
#### clipLo:     [1] -0.06456431
#### clipUp:     [1] 0.05412999
#### stand:
#            prob
#prob 0.003544421
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, 5]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.307"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.181"
#steps:
#[1] 1

checkIC(pIC(est3v1))

#precision of centering:  5.876106e-18
#precision of Fisher consistency:
#            prob
#prob 6.18579e-05
#maximum deviation
#      6.18579e-05


## comparison of estimates
estimate(est3c)
#     prob
#0.2526828
estimate(est3v)
#     prob
#0.2522425
estimate(est3c1)
#     prob 
#0.2455461
estimate(est3v1)
#     prob 
#0.25237905

## confidence intervals
confint(est3c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2432849 0.2620807
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#oneStepEstimator(x = x, IC = IC11, start = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est3v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2393345 0.2651505
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#oneStepEstimator(x = x, IC = IC12, start = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est3c1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2372651 0.2557790
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est3v1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2394749 0.2652832
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 4.2. k-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=NbinomFamily(size=25, prob=estimate(est0)),
                 neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est4c <- kStepEstimator(x, IC=IC11, start=est0, steps = 3L))
#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC11, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.252744439
# (0.004763277)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002268880
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03651059
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.252743351251028
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.6077586
#
#### clip:       [1] 0.06007416
#### cent:       [1] 0.004611258
#### stand:
#            prob
#prob 0.003601903
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3


checkIC(pIC(est4c))
#precision of centering:  7.90373e-11
#precision of Fisher consistency:
#             prob
#prob 6.837836e-05
#maximum deviation
#     6.837836e-05

IC12 <- radiusMinimaxIC(L2Fam=NbinomFamily(size=25, prob=estimate(est0)),
                 neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC12, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.252684415
# (0.006539431)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.004276416
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Negative Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03691656
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.252646956249659
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3088797
#
#### clipLo:     [1] -0.06500297
#### clipUp:     [1] 0.05451465
#### stand:
#            prob
#prob 0.003619062
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4v))
#precision of centering:  1.975176e-14
#precision of Fisher consistency:
#             prob
#prob 6.190507e-05
#maximum deviation
#     6.190507e-05

# maximum radius for given sample size n: sqrt(n)*0.5
(est4c1 <- roptest(x, NbinomFamily(size = 25), eps.upper = 0.5, steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.lower = 0.001,
#    eps.upper = 0.5, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.252656191
# (0.004751824)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.002257983
#Infos:
#     method    message
#[1,] "roptest" "3-step estimate for Negative Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03572574
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.252655619320179
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5926806
#
#### clip:       [1] 0.06027824
#### cent:       [1] 0.004489273
#### stand:
#            prob
#prob 0.003534312
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4c1))

#precision of centering:  1.796719e-13
#precision of Fisher consistency:
#             prob
#prob 6.941785e-05
#maximum deviation
#     6.941785e-05

(est4v1 <- roptest(x, NbinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood(),
         steps = 3L))
#maximum iterations reached!
# achieved precision:     0.02142607
#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#samplesize:   100
#estimate:
#     prob
#  0.2526652
# (0.0065150)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.004244522
#Infos:
#     method    message
#[1,] "roptest" "3-step estimate for Negative Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.03670764
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Negative Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.252641193060971
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3068304
#
#### clipLo:     [1] -0.06507884
#### clipUp:     [1] 0.05455611
#### stand:
#            prob
#prob 0.003600885
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4v1))
#precision of centering:  1.138307e-14
#precision of Fisher consistency:
#             prob
#prob 6.200529e-05
#maximum deviation
#     6.200529e-05

## comparison of estimates
estimate(est4c)
#     prob
#0.2527444
estimate(est4v)
#     prob
#0.2526844
estimate(est4c1)
#     prob
#0.2526562
estimate(est4v1)
#     prob
#0.2526652

## confidence intervals
confint(est4c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2432347 0.2622542
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#kStepEstimator(x = x, IC = IC11, start = est0, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est4v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2396260 0.2657429
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#kStepEstimator(x = x, IC = IC12, start = est0, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est4c1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2431730 0.2621394
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.lower = 0.001,
#    eps.upper = 0.5, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est4v1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2396568 0.2656735
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = NbinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
