###############################################################################
## Example: Binomial Family
###############################################################################

require(ROptEst)
options("newDevice"=TRUE)

#-------------------------------------------------------------------------------
## Preparations
#-------------------------------------------------------------------------------
## generates Binomial Family with
## m = 25 and probability of success theta = 0.25
B <- BinomFamily(size = 25, prob = 0.25) 
B       # show B 

#An object of class "BinomFamily"
#### name:       Binomial family
#
#### distribution:       Distribution Object of Class: Binom
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
#[1] "The Binomial family is symmetric with respect to prob = 0.5;"
#[2] "i.e., d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)"

plot(B) # plot of Binom(size = 25, prob = 0.25) and L_2 derivative
checkL2deriv(B)

#precision of centering:  -3.103725e-15
#precision of Fisher information:
#              prob
#prob -2.842171e-14
#$maximum.deviation
#[1] 2.842171e-14

#-------------------------------------------------------------------------------
## classical optimal IC
#-------------------------------------------------------------------------------
IC0 <- optIC(model = B, risk = asCov())
IC0       # show IC

#An object of class “IC”
#### name:        Classical optimal influence curve for Binomial family
#### L2-differentiable parametric family:         Binomial family
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

#precision of centering:  -2.403827e-17
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.110223e-16

Risks(IC0)

#$asCov
#       prob
#prob 0.0075
#
#$trAsCov
#[1] 0.0075

#-------------------------------------------------------------------------------
## lower case radius
#-------------------------------------------------------------------------------
lowerCaseRadius(L2Fam = B, neighbor = ContNeighborhood(), risk = asMSE())

#lower case radius
#         1.197809

lowerCaseRadius(L2Fam = B, neighbor = TotalVarNeighborhood(), risk = asMSE())

#lower case radius
#         1.130615
         
#-------------------------------------------------------------------------------
## L_2 family + infinitesimal neighborhood
#-------------------------------------------------------------------------------
RobB1 <- InfRobModel(center = B, neighbor = ContNeighborhood(radius = 0.5))
RobB1     # show RobB1

#An object of class “InfRobModel”
####### center:  An object of class "BinomFamily"
#### name:       Binomial family
#
#### distribution:       Distribution Object of Class: Binom
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
#[1] "The Binomial family is symmetric with respect to prob = 0.5;"
#[2] "i.e., d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)"
#
####### neighborhood:    An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.5

(RobB2 <- InfRobModel(center = B, neighbor = TotalVarNeighborhood(radius = 0.5)))

#An object of class “InfRobModel”
####### center:  An object of class "BinomFamily"
#### name:       Binomial family
#
#### distribution:       Distribution Object of Class: Binom
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
#[1] "The Binomial family is symmetric with respect to prob = 0.5;"
#[2] "i.e., d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)"
#
####### neighborhood:    An object of class “TotalVarNeighborhood”
#type:    (uncond.) total variation neighborhood
#radius:  0.5

#-------------------------------------------------------------------------------
## OBRE solution
#-------------------------------------------------------------------------------

system.time(ICA <-  optIC(model=RobB1, risk=asAnscombe(),
            verbose=TRUE,lower=NULL,upper=10))


#-------------------------------------------------------------------------------
## MSE solution
#-------------------------------------------------------------------------------
system.time(IC1 <- optIC(model=RobB1, risk=asMSE()))

#   user  system elapsed
#   3.62    0.00    3.62

IC1

#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
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
#### clip:       [1] 0.1213661
#### cent:       [1] -0.003061948
#### stand:
#           prob
#prob 0.01222674
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"

checkIC(IC1)

#precision of centering:  -3.17283e-17
#precision of Fisher consistency:
#              prob
#prob -3.330669e-16
#maximum deviation
#     3.330669e-16

Risks(IC1)

#$asCov
#            prob
#prob 0.008544305
#
#$asBias
#$asBias$value
#[1] 0.1213661
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#prob 0.008544305
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.01222674
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
#[1] "Binom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) convex contamination neighborhood"
#
#$asBias$value
#[1] 0.1213661


getRiskIC(IC1, asMSE(), ContNeighborhood(radius = 0.5))

#$asMSE
#$asMSE$distribution
#[1] "Binom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) convex contamination neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5


(Cov1 <- getRiskIC(IC1, asCov()))

#$asCov
#$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$value
#$asCov$value[[1]]
#[1] 0.008544305
#
#$asCov$value$asCov
#$asCov$value$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$value$asCov$value
#            prob
#prob 0.008544305


(mse1 <- getRiskIC(IC1, asMSE(), TotalVarNeighborhood(radius = 0.5)))

#$asMSE
#$asMSE$distribution
#[1] "Binom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5
#
#$asMSE$value
#[1] 0.01222674

(bias1 <- getRiskIC(IC1, asBias(), TotalVarNeighborhood()))
#
#$asBias
#$asBias$distribution
#[1] "Binom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$value
#[1] 0.1213661

## only suboptimal -> ToDo-List
addRisk(IC1) <- list(Cov1, mse1, bias1)
Risks(IC1)

#$asCov
#$asCov[[1]]
#[1] 0.008544305
#
#$asCov$asCov
#$asCov$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$asCov$value
#            prob
#prob 0.008544305
#
#
#$asCov$asCov
#$asCov$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$asCov$value
#$asCov$asCov$value[[1]]
#[1] 0.008544305
#
#$asCov$asCov$value$asCov
#$asCov$asCov$value$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$asCov$value$asCov$value
#            prob
#prob 0.008544305
#
#
#$asBias
#$asBias$value
#[1] 0.1213661
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#[1] "Binom(25, 0.25)"
#
#$asBias$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$asBias$value
#[1] 0.1213661
#
#
#$asBias$asBias
#$asBias$asBias$distribution
#[1] "Binom(25, 0.25)"
#
#$asBias$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$asBias$value
#[1] 0.1213661
#
#
#
#$trAsCov
#$trAsCov$value
#            prob
#prob 0.008544305
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.01222674
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
#[1] "Binom(25, 0.25)"
#
#$asMSE$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$asMSE$radius
#[1] 0.5
#
#$asMSE$asMSE$value
#[1] 0.01222674
#
#
#$asMSE$asMSE
#$asMSE$asMSE$distribution
#[1] "Binom(25, 0.25)"
#
#$asMSE$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$asMSE$radius
#[1] 0.5
#
#$asMSE$asMSE$value
#[1] 0.01222674

plot(IC1)

system.time(IC2 <- optIC(model=RobB2, risk=asMSE()))

#   user  system elapsed
#   9.46    0.02    9.47

IC2

#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
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
#### clipLo:     [1] -0.1081481
#### clipUp:     [1] 0.1158377
#### stand:
#           prob
#prob 0.02208611
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"

checkIC(IC2)

#precision of centering:  -1.745301e-13
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.745301e-13

Risks(IC2)

#$asCov
#           prob
#prob 0.03342380
#
#$asBias
#$asBias$value
#[1] 0.2239858
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#           prob
#prob 0.03342380
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.04596621
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
#[1] "Binom(25, 0.25)"
#
#$asMSE$neighborhood
#[1] "(uncond.) total variation neighborhood with radius 0.5"
#
#$asMSE$radius
#[1] 0.5
#
#$asMSE$value
#[1] 0.04596621

getRiskIC(IC2, asBias(), TotalVarNeighborhood())

#$asBias
#$asBias$distribution
#[1] "Binom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) total variation neighborhood"
#
#$asBias$value
#[1] 0.2239858

getRiskIC(IC2, asBias(), ContNeighborhood())

#$asBias
#$asBias$distribution
#[1] "Binom(25, 0.25)"
#
#$asBias$neighborhood
#[1] "(uncond.) convex contamination neighborhood"
#
#$asBias$value
#[1] 0.2239858

Cov2 <- getRiskIC(IC2, asCov())
addRisk(IC2) <- Cov2
Risks(IC2)

#$asCov
#$asCov[[1]]
#[1] 0.03342380
#
#$asCov$asCov
#$asCov$asCov$distribution
#[1] "Binom(25, 0.25)"
#
#$asCov$asCov$value
#           prob
#prob 0.03342380
#
#$asBias
#$asBias$value
#[1] 0.2239858
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#           prob
#prob 0.03342380
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.04596621
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
(IC3 <- optIC(model=RobB1, risk=asBias()))

#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
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
#### clip:       [1] 0.1098910
#### cent:       [1] -1.333333
#### stand:
#     [,1]
#[1,]    1
#### lowerCase:  [1] -0.3316026
#
#### Infos:
#     method  message
#[1,] "optIC" "minimum asymptotic bias (lower case) solution"

checkIC(IC3)

#precision of centering:  -1.714470e-17
#precision of Fisher consistency:
#             prob
#prob 2.220446e-16
#maximum deviation
#     2.220446e-16

Risks(IC3)

#$asBias
#$asBias$value
#[1] 0.1098910
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#[1] 0.01011106
#
#$trAsCov
#$trAsCov$value
#[1] 0.01011106
#
#$trAsCov$normtype
#An object of class “NormType”
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
#[1] 0.01086581
#
#$asMSE$r
#[1] 0.25
#
#$asMSE$at
#An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.5

plot(IC3)

(IC4 <- optIC(model=RobB2, risk=asBias()))

#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
### param:      An object of class "ParamFamParameter"
##name:   probability of success
#prob:   0.25
##fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5
#
#### clipLo:     [1] -0.094766
#### clipUp:     [1] 0.1211501
#### stand:
#     [,1]
#[1,]    1
#
#### Infos:
#     method  message
#[1,] "optIC" "minimum asymptotic bias (lower case) solution"

checkIC(IC4)

#precision of centering:  -1.076010e-17
#precision of Fisher consistency:
#              prob
#prob -3.330669e-16
#maximum deviation
#     3.330669e-16

Risks(IC4)

#$asBias
#$asBias$value
#[1] 0.2159161
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#[1] 0.01148091
#
#$trAsCov
#$trAsCov$value
#[1] 0.01148091
#
#$trAsCov$normtype
#An object of class “NormType”
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
#[1] 0.01439465
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
(IC5 <- optIC(model=RobB1, risk=asHampel(bound=clip(IC1))))

#minimal bound:   0.1098910
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
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
#### clip:       [1] 0.1213661
#### cent:       [1] -0.003061961
#### stand:
#           prob
#prob 0.01222672
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for 'asHampel' with bound = 0.121"

checkIC(IC5)

#precision of centering:  -3.524202e-17
#precision of Fisher consistency:
#             prob
#prob -3.56339e-07
#maximum deviation
#      3.56339e-07

Risks(IC5)

#$asCov
#            prob
#prob 0.008544296
#
#$asBias
#$asBias$value
#[1] 0.1213661
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#prob 0.008544296
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.01222673
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.5


plot(IC5)

(IC6 <- optIC(model=RobB2, risk=asHampel(bound=Risks(IC2)$asBias$value), maxiter = 200))

#minimal bound:   0.2159161
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
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
#### clipLo:     [1] -0.1081481
#### clipUp:     [1] 0.1158378
#### stand:
#           prob
#prob 0.02208608
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for 'asHampel' with bound = 0.224"

checkIC(IC6)

#precision of centering:  -1.745225e-13
#precision of Fisher consistency:
#              prob
#prob -1.166564e-07
#maximum deviation
#     1.166564e-07

Risks(IC6)

#$asCov
#           prob
#prob 0.03342372
#
#$asBias
#$asBias$value
#[1] 0.2239858
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#           prob
#prob 0.03342372
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.04596613
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
system.time(IC7 <- radiusMinimaxIC(L2Fam=B, neighbor=ContNeighborhood(),
                        risk=asMSE(), loRad=0, upRad=1))
#   user  system elapsed
#  39.39    0.02   39.45

IC7

#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3913008
#
#### clip:       [1] 0.1269559
#### cent:       [1] -0.004305063
#### stand:
#           prob
#prob 0.01074017
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, 1]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.391"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.103"

checkIC(IC7)

#precision of centering:  -3.056847e-17
#precision of Fisher consistency:
#              prob
#prob -3.330669e-16
#maximum deviation
#     3.330669e-16

Risks(IC7)

#$asCov
#            prob
#prob 0.008272273
#
#$asBias
#$asBias$value
#[1] 0.1269559
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#prob 0.008272273
#
#$trAsCov$normtype
#An object of class “NormType”
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
#prob 0.01074017
#
#$asMSE$r
#[1] 0.3913008
#
#$asMSE$at
#An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.3913008

plot(IC7)

system.time(IC8 <- radiusMinimaxIC(L2Fam=B, neighbor=TotalVarNeighborhood(),
                        risk=asMSE(), loRad=0, upRad=1))
#   user  system elapsed
# 163.20    0.12  168.21

IC8

#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.25
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.2661058
#
#### clipLo:     [1] -0.1157700
#### clipUp:     [1] 0.1248992
#### stand:
#           prob
#prob 0.01269835
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, 1]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.266"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.146"

checkIC(IC8)

#precision of centering:  6.021289e-10
#precision of Fisher consistency:
#              prob
#prob -2.220446e-16
#maximum deviation
#     6.021289e-10

Risks(IC8)

#$asCov
#           prob
#prob 0.01546324
#
#$asBias
#$asBias$value
#[1] 0.2406692
#
#$asBias$biastype
#An object of class “symmetricBias”
#Slot "name":
#[1] "symmetric Bias"
#
#
#$asBias$normtype
#An object of class “NormType”
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
#           prob
#prob 0.01546324
#
#$trAsCov$normtype
#An object of class “NormType”
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
#          prob
#prob 0.0195648
#
#$asMSE$r
#[1] 0.2661058
#
#$asMSE$at
#An object of class “TotalVarNeighborhood”
#type:    (uncond.) total variation neighborhood
#radius:  0.2661058

plot(IC8)


#-------------------------------------------------------------------------------
## least favorable radius
#-------------------------------------------------------------------------------
## (may take quite some time!)
system.time(r.rho1 <- leastFavorableRadius(L2Fam=B, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))

#current radius:  0.3820278      inefficiency:    1.038184
#current radius:  0.6180722      inefficiency:    1.044504
#current radius:  0.7639556      inefficiency:    1.042432
#current radius:  0.6248193      inefficiency:    1.044545
#current radius:  0.6423138      inefficiency:    1.044575
#current radius:  0.6887768      inefficiency:    1.044172
#current radius:  0.6381993      inefficiency:    1.044581
#current radius:  0.6371031      inefficiency:    1.044575
#current radius:  0.6397095      inefficiency:    1.044575
#current radius:  0.6387762      inefficiency:    1.044584
#current radius:  0.6387355      inefficiency:    1.044584
#current radius:  0.6391327      inefficiency:    1.044573
#current radius:  0.6389123      inefficiency:    1.044572
#current radius:  0.6388282      inefficiency:    1.044571
#current radius:  0.6387762      inefficiency:    1.044584
#   user  system elapsed
# 437.47    0.09  438.58

r.rho1

#$rho
#[1] 0.5
#
#$leastFavorableRadius
#[1] 0.6387762
#
#$`asMSE-inefficiency`
#[1] 1.044584

system.time(r.rho2 <- leastFavorableRadius(L2Fam=B, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=0.5))

#current radius:  0.3820278      inefficiency:    1.040750
#current radius:  0.6180722      inefficiency:    1.028895
#current radius:  0.2361444      inefficiency:    1.043035
#current radius:  0.2225705      inefficiency:    1.042438
#current radius:  0.2881511      inefficiency:    1.043182
#current radius:  0.3240088      inefficiency:    1.043234
#current radius:  0.3524144      inefficiency:    1.042468
#current radius:  0.3103124      inefficiency:    1.043323
#current radius:  0.3081042      inefficiency:    1.043320
#current radius:  0.3106103      inefficiency:    1.043323
#current radius:  0.3105373      inefficiency:    1.043323
#current radius:  0.3104966      inefficiency:    1.043323
#current radius:  0.3105373      inefficiency:    1.043323
#   user  system elapsed
#2211.92    1.13 2235.70

r.rho2

#$rho
#[1] 0.5
#
#$leastFavorableRadius
#[1] 0.3105373
#
#$`asMSE-inefficiency`
#[1] 1.043323


###############################################################################
## k-step (k >= 1) estimation
###############################################################################

## one-step estimation
## 1. generate a contaminated sample
ind <- rbinom(100, size=1, prob=0.05)
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.75)

## 2. Kolmogorov(-Smirnov) minimum distance estimator
(est0 <- MDEstimator(x=x, BinomFamily(size=25)))

#Evaluations of Minimum Kolmogorov distance estimate:
#----------------------------------------------------
#An object of class “Estimate”
#generated by call
#  MDEstimator(x = x, ParamFamily = BinomFamily(size = 25))
#samplesize:   100
#estimate:
#     prob
#0.2494779
#fixed part of the parameter:
#size
#  25
#Criterion:
#Kolmogorov distance
#         0.05944897

## 3.1. one-step estimation: radius known
## ac) Define infinitesimal robust model
RobB3 <- InfRobModel(center=BinomFamily(size=25, prob=estimate(est0)),
                     neighbor=ContNeighborhood(radius=0.5))
## bc) Compute optimally robust IC
IC9 <- optIC(model=RobB3, risk=asMSE())
checkIC(IC9)

#precision of centering:  -2.971365e-17
#precision of Fisher consistency:
#              prob
#prob -2.220446e-16
#maximum deviation
#     2.220446e-16

## cc) Determine 1-step estimate
(est1c <- oneStepEstimator(x, IC=IC9, start=est0))

#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  oneStepEstimator(x = x, IC = IC9, start = est0)
#samplesize:   100
#estimate:
#      prob
#  0.24253514
# (0.00923989)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.008537558
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.06058743
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249477874158289
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5
#
#### clip:       [1] 0.1211749
#### cent:       [1] -0.003012083
#### stand:
#           prob
#prob 0.01220839
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 1

## instead of ac)-cc) you can also use function roptest
est1c1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, initial.est = est0)
checkIC(pIC(est1c1))

#precision of centering:  4.569788e-18
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.440892e-16

## you can also omit step 2
est1c2 <- roptest(x, BinomFamily(size = 25), eps = 0.05, distance = KolmogorovDist)
checkIC(pIC(est1c2))

#precision of centering:  4.569788e-18
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.440892e-16

## Using Cramer-von-Mises MD estimator (default)
est1c3 <- roptest(x, BinomFamily(size = 25), eps = 0.025)
checkIC(pIC(est1c3))

#precision of centering:  -7.770327e-10
#precision of Fisher consistency:
#             prob
#prob 2.220446e-16
#maximum deviation
#     7.770327e-10

## comparison of estimates
estimate(est1c)
#     prob
#0.2425351
estimate(est1c1)
#     prob
#0.2425351
estimate(est1c2)
#     prob
#0.2425351
estimate(est1c3)
#     prob
#0.2427226

## confidence intervals
confint(est1c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2238655 0.2612048
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
#prob 0.2238655 0.2612048
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.05, initial.est = est0)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1c2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2238655 0.2612048
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.05, distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1c3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2251586 0.2602866
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025)
#Fixed part of the parameter at which estimate was produced:
#size
#  25

## av) Define infinitesimal robust model
RobB4 <- InfRobModel(center=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(radius=0.25))
## bv) Compute optimally robust IC
IC10 <- optIC(model=RobB4, risk=asMSE())
checkIC(IC10)
## cv) Determine 1-step estimate
(est1v <- oneStepEstimator(x, IC=IC10, start=est0))

## instead of av)-cv) you can also use function roptest
est1v1 <- roptest(x, BinomFamily(size = 25), eps = 0.025, initial.est = est0,
                  neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v1))
#precision of centering:  -9.843168e-15
#precision of Fisher consistency:
#              prob
#prob -2.220446e-16
#maximum deviation
#
     9.843168e-15
## you can also omit step 2
est1v2 <- roptest(x, BinomFamily(size = 25), eps = 0.025,
                  neighbor = TotalVarNeighborhood(), distance = KolmogorovDist)
checkIC(pIC(est1v2))
#precision of centering:  -9.843168e-15
#precision of Fisher consistency:
#              prob
#prob -2.220446e-16
#maximum deviation
#     9.843168e-15

## Using Cramer-von-Mises MD estimator (default)
est1v3 <- roptest(x, BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood())
checkIC(pIC(est1v3))
#precision of centering:  6.67169e-18
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.110223e-16

## comparison of estimates
estimate(est1v)
#     prob
#0.2429921
estimate(est1v1)
#     prob
#0.2429921
estimate(est1v2)
#     prob
#0.2429921
estimate(est1v3)
#     prob
#0.2424179

## confidence intervals
confint(est1v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2186101 0.2673741
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
#prob 0.2186101 0.2673741
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, initial.est = est0,
#    neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1v2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2186101 0.2673741
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood(),
#    distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est1v3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2182528 0.2665830
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 3.2. k-step estimation: radius known
IC9 <- optIC(model=RobB3, risk=asMSE())
(est2c <- kStepEstimator(x, IC=IC9, start=est0, steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC9, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.241900493
# (0.009141218)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.008356187
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.05992703
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.241953383964985
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5
#
#### clip:       [1] 0.1198541
#### cent:       [1] -0.004676326
#### stand:
#           prob
#prob 0.01194744
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

est2c1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, initial.est = est0, steps = 3L)
checkIC(pIC(est2c1))
#precision of centering:  -4.788156e-12
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.788156e-12
est2c2 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3L,
                  distance = KolmogorovDist)
checkIC(pIC(est2c2))
#precision of centering:  -4.788156e-12
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.788156e-12
# Using Cramer-von-Mises MD estimator
est2c3 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3L)
checkIC(pIC(est2c3))
#precision of centering:  -6.068637e-12
#precision of Fisher consistency:
#             prob
#prob 4.440892e-16
#maximum deviation
#     6.068637e-12

## comparison of estimates
estimate(est2c)
#     prob
#0.2419005
estimate(est2c1)
#     prob
#0.2419005
estimate(est2c2)
#     prob
#0.2419005
estimate(est2c3)
#     prob
#0.2418961

## confidence intervals
confint(est2c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2234362 0.2603648
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
#prob 0.2234362 0.2603648
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.05, initial.est = est0,
#    steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2c2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2234362 0.2603648
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.05, steps = 3L,
#    distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2c3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2234331 0.2603591
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.05, steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
  
IC10 <- optIC(model=RobB4, risk=asMSE())
(est2v <- kStepEstimator(x, IC=IC10, start=est0, steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC10, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.24237558
# (0.01195027)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.01428089
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.05982794
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.242410409425593
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.25
#
#### clipLo:     [1] -0.1148244
#### clipUp:     [1] 0.1244874
#### stand:
#           prob
#prob 0.01194066
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est2v))
#precision of centering:  -2.680039e-18
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.110223e-16

est2v1 <- roptest(x, BinomFamily(size = 25), eps = 0.025, initial.est = est0, 
                  steps = 3L, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v1))
#precision of centering:  7.299035e-18
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.440892e-16
     
est2v2 <- roptest(x, BinomFamily(size = 25), eps = 0.025, steps = 3L,
                  distance = KolmogorovDist, neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v2))
#precision of centering:  7.299035e-18
#precision of Fisher consistency:
#              prob
#prob -4.440892e-16
#maximum deviation
#     4.440892e-16

## Using Cramer-von-Mises MD estimator
est2v3 <- roptest(x, BinomFamily(size = 25), eps = 0.025, steps = 3L, 
                  neighbor = TotalVarNeighborhood())
checkIC(pIC(est2v3))
#precision of centering:  3.543879e-18
#precision of Fisher consistency:
#              prob
#prob -3.330669e-16
#maximum deviation
#     3.330669e-16

## comparison of estimates
estimate(est2v)
#     prob
#0.2423756
estimate(est2v1)
#     prob
#0.2423756
estimate(est2v2)
#     prob
#0.2423756
estimate(est2v3)
#     prob
#0.2423735

## confidence intervals
confint(est2v, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2182385 0.2665126
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
#         2.5 %    97.5 %
#prob 0.2182385 0.2665126
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, initial.est = est0,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2v2, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2182385 0.2665126
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood(),
#    steps = 3L, distance = KolmogorovDist)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est2v3, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2182378 0.2665092
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps = 0.025, neighbor = TotalVarNeighborhood(),
#    steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 4.1. one-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
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
#  0.242379282
# (0.009339235)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.00872213
#Infos:
#     method
#[1,] "oneStepEstimator"
#[2,] "oneStepEstimator"
#     message
#[1,] "1-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07121993
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.249477874158289
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.6006922
#
#### clip:       [1] 0.1185631
#### cent:       [1] -0.00453307
#### stand:
#           prob
#prob 0.01379441
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, Inf]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.601"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.165"
#steps:
#[1] 1

checkIC(pIC(est3c))

#precision of centering:  9.537084e-18
#precision of Fisher consistency:
#              prob
#prob -2.220446e-16
#maximum deviation
#     2.220446e-16

IC12 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est3v <- oneStepEstimator(x, IC=IC12, start=est0))
checkIC(pIC(est3v))

## maximum radius for given sample size n: sqrt(n)*0.5
(est3c1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5))
#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5)
#samplesize:   100
#estimate:
#      prob
#  0.241842958
# (0.009234894)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.008528328
#Infos:
#     method    message
#[1,] "roptest" "1-step estimate for Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07044238
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.24311883705201
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5967953
#
#### clip:       [1] 0.1180344
#### cent:       [1] -0.005169814
#### stand:
#           prob
#prob 0.01349046
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, 5]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.597"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.159"
#steps:
#[1] 1

checkIC(pIC(est3c1))

#precision of centering:  -1.335405e-18
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.110223e-16

(est3v1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood()))

#Evaluations of 1-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood())
#samplesize:   100
#estimate:
#      prob
#  0.24231510
# (0.01333293)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.01777669
#Infos:
#     method    message
#[1,] "roptest" "1-step estimate for Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07270451
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.24311883705201
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3104128
#
#### clipLo:     [1] -0.1113605
#### clipUp:     [1] 0.1228583
#### stand:
#           prob
#prob 0.01386519
#
#### Infos:
#     method            message
#[1,] "radiusMinimaxIC" "radius minimax IC for radius interval [0, 5]"
#[2,] "radiusMinimaxIC" "least favorable radius: 0.31"
#[3,] "radiusMinimaxIC" "maximum ‘asMSE’-inefficiency: 1.166"
#steps:
#[1] 1

checkIC(pIC(est3v1))

#precision of centering:  -1.179940e-11
#precision of Fisher consistency:
#     prob
#prob    0
#maximum deviation
#     1.179940e-11


## comparison of estimates
estimate(est3c)
#     prob
#0.2423793
estimate(est3v)
#     prob
#0.2428346
estimate(est3c1)
#     prob
#0.2418430
estimate(est3v1)
#     prob
#0.2423151

## confidence intervals
confint(est3c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %   97.5 %
#prob 0.2234096 0.261349
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
#prob 0.2154232 0.2702460
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
#prob 0.2230924 0.2605935
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
#prob 0.2152137 0.2694165
#Type of estimator: 1-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood())
#Fixed part of the parameter at which estimate was produced:
#size
#  25


## 4.2. k-step estimation: radius interval
IC11 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
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
#  0.24171362
# (0.00922172)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.008504013
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07071585
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.241769273363845
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.6006922
#
#### clip:       [1] 0.1177239
#### cent:       [1] -0.005337597
#### stand:
#           prob
#prob 0.01350474
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4c))
#precision of centering:  -1.717330e-17
#precision of Fisher consistency:
#             prob
#prob 2.220446e-16
#maximum deviation
#     2.220446e-16
     
IC12 <- radiusMinimaxIC(L2Fam=BinomFamily(size=25, prob=estimate(est0)),
                neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(est4v <- kStepEstimator(x, IC=IC12, start=est0, steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  kStepEstimator(x = x, IC = IC12, start = est0, steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.24224203
# (0.01335018)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.01782272
#Infos:
#     method
#[1,] "kStepEstimator"
#[2,] "kStepEstimator"
#     message
#[1,] "3-step estimate for Binomial family"
#[2,] "computation of IC, trafo, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07293007
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.242288609339931
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3116990
#
#### clipLo:     [1] -0.1111696
#### clipUp:     [1] 0.1228064
#### stand:
#           prob
#prob 0.01387656
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4v))
#precision of centering:  -1.454523e-11
#precision of Fisher consistency:
#             prob
#prob 4.440892e-16
#maximum deviation
#     1.454523e-11

## maximum radius for given sample size n: sqrt(n)*0.5
(est4c1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.241716814
# (0.009217234)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.00849574
#Infos:
#     method    message
#[1,] "roptest" "3-step estimate for Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07030448
#(partial) influence curve:
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.241727413448719
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.5967953
#
#### clip:       [1] 0.1178033
#### cent:       [1] -0.005316906
#### stand:
#           prob
#prob 0.01343846
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4c1))
#precision of centering:  -1.230939e-17
#precision of Fisher consistency:
#              prob
#prob -3.330669e-16
#maximum deviation
#     3.330669e-16

(est4v1 <- roptest(x, BinomFamily(size = 25), eps.upper = 0.5, neighbor = TotalVarNeighborhood(),
        steps = 3L))

#Evaluations of 3-step estimate:
#-------------------------------
#An object of class “Estimate”
#generated by call
#  roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#samplesize:   100
#estimate:
#      prob
#  0.24224135
# (0.01331745)
#fixed part of the parameter:
#size
#  25
#asymptotic (co)variance (multiplied with samplesize):
#[1] 0.01773545
#Infos:
#     method    message
#[1,] "roptest" "3-step estimate for Binomial family"
#[2,] "roptest" "computation of IC, asvar and asbias via useLast = FALSE"
#asymptotic bias:
#[1] 0.07266429
#(partial) influence curve:
#An object of class “TotalVarIC”
#### name:        IC of total variation type
#
#### L2-differentiable parametric family:         Binomial family
#### param:      An object of class "ParamFamParameter"
#name:   probability of success
#prob:   0.242247135068009
#fixed part of param.:
#        size:   25
#trafo:
#     prob
#prob    1
#
#### neighborhood radius:         0.3104128
#
#### clipLo:     [1] -0.1112452
#### clipUp:     [1] 0.122844
#### stand:
#           prob
#prob 0.01383111
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
#steps:
#[1] 3

checkIC(pIC(est4v1))
#precision of centering:  -1.224149e-11
#precision of Fisher consistency:
#              prob
#prob -1.110223e-16
#maximum deviation
#     1.224149e-11

## comparison of estimates
estimate(est4c)
#     prob
#0.2417136
estimate(est4v)
#     prob
#0.2422420
estimate(est4c1)
#     prob
#0.2417168
estimate(est4v1)
#     prob
#0.2422414

## confidence intervals
confint(est4c, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %  97.5 %
#prob 0.2229873 0.26044
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
#prob 0.2151025 0.2693815
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
#prob 0.2230034 0.2604303
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25
confint(est4v1, symmetricBias())
#A[n] asymptotic (LAN-based), uniform (bias-aware)
# confidence interval:
#for symmetric Bias
#         2.5 %    97.5 %
#prob 0.2151719 0.2693108
#Type of estimator: 3-step estimate
#samplesize:   100
#Call by which estimate was produced:
#roptest(x = x, L2Fam = BinomFamily(size = 25), eps.upper = 0.5,
#    neighbor = TotalVarNeighborhood(), steps = 3L)
#Fixed part of the parameter at which estimate was produced:
#size
#  25