###############################################################################
## Example: Exponential Scale Family
###############################################################################
require(ROptEst)
options("newDevice"=TRUE)

## generates Exponential Scale Family with scale = 0.5 (rate = 2)
E1 <- ExpScaleFamily(scale = 0.5) 
E1        # show E1
#An object of class "ExpScaleFamily"
#### name:       Exponential scale family
#
#### distribution:       Distribution Object of Class: Exp
# rate: 2
#
#### param:      An object of class "ParamFamParameter"
#name:   scale
#scale:  0.5
#fixed part of param.:
#        :       0
#trafo:
#      scale
#scale     1
#
#### props:
#[1] "The Exponential scale family is invariant under"
#[2] "the group of transformations 'g(y) = scale*y'"
#[3] "with scale parameter 'scale'"
plot(E1)  # plot of Exp(rate = 1) and L_2 derivative
checkL2deriv(E1)
#precision of centering:  -3.023619e-06
#precision of Fisher information:
#              scale
#scale -0.0001047172
#$maximum.deviation
#[1] 0.0001047172

## classical optimal IC
E1.IC0 <- optIC(model = E1, risk = asCov())
E1.IC0       # show IC
#An object of class “IC”
#### name:        Classical optimal influence curve for Exponential scale family
#### L2-differentiable parametric family:         Exponential scale family
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
checkIC(E1.IC0)
#precision of centering:  -7.559048e-07
#precision of Fisher consistency:
#             scale
#scale -2.61793e-05
#maximum deviation
#      2.61793e-05
Risks(E1.IC0)
#$asCov
#      scale
#scale  0.25
#
#$trAsCov
#[1] 0.25
plot(E1.IC0) # plot IC

## L_2 family + infinitesimal neighborhood
E1.Rob1 <- InfRobModel(center = E1, neighbor = ContNeighborhood(radius = 0.5))
E1.Rob1     # show E1.Rob1
#An object of class “InfRobModel”
####### center:  An object of class "ExpScaleFamily"
#### name:       Exponential scale family
#
#### distribution:       Distribution Object of Class: Exp
# rate: 2
#
#### param:      An object of class "ParamFamParameter"
#name:   scale
#scale:  0.5
#fixed part of param.:
#        :       0
#trafo:
#      scale
#scale     1
#
#### props:
#[1] "The Exponential scale family is invariant under"
#[2] "the group of transformations 'g(y) = scale*y'"
#[3] "with scale parameter 'scale'"
#
####### neighborhood:    An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.5
(E1.Rob2 <- InfRobModel(center = E1, neighbor = TotalVarNeighborhood(radius = 0.5)))
#An object of class “InfRobModel”
####### center:  An object of class "ExpScaleFamily"
#### name:       Exponential scale family
#
#### distribution:       Distribution Object of Class: Exp
# rate: 2
#
#### param:      An object of class "ParamFamParameter"
#name:   scale
#scale:  0.5
#fixed part of param.:
#        :       0
#trafo:
#      scale
#scale     1
#
#### props:
#[1] "The Exponential scale family is invariant under"
#[2] "the group of transformations 'g(y) = scale*y'"
#[3] "with scale parameter 'scale'"
#
####### neighborhood:    An object of class “TotalVarNeighborhood”
#type:    (uncond.) total variation neighborhood
#radius:  0.5

## MSE solution
(E1.IC1 <- optIC(model=E1.Rob1, risk=asMSE()))
#An object of class “ContIC”
#### name:        IC of contamination type
#
#### L2-differentiable parametric family:         Exponential scale family
#### param:      An object of class "ParamFamParameter"
#name:   scale
#scale:  0.5
#fixed part of param.:
#        :       0
#trafo:
#      scale
#scale     1
#
### neighborhood radius:         0.5
#
#### clip:       [1] 0.8449229
#### cent:       [1] -0.2112307
#### stand:
#          scale
#scale 0.5249744
#
#### Infos:
#     method  message
#[1,] "optIC" "optimally robust IC for ‘asMSE’"
checkIC(E1.IC1)
#precision of centering:  -2.924022e-05
#precision of Fisher consistency:
#              scale
#scale -1.068126e-05
#maximum deviation
#     2.924022e-05
Risks(E1.IC1)
#$asCov
#          scale
#scale 0.3465008
#
#$asBias
#$asBias$value
#[1] 0.8449229
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
#          scale
#scale 0.3465008
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
#          scale
#scale 0.5249744
#
#$asMSE$r
#[1] 0.5
#
#$asMSE$at
#An object of class “ContNeighborhood”
#type:    (uncond.) convex contamination neighborhood
#radius:  0.5

plot(E1.IC1)
(E1.IC2 <- optIC(model=E1.Rob2, risk=asMSE()))
checkIC(E1.IC2)
Risks(E1.IC2)
plot(E1.IC2)

## lower case solutions
(E1.IC3 <- optIC(model=E1.Rob1, risk=asBias()))
checkIC(E1.IC3)
Risks(E1.IC3)
plot(E1.IC3)
(E1.IC4 <- optIC(model=E1.Rob2, risk=asBias()))
checkIC(E1.IC4)
Risks(E1.IC4)
plot(E1.IC4)

## Hampel solution
(E1.IC5 <- optIC(model=E1.Rob1, risk=asHampel(bound=clip(E1.IC1))))
checkIC(E1.IC5)
Risks(E1.IC5)
plot(E1.IC5)
(E1.IC6 <- optIC(model=E1.Rob2, risk=asHampel(bound=Risks(E1.IC2)$asBias$value), maxiter = 200))
checkIC(E1.IC6)
Risks(E1.IC6)
plot(E1.IC6)

## radius minimax IC
(E1.IC7 <- radiusMinimaxIC(L2Fam=E1, neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=0.5))
checkIC(E1.IC7)
Risks(E1.IC7)
  (E1.IC8 <- radiusMinimaxIC(L2Fam=E1, neighbor=TotalVarNeighborhood(), risk=asMSE(), loRad=0, upRad=0.5))
checkIC(E1.IC8)
Risks(E1.IC8)

## least favorable radius
## (may take quite some time!)
(E1.r.rho1 <- leastFavorableRadius(L2Fam=E1, neighbor=ContNeighborhood(),
                    risk=asMSE(), rho=0.5))
(E1.r.rho2 <- leastFavorableRadius(L2Fam=E1, neighbor=TotalVarNeighborhood(),
                    risk=asMSE(), rho=1/3))

## one-step estimation
## We use contamination neighborhoods but could also use total variation 
## neighborhoods
## 1. generate a contaminated sample
ind <- rbinom(1e2, size=1, prob=0.05) 
E1.x <- rexp(1e2, rate=(1-ind)*2+ind*10)

## 2.1 Kolmogorov(-Smirnov) minimum distance estimator
(E1.est01 <- MDEstimator(x=E1.x, ExpScaleFamily()))

## 2.2 Cramer-von-Mises minimum distance estimator
(E1.est02 <- MDEstimator(x=E1.x, ExpScaleFamily(), distance = CvMDist))

## 3. k-step estimation: radius known
E1.Rob31 <- InfRobModel(center=ExpScaleFamily(scale=estimate(E1.est01)),
                neighbor=ContNeighborhood(radius=0.5))
E1.IC9 <- optIC(model=E1.Rob31, risk=asMSE())
(E1.est11 <- oneStepEstimator(E1.x, IC=E1.IC9, start=E1.est01))
(E1.est12 <- kStepEstimator(E1.x, IC=E1.IC9, start=E1.est01, steps = 3))

## its simpler to use function roptest
(E1.est13 <- roptest(E1.x, ExpScaleFamily(), eps = 0.05, distance = KolmogorovDist,
                     steps = 3))

E1.Rob32 <- InfRobModel(center=ExpScaleFamily(scale=estimate(E1.est02)),
                neighbor=ContNeighborhood(radius=0.5))
E1.IC10 <- optIC(model=E1.Rob32, risk=asMSE())
(E1.est21 <- oneStepEstimator(E1.x, IC=E1.IC10, start=E1.est02))
(E1.est22 <- kStepEstimator(E1.x, IC=E1.IC10, start=E1.est02, steps = 3))

## its simpler to use function roptest
(E1.est23 <- roptest(E1.x, ExpScaleFamily(), eps = 0.05, steps = 3))

## comparison
estimate(E1.est11)
estimate(E1.est12)
estimate(E1.est13)
estimate(E1.est21)
estimate(E1.est22)
estimate(E1.est23)

## confidence intervals
confint(E1.est11, symmetricBias())
confint(E1.est12, symmetricBias())
confint(E1.est13, symmetricBias())
confint(E1.est21, symmetricBias())
confint(E1.est22, symmetricBias())
confint(E1.est23, symmetricBias())


## 4. one-step estimation: radius interval
E1.IC11 <- radiusMinimaxIC(L2Fam=ExpScaleFamily(scale=estimate(E1.est01)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(E1.est31 <- oneStepEstimator(E1.x, IC=E1.IC11, start=E1.est01))
(E1.est32 <- kStepEstimator(E1.x, IC=E1.IC11, start=E1.est01, steps = 3))

## its simpler to use function roptest
(E1.est33 <- roptest(E1.x, ExpScaleFamily(), eps.upper = 0.5, distance = KolmogorovDist,
                     steps = 3))

E1.IC12 <- radiusMinimaxIC(L2Fam=ExpScaleFamily(scale=estimate(E1.est02)),
                neighbor=ContNeighborhood(), risk=asMSE(), loRad=0, upRad=Inf)
(E1.est41 <- oneStepEstimator(E1.x, IC=E1.IC12, start=E1.est02))
(E1.est42 <- kStepEstimator(E1.x, IC=E1.IC12, start=E1.est02, steps = 3))

## its simpler to use function roptest
(E1.est43 <- roptest(E1.x, ExpScaleFamily(), eps.upper = 0.5, steps = 3))

## comparison
estimate(E1.est31)
estimate(E1.est32)
estimate(E1.est33)
estimate(E1.est41)
estimate(E1.est42)
estimate(E1.est43)

## confidence intervals
confint(E1.est31, symmetricBias())
confint(E1.est32, symmetricBias())
confint(E1.est33, symmetricBias())
confint(E1.est41, symmetricBias())
confint(E1.est42, symmetricBias())
confint(E1.est43, symmetricBias())
