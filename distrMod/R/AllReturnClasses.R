#####################################################################
###
## pure return value classes as subclasses of
## L2ParamFamily, L2LocationFamily, L2ScaleFamily, L2LocationScaleFamily
## used for method dispatch later
###
#####################################################################



## Binomial family
setClass("BinomFamily",
          contains = "L2ParamFamily")

## Poisson family
setClass("PoisFamily",
          contains = "L2ParamFamily")

## neg.Binomial family
setClass("NbinomFamily",
          contains = "L2ParamFamily")

## neg.Binomial family with size
setClass("NbinomwithSizeFamily",
          contains = "L2ParamFamily")
## neg.Binomial family in different parametrization
setClass("NbinomMeanSizeFamily",
          contains = "L2ParamFamily")

## Gamma family
setClass("GammaFamily", prototype=prototype(scaleshapename=c("scale"="scale",
                                  "shape"="shape")),
          contains = "L2ScaleShapeUnion")

## Beta family
setClass("BetaFamily",
          contains = "L2ParamFamily")

## Normal location family
setClass("NormLocationFamily",
          contains = "L2LocationFamily")

## Normal scale family
setClass("NormScaleFamily",
          contains = "L2ScaleFamily")

## Exponential scale family
setClass("ExpScaleFamily",
          contains = "L2ScaleFamily")

## Lognormal scale family
setClass("LnormScaleFamily",
          contains = "L2ScaleFamily")

## Normal location and scale family
setClass("NormLocationScaleFamily",
          contains = "L2LocationScaleFamily")

## Cauchy location scale family
setClass("CauchyLocationScaleFamily",
          contains = "L2LocationScaleFamily")


