#
# vim:set ff=unix expandtab ts=2 sw=2:
TimeMapWarningOperators <- function(){
            "The use of object of class TimeMap is deprecated.
            At the moment we cast TimeMap objects to the new class BoundLinDecompOp.
            which replaces TimeMap as class of the Decomposition Operator argument 
            To get rid of this warning adapt your code to use a BoundLinDecompOp instead of a TimeMap.
            You can als use an object of class ConstLinDecompOp.
            Other classes may be implemented in the future." 
}
TimeMapWarningInFluxes <- function(){
            "The use of object of class TimeMap is deprecated for the InFluxes argument.
            At the moment we cast TimeMap objects to the new class BoundInFlux
            which replaces TimeMap as input to this function.
            To get rid of this warning adapt your code to use BoundInFlux yourself instead of TimeMap.
            Other classes describing InFluxes may be implemented in the future as need be." 
}
TimeMapWarningBoundFc<- function(){
            "The use of object of class TimeMap is deprecated for the BoundFc argument.
            At the moment we cast TimeMap objects to the new class BoundFc
            which replaces TimeMap as input to this function.
            To get rid of this warning adapt your code to use BoundFc yourself instead of TimeMap.
            Other classes describing the atmospheric 14C fraction may be implemented in the future as need be." 
}
WarningConstFc<- function(){
            "The use of object of class SoilR.F0 is deprecated 
            At the moment we cast SoilR.Fo objects to the new class ConstFc
            which replaces it as input to this function.
            To get rid of this warning adapt your code to use CostFc in the first place. e.g replace SoilR.F0.new(some,args )  by ConstFC(same,args)" 
}


