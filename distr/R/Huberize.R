
setMethod("Huberize", "AcDcLcDistribution",
          function(object, lower, upper, withSimplify = getdistrOption("simplifyD"))
                {Mi <- Minimum( object, Dirac(upper),
                                withSimplify = withSimplify)
                 M0 <- Maximum(Dirac(lower),Mi,
                        withSimplify = withSimplify)
                 if(is(object@Symmetry,"SphericalSymmetry"))
                    if(.isEqual(lower+upper,2*SymmCenter(object@Symmetry))) 
                       M0@Symmetry <- SphericalSymmetry(SymmCenter(object@Symmetry))        
                 M0})
