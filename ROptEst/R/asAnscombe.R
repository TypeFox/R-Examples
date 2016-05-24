## generating function
asAnscombe <- function(eff = .95, biastype = symmetricBias(),
                     normtype = NormType()){
   new("asAnscombe", eff = eff, biastype = biastype,
                   normtype = normtype) }

## access method
setMethod("eff", "asAnscombe", function(object) object@eff)
