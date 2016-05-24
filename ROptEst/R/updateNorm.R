#setMethod("updateNorm", "NormType", function(normtype, ...) normtype)
#setMethod("updateNorm", "InfoNorm", function(normtype, FI, ...)
#           {QuadForm(normtype) <- PosSemDefSymmMatrix(FI); normtype})
setMethod("updateNorm", "SelfNorm", function(normtype, L2, neighbor, biastype, 
                         Distr, V.comp, cent, stand,  w)
           {Cv <- getInfV(L2deriv = L2, neighbor = neighbor, 
                       biastype = biastype, Distr = Distr, 
                       V.comp = V.comp, cent = cent, stand = stand,  w = w)
            QuadForm(normtype) <- PosSemDefSymmMatrix(solve(Cv)) 
            normtype})

                                                       