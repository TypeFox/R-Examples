.getComp <- function(L2deriv, DistrSymm, L2derivSymm,
             L2derivDistrSymm){
        nrvalues <- length(L2deriv)
        z.comp <- rep(TRUE, nrvalues)
        A.comp <- matrix(TRUE, ncol = nrvalues, nrow = nrvalues)
        for(i in 1:nrvalues){
            if(is(L2derivDistrSymm[[i]], "SphericalSymmetry"))
                if(L2derivDistrSymm[[i]]@SymmCenter == 0)
                    z.comp[i] <- FALSE
        }
        for(i in 1:(nrvalues-1))
            for(j in (i+1):nrvalues){
                if(is(DistrSymm, "SphericalSymmetry")){
                    if((is(L2derivSymm[[i]], "OddSymmetric") & is(L2derivSymm[[j]], "EvenSymmetric"))
                       | (is(L2derivSymm[[j]], "OddSymmetric") & is(L2derivSymm[[i]], "EvenSymmetric")))
                        if((L2derivSymm[[i]]@SymmCenter == L2derivSymm[[j]]@SymmCenter)
                           & (L2derivSymm[[i]]@SymmCenter == DistrSymm@SymmCenter))
                            A.comp[i,j] <- FALSE
                }
            }
        A.comp[col(A.comp) < row(A.comp)] <- A.comp[col(A.comp) > row(A.comp)]
        return(list(A.comp = A.comp, z.comp = z.comp))

}
