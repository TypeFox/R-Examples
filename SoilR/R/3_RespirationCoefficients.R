#
# vim:set ff=unix expandtab ts=2 sw=2:
RespirationCoefficients=function #helper function to compute respiration coefficients
###This function computes the respiration coefficients as function of time for all pools
### according to the given matrix A
(A ##<< A matrix valued function representing the model.
 ){
   nr=nrow(A(1))
    testvec=matrix(nrow=1,ncol=nr,1)
    ### The respiration coefficient of the nth pool is the negative 
    ### sum of the elements of the n th column of the model matrix
    rcoeffs= function(t){-testvec%*%A(t)}
    return(rcoeffs)
    ### A vector valued function of time containing the respiration coefficients for all pools.
}
 
