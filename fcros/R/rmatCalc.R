rmatCalc = function(fvect, n, m1, m2, rvect, FC) {
    .C("rmat", fvectC = as.double(fvect),
               nC = as.integer(n),
               m1C = as.integer(m1),
               m2C = as.integer(m2),
               rvectC = as.double(rvect),
               FCC = as.double(FC));
}
