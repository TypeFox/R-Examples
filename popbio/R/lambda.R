lambda<-function(A)
{
    ev <- eigen(A)
    # R sorts eigenvalues in decreasing order, according to Mod(values)
    #  ususally dominant eigenvalue is first (ev$values[1]), except for
    #  imprimitive matrices with d eigenvalues of equal modulus
    # this should work for most cases
    lmax <- which.max(Re(ev$values))
    lambda <- Re(ev$values[lmax])
    lambda
}
