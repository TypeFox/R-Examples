zerobasis = function(k) {
#  sets up a k by k-1 matrix with orthonormal columns
#  using the first k non-constant fourier basis function
#  values at 0.5, ..., k-0.5
fourierbasis = create.fourier.basis(k,k);
tk           = 0:(k-1) + 0.5
fbasmat      = eval.basis(tk, fourierbasis)
zerobasmat   = fbasmat[,2:k]
return(zerobasmat)
}



