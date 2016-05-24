Oversample = function(x, n){
    # function oversamples trace x by a factor of n.
    # use nearest-neighbor interpolation to oversample.  this works, and is
    # causal (unlike sinc interpolation, which is achieved by zero-padding
    # the trace's fft).
    # n must be a positive integer.  n = 1 means the original trace is returned.

    if(n != floor(n) || n < 1){
        return(NaN)
    }
    
    m = length(x)
    y = rep(0, n * m)

    for(i in 1:n){
        y[(1:m - 1) * n + i] = x
    }

    return(y)
}
