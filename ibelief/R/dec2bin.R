dec2bin  <-  function(nb, n) {
    v = nb
    if (v == 0) {
        pos = 1
    } else {
        pos = floor(log(v, 2)) + 1
    }
    
    if (!missing(n)) {
        if (nb >= 2^n) {
            warning("n is too small \n")
        } else {
            pos = n
        }
    }
    bin = rep(0, pos)
    i = 1
    while (v != 0) {
        p = v%/%2
        if (2 * p == v) {
            bin[i] = 0
        } else {
            bin[i] = 1
        }
        v = p
        i = i + 1
    }
    return(bin)  #return binary encoding
} 

