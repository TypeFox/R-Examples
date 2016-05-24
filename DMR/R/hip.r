hip = function(n.levels, n.cont){
    B = 0
    for (m in 1:length(n.levels)){
        k = n.levels[m] - 1
        k2 = k*(k + 1)/2
        if (k == 1){
            A = 1
        } else {
            A = diag(1, k2, k2)[, 1:k]
            out = matrix(0, 2, k2)
            out[1, 1:k] = 1:k
            ind = 0
            for (i in 2:k){
                j = 1:(i - 1)
                out[1, (k + ind + 1):(k + ind + i - 1)] = rep(i, times = length(j))
                out[2, (k + ind + 1):(k + ind + i - 1)] = j
                ind = ind + i - 1
            }
            for (i in (k + 1):ncol(out)){
                A[i, out[1, i]] = -1
                A[i, out[2, i]] = 1
            }
        }
        B = adiag(B, A)
    }
    B = adiag(B, diag(n.cont))
    return(B[-1, ])
}
