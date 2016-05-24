hp.nqx <-
function(H.out, age = seq(0, 85, 1)) 
{
    H.new.hat <- matrix(NA, nrow = nrow(H.out), ncol = length(age))
    for (i in 1:nrow(H.out)) {
        H.new.hat[i, ] <- mod8p(theta = H.out[i, ], x = age)
    }
    ans <- H.new.hat
    return(ans)
}

