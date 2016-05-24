`pimmlawton` <-
function (mat, N = 1, omni.i = NA, omni.j = NA, omega = NULL) 
{
    S <- nrow(mat)
    if (is.na(omni.i)) {
        out <- matrix(NA, nrow = N, ncol = 4)
        colnames(out) <- c("DomEig", "Im", "IntraDD", "I")
        for (n in 1:N) out[n, ] <- {
            M = runif(S^2) * mat
            eigs <- eigen(M)[["values"]]
            mx <- which.max(Re(eigs))
            deM = Re(eigs)[mx]
            deMi = Im(eigs)[mx]
            intraNDD <- sqrt(sum(diag(M)^2)/S)
            diag(M) <- 0
            IS = sqrt(sum(M^2)/(S * (S - 1)))
            c(deM, deMi, intraNDD, IS)
        }
    }
    else {
        out <- matrix(NA, nrow = N, ncol = 5)
        colnames(out) <- c("DomEig", "Im", "IntraDD", "I", "I.omni")
        for (n in 1:N) out[n, ] <- {
            M = runif(S^2) * mat
            if (!is.null(omega)) {
                M[omni.i, omni.j] <- omega * M[omni.i + 1, omni.j]
                M[omni.i + 1, omni.j] <- (1 - omega) * M[omni.i + 
                  1, omni.j]
                M[omni.j, omni.i] <- omega * M[omni.j, omni.i + 
                  1]
                M[omni.j, omni.i + 1] <- (1 - omega) * M[omni.j, 
                  omni.i + 1]
            }
            eigs <- eigen(M)[["values"]]
            mx <- which.max(Re(eigs))
            deM = Re(eigs)[mx]
            deMi = Im(eigs)[mx]
            intraNDD <- sqrt(sum(diag(M)^2)/S)
            diag(M) <- 0
            IS = sqrt(sum(M^2)/(S * (S - 1)))
            omnivory <- sqrt(mean(c(M[omni.i, omni.j], M[omni.j, 
                omni.i])^2))
            c(deM, deMi, intraNDD, IS, omnivory)
        }
    }
    return(as.data.frame(out))
}
