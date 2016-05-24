HIsurf <-
function (G, P, type, size) 
{
    S <- H <- seq(from = 0, to = 1, length.out = size)
    surf <- matrix(NA, ncol = size, nrow = size)
    for (i in 1:size) {
        for (j in 1:size) {
            surf[i, j] <- HILL(c(S[i], H[j]), G, P, type)
        }
    }
    for (i in 1:size) {
        for (j in 1:size) {
            if (H[j] > min(2 * S[i], 2 - 2 * S[i])) 
                surf[i, j] <- NA
        }
    }
    surf
}
