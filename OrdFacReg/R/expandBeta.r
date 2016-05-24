expandBeta <- function(beta.col, sums, JJs.A){

eps <- 0
mini <- min(JJs.A[, 4])

## expand estimate to receive estimate in original dimension back
sums <- cbind(sums, sums[, 3] - sums[, 2] + 1)
if (mini == 1){beta <- NULL} else {beta <- beta.col[1:(mini - 1)]}
ind <- mini

for (s in 1:length(unique(JJs.A[, 2]))){

    JJs.A.s <- matrix(JJs.A[JJs.A[, 2] == s, ], ncol = 5)
    
    if (JJs.A.s[1, 5] == 1){
    
        if (nrow(JJs.A.s) == 1){n.act <- 1} else {
    
        n.act <- cbind(JJs.A.s[, 5], c(JJs.A.s[2:(nrow(JJs.A.s)), 5], NA))
        n.act <- min(nrow(JJs.A.s), (1:nrow(JJs.A.s))[(n.act[, 1] == n.act[, 2]) == 0], na.rm = TRUE)
        }    
        beta <- c(beta, rep(eps, n.act))
        JJs.A.s <- JJs.A.s[-c(1:n.act), ]
    }

    if (length(JJs.A.s) > 0){
    sums.tmp <- matrix(sums[sums[, 1] == s, ], ncol = 4)
    n.col <- length(sums.tmp[, 1])   

    if (n.col > 0){    
        beta.col.tmp <- beta.col[ind - 1 + 1:n.col]            
        for (i in 1:n.col){
            beta <- c(beta, rep(beta.col.tmp[i], sums.tmp[i, 4]))
            } # end for i
        
        ind <- ind + n.col} else {
            nr.tmp <- length(JJs.A.s[, 1])
            beta <- c(beta, beta.col[ind + 0:(nr.tmp - 1)])
            ind <- ind + nr.tmp           
        }
    } # if length(JJs.A.s) > 0
} # end for s

beta <- matrix(beta, ncol = 1)

return(list("beta" = beta))
}
