getsizeratios <-
function (bycx, ac, method = "fast") 
{
    K1 = nrow(ac)
    nc = length(bycx)
    if (method == "fast") {
        bwx = bycx %*% t(ac) %*% solve(ac %*% t(ac))
        sizeratios = rep(0, K1)
        for (i in 1:K1) {
            divisor = sqrt(1 + sum(bwx[1:i]^2))
            betachat = bycx - bwx[, 1:i, drop = FALSE] %*% ac[1:i, 
                , drop = FALSE]
            scaledbetac = sqrt(nc/(nc - i)) * betachat/divisor
            sizeratios[i] = median(abs(as.vector(ac[i, ])/as.vector(scaledbetac)))
        }
        return(sizeratios)
    }
    else if (method == "leave1out") {
        bycx = as.vector(bycx)
        A = Aj = ajtB = bwx = t(bycx) %*% t(ac)
        B = Bj = ac %*% t(ac)
        sizeratios = rep(0, K1)
        scaledbetac = rep(0, nc)
        temp = .C("getsizeratios", as.double(A), as.double(B), 
            as.double(Aj), as.double(Bj), as.double(bycx), as.double(ac), 
            as.double(ajtB), as.double(sizeratios), as.double(scaledbetac), 
            as.double(bwx), as.integer(K1), as.integer(nc))
        return(temp[[8]])
    }
    else return(0)
}
