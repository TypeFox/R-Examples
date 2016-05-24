turnPts <- local({

    goodZeroes <- function(dp,xlo,xhi){
        zzz  <- polyroot(dp)
        rrr  <- sapply(zzz,function(z){isTRUE(all.equal(Im(z),0))})
        if(!any(rrr)) return(numeric(0))
        zzz  <- unique(Re(zzz[rrr]))
        ok   <- xlo <= zzz & zzz <= xhi
        zzz[ok]
    }

function(a,b,v,Kpa,xlo,xhi,type) {
#
# Construct the polynomial.
    q <- length(v)
    rmax  <- suppressWarnings(max(which(Kpa > 0)))
    if(rmax < 1) return(if(type=="sip") NULL else rep(list(NULL),q))
    Kpa <- Kpa[1:rmax]
    vq <- v[q]
    if(type=="sip") ply <- 0 else rslt <- vector("list",q)
    for(r in 1:rmax) {
        vqmr <- if(r < q) v[q-r] else 0
        c1 <- vqmr - vq
        if(type=="dip") {
            d1 <- if(r>1) polynomial(c(a,b)) else 1
            dply <- d1*polynomial(c(a+b*c1,b*(r+1)))
            rslt[[r]] <- goodZeroes(dply,xlo,xhi)
        } else {
            p1  <- polynomial(c(c1,r))
            p2  <- polynomial(c(a,b))
            ply <- ply + Kpa[r]*p1*p2^r
        }
    }
    if(type=="dip") return(rslt)
#
# Take the derivative of ply.
    dply <- deriv(ply)
#
# Return the "good" zeroes.
    goodZeroes(dply,xlo,xhi)
}})
