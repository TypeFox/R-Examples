
.sign.max <- function(x)
 {
 ifelse(identical(max(x), max(abs(x))), 1, -1)
 }


.standard.B <- function(B)
    {
    B.norm <- sqrt(rowSums((B)^2))
    B.r <- sweep(B, 1, B.norm, "/")
    row.signs <- apply(B.r, 1, .sign.max)
    B.sr <- sweep(B.r, 1, row.signs,"*")
    row.order <- apply(B.sr,1,max)
    B.osr <- B.sr[order(row.order, decreasing=TRUE), ]
    
    return(B.osr)
    }
