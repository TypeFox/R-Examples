"aslopect" <-
function(asp, slo, names=rownames(asp), fc=FALSE, listout=FALSE) {
    if(fc) {
    	asp[asp > 180] <- -(180 - (asp[asp > 180] - 180))
    	}
    asp[is.na(asp)] <- 0
    slo <- 90-slo
    anz <- length(asp)
    t1 <- function(x, y){sin(x*pi/180)*sin(y*pi/180)}
    t2 <- function(x, y){cos(x*pi/180)*cos(y*pi/180)}
    t3 <- function(x, y){cos(abs(x-y)*pi/180)}
    t1e <- outer(slo, slo, FUN = "t1")
    t2e <- outer(slo, slo, FUN = "t2")
    t3e <- outer(asp, asp, FUN = "t3")
    erg <- acos(t1e + t2e * t3e)/pi
    erg <- as.dist(erg)
    attr(erg, "Size") <- anz
    attr(erg, "Labels") <- names
    attr(erg, "call") <- match.call()
    attr(erg, "method") <- "aslopect"
    class(erg) <- "dist"
    if (listout) {
    	erg <- liste(erg, entry="aslopect")
    	}
    return(erg)
    }