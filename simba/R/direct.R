"direct" <-
	function(coord, listout=FALSE) {
	   deg <- function(x){(x * 180)/pi}
	x <- coord[,1]
	y <- coord[,2]
	anz <- nrow(coord)
	rn <- rownames(coord)
	tmp1 <- outer(x,x,"-")
	tmp2 <- outer(y,y,"-")
	tmp3 <- round(deg(atan(tmp2/tmp1)), 1)
	tmp3[(tmp3 >= -22.5)&(tmp3 <= 22.5)] <- 100
	tmp3[(tmp3 < -22.5)&(tmp3 >= -67.5)] <- 200
	tmp3[(tmp3 <= -67.5)&(tmp3 >= -90)] <- 300
	tmp3[(tmp3 <= 90)&(tmp3 >= 67.5)] <- 300
	tmp3[(tmp3 <= 67.5)&(tmp3 > 22.5)] <- 400
	tmp3[(tmp3 <= 22.5)&(tmp3 >= -22.5)] <- 100
	tmp3[(tmp3 >= -67.5)&(tmp3 < -22.5)] <- 200
	tmp3[(tmp3 >= -90)&(tmp3 <= -67.5)] <- 300
	tmp3[(tmp3 >= 67.5)&(tmp3 <= 90)] <- 300
	tmp3[(tmp3 > 22.5)&(tmp3 <= 67.5)] <- 400
	d <- as.dist(tmp3)
	attr(d, "Size") <- anz
    attr(d, "Labels") <- rn
    attr(d, "call") <- match.call()
    attr(d, "method") <- "direct"
    class(d) <- "dist"
	if (listout) {
		d <- liste(d, entry="direct")
		}
	return(d)
}