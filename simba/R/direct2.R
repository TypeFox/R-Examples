"direct2" <-
	function(coord, listout=FALSE) {
	   deg <- function(x){(x * 180)/pi}
	x <- coord[,1]
	y <- coord[,2]
	anz <- nrow(coord)
	rn <- rownames(coord)
	tmp1 <- outer(x,x,"-")
	tmp2 <- outer(y,y,"-")
	tmp3 <- round(deg(atan(tmp2/tmp1)), 1)
	tmp3[(tmp3 <= 0)&(tmp3 >= -30)] <- 100
	tmp3[(tmp3 <= -30)&(tmp3 >= -60)] <- 200
	tmp3[(tmp3 <= -60)&(tmp3 >= -90)] <- 300
	tmp3[(tmp3 <= 90)&(tmp3 >= 60)] <- 400
	tmp3[(tmp3 <= 60)&(tmp3 >= 30)] <- 500
	tmp3[(tmp3 <= 30)&(tmp3 >= 0)] <- 600
	tmp3[(tmp3 >= -30)&(tmp3 <= 0)] <- 100
	tmp3[(tmp3 >= -60)&(tmp3 <= -30)] <- 200
	tmp3[(tmp3 >= -90)&(tmp3 <= -60)] <- 300
	tmp3[(tmp3 >= 60)&(tmp3 <= 90)] <- 400
	tmp3[(tmp3 >= 30)&(tmp3 <= 60)] <- 500
	tmp3[(tmp3 >= 0)&(tmp3 <= 30)] <- 600
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
	