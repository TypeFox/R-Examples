NCEP.uv.revert <- function(spd, dir, radians=FALSE){
	if(radians == FALSE) {
		dir <- dir * (pi/180)
		}
		U <- sin(dir)*spd
		V <- cos(dir)*spd
		UV <- data.frame(U, V)
return(UV)
}