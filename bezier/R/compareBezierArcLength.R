compareBezierArcLength <- function(p, l, t1 = 0, t2 = NULL, deg = NULL, relative.min.slope = 1e-6, absolute.min.slope = 0){
	abs(bezierArcLength(p=p, t1=t1, t2=t2, deg=deg, relative.min.slope=relative.min.slope, absolute.min.slope=absolute.min.slope)$arc.length - l)
}