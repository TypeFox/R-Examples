elevateBezierDegree <- function(p, deg){
	# Math explained at http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-elev.html

	# GET CURRENT DEGREE
	cdeg <- length(p) - 1

	# CHECK IF CURRENT DEGREE ALREADY MATCHES NEW DEGREE
	if(cdeg == deg) return(p)

	# CHECK IF NEW DEGREE IS LESS THAN CURRENT DEGREE
	if(deg < cdeg) stop(paste0("Input degree (", deg, ") cannot be less than the degree of the input control points (", cdeg, ")"))

	# ELEVATE BY ONE DEGREE SEQUENTIALLY UNTIL REACHES NEW DEGREE
	for(ndeg in (cdeg + 1):deg){i <- 1:(ndeg - 1);p <- c(p[1], (i / ndeg)*p[i] + (1 - (i/ndeg))*p[i+1], p[length(p)])}

	return(p)	
}