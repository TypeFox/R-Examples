nodeSpacer <- function(K, leftBound, rightBound, nodes, density, kShift) {

  # check inputs
  if (nargs() != 6) {
    stop("incorrect number of arguments")
  }
  if (length(K) != 1 | !is.numeric(K) | length(leftBound) != 1 | 
        !is.numeric(leftBound) | length(rightBound) != 1 | 
        !is.numeric(rightBound) | length(nodes) != 1 | !is.numeric(nodes) 
        | length(density) != 1 | !is.numeric(density) | length(kShift) 
        != 1 | !is.numeric(kShift)) {
    stop("arguments must be numeric and of length 1")
  }
  if (!(leftBound < K && K < rightBound)) {
    stop("must have that leftBound < K < rightBound")
  }
  if (nodes != round(nodes) | nodes <= 2) {
    stop("nodes must be an integer greater than 2")
  }

	# create node vector
	if (density == 0) {
		nodeSpacer <- seq(leftBound, rightBound, by=((rightBound - leftBound) 
                    / (nodes-1)))
	} else if (density > 0) {
	  # (see Hout and Foulon (2010), pp. 306)
		c <- (K - leftBound) / density
		dxi <- (1 / (nodes-1)) * (asinh((rightBound - K) / c) - asinh(-(K - 
                    leftBound) / c))
		xi <- asinh(-(K - leftBound) / c) + seq(0, nodes-1, by=1) * dxi
		nodeSpacer <- K + c * sinh(xi)
	} else {
		stop("invalid density parameter. Must be >= 0")
	}

	# mesh shifting
	if (kShift == 0) {
		# do nothing
	} else if (kShift == 1) {
	  # (see Pooley, Vetzal, and Forsyth (2002), pp. 30) 
		for (i in 1:length(nodeSpacer)) {
			if (nodeSpacer[i] < K & nodeSpacer[i+1] >= K) {
				mid <- (nodeSpacer[i+1] + nodeSpacer[i]) / 2
				nodeSpacer <- nodeSpacer * (1 + (K - mid) / mid)
				break
			} 
		}
	} else if (kShift == 2) {
		for (i in 1:length(nodeSpacer)) {
			if (nodeSpacer[i] < K & nodeSpacer[i+1] >= K) {
				nodeSpacer <- nodeSpacer * (1 + (K - nodeSpacer[i]) / nodeSpacer[i])
				break
			} 
		}
	} else {
		stop("invalid centering parameter. Must be 0, 1, or 2")
	}

	# return node vector
	nodeSpacer

}