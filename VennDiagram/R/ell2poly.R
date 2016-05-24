# The VennDiagram package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION TO OBTAIN POLYGON COORDINATES SIMILAR TO A CIRCLE ####################################
ell2poly <- function(x, y, a, b, rotation, n.sides) {
	# draw an n-sided polygon that resembles an ellipse

	rotation <- rotation * pi / 180;
	# calculate the angle corresponding to each "section" of the polygon
	# (there are as many sections as there are sides in the polygon)
	theta <- 2 * pi / n.sides;
	angles <- seq(0, 2 * pi, theta);

	# initialize vectors to hold the x and y coordinates of each vertex of the polygon
	x.coord <- vector(length = n.sides + 1, mode = 'numeric');
	x.coord[1] <- x + a * cos(rotation);
	y.coord <- vector(length = n.sides + 1, mode = 'numeric');
	y.coord[1] <- y + a * sin(rotation);

	# starting from the initial point, sequentially obtain the coordinates of each vertex of the polygon and store them
	for (i in 1:n.sides) {
		x.coord[i + 1] <- x + a * cos(angles[i + 1]) * cos(rotation) - b * sin(angles[i + 1]) * sin(rotation);
		y.coord[i + 1] <- y + a * cos(angles[i + 1]) * sin(rotation) + b * sin(angles[i + 1]) * cos(rotation);
		}

	return(
		list(
			x = x.coord,
			y = y.coord
			)
		);
	}
