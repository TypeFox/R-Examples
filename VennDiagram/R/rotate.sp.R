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

### ROTATION WITHOUT NEED OF OTHER ARGUMENTS ######################################################
rotate.sp <- function(area.vector, rotation, reverse, additional.rot = FALSE, additional.o7 = 1:7, additional.o3 = 1:3) {

	rot.f.7 <- list(1:7,c(3,6,7,2,5,4,1),c(7,4,1,6,5,2,3));
	rot.f.3 <- list(1:3,c(2,3,1),c(3,1,2));
	rot.r.7 <- list(c(3,2,1,6,5,4,7),c(7,6,3,4,5,2,1),c(1,4,7,2,5,6,3));
	rot.r.3 <- list(c(2,1,3),3:1,c(1,3,2));
	
	#Add an additional rotation to the orders but not the areas for cases 011A and 111A
	#Permutations are associative, so this allows for the chaining of the permutations (think of boxes with pointers above indicating the permutation)
	if(additional.rot)
	{
		if (reverse) {
			area.rot <- rot.r.7[[rotation]];
			order.7 <- additional.o7[rot.r.7[[rotation]]];
			order.3 <- additional.o3[rot.r.3[[rotation]]];
		}
		else {
			area.rot <- rot.f.7[[rotation]];
			order.7 <- additional.o7[rot.f.7[[rotation]]];
			order.3 <- additional.o3[rot.f.3[[rotation]]];
		}
		return(
			list(
				areas = area.vector[area.rot],
				o7 = order.7,
				o3 = order.3
				)
			);
		}
	
	#If not adding an additional rotation to only the orders but not the areas, then continue as usual

	if (reverse) {
			order.7 <- rot.r.7[[rotation]];
			order.3 <- rot.r.3[[rotation]];
		}
	else {
			order.7 <- rot.f.7[[rotation]];
			order.3 <- rot.f.3[[rotation]];
		}

	return(
		list(
			areas = area.vector[order.7],
			o7 = order.7,
			o3 = order.3
			)
		);

	}
