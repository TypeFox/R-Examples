#Called if the sp.case.name is not one of these: 022AAAO, 022AAOO, 023, 032, 120, 121AO, 122AAOO, 130
#if(!sp.case.name %in% c("022AAAO", "022AAOO", "023", "032", "120", "121AO", "122AAOO", "130"))

draw.sp.case.preprocess <- function(
	sp.case.name,
	a1,
	a2,
	a3,
	a4,
	a5,
	a6,
	a7,
	category = rep("", 3),
	reverse = FALSE,
	cat.default.pos = "outer",
	lwd = rep(2, 3),
	lty = rep("solid", 3),
	col = rep("black", 3),
	label.col = rep("black", 7),
	cex = rep(1, 7),
	fontface = rep("plain", 7),
	fontfamily = rep("serif", 7),
	cat.pos = c(-40, 40, 180),
	cat.dist = c(0.05, 0.05, 0.025),
	cat.col = rep("black", 3),
	cat.cex = rep(1, 3),
	cat.fontface = rep("plain", 3),
	cat.fontfamily = rep("serif", 3),
	cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
	cat.prompts = FALSE,
	fill = NULL,
	alpha = rep(0.5, 3),
	scaled = TRUE,
	offset = 0,
	sep.dist = rep(0.05, 3),
	print.mode="raw",
	sigdigs=3,
	...
	) {

	if(sp.case.name == "130")
	{
		sep.dist = rep(0.05,3);
	}
	
	area.zeroes <- list("023" = c(1,3,4,6,7), "022AAOO" = c(1,3,4,6), "001" = c(7), "010" = c(2), "011O" = c(2,7), "100" = c(5), "112AA" = c(1,4,5,7), "021AA" = c(4,6,7), "012AA" = c(1,4,7), "122AAOO" = c(1,3,4,5,6), "033" = c(1,2,3,4,6,7), "120" = c(4,5,6), "022AAAO" = c(3,4,6,7), "111A" = c(4,5,7), "011A" = c(4,7), "130" = c(2,4,5,6), "032" = c(2,3,4,6,7), "121AO" = c(3,4,5,6), "110" = c(4,5));
	
	area.labels <- list("012AA" = c(2,3,6), "122AAOO" = c(2,2,7), "033" = c(5,5,5), "120" = c(1,3,7), "022AAAO" = c(1,2,5), "111A" = c(1,3,6), "011A" = c(1,3,6), "130" = c(1,3,7), "032" = c(1,5,5), "121AO" = c(1,2,7), "110" = c(1,3,7), "023" = c(2,2,5), "022AAOO" = c(2,2,7), "001" = c(1,3,5), "010" = c(1,3,5), "011O" = c(1,3,5), "100" = c(1,3,7), "112AA" = c(2,3,6), "021AA" = c(1,3,5));

	a.x <- list("012AA" = c(0,0.3,0.5,0,0.5,0.7,0), "111A" = c(0.27,0.5,0.68,0,0,0.76,0), "011O" = c(0.175,0,0.825,0.325,0.5,0.675,0), "021AA" = c(0.2,0.5,0.8,0,0.5,0,0), "033" = c(0,0,0,0,0.5,0,0), "010" = c(0.18,0,0.82,0.32,0.5,0.68,0.5), "110" = c(0.175,0.375,0.5,0,0,0.625,0.825), "112AA" = c(0,0.34,0.5,0,0,0.66,0), "011A" = c(0.2,0.425,0.875,0,0.575,0.725,0), "100" = c(0.31,0.495,0.68,0.41,0,0.59,0.5), "001" = c(0.3,0.5,0.7,0.32,0.5,0.68,0));

	a.y <- list("012AA" = c(0,0.5,0.775,0,0.5,0.5,0), "111A" = c(0.5,0.5,0.705,0,0,0.5,0), "011O" = c(0.5,0,0.5,0.5,0.5,0.5,0), "021AA" = c(0.5,0.675,0.5,0,0.5,0,0), "033" = c(0,0,0,0,0.5,0,0), "010" = c(0.5,0,0.5,0.54,0.51,0.54,0.75), "110" = c(0.5,0.5,0.5,0,0,0.5,0.5), "112AA" = c(0,0.5,0.725,0,0,0.5,0), "011A" = c(0.5,0.5,0.5,0,0.5,0.5,0), "100" = c(0.66,0.66,0.66,0.5,0,0.5,0.335), "001" = c(0.35,0.33,0.35,0.58,0.55,0.58,0));

	centre.x <- list("012AA" = c(0.4,0.5,0.6), "111A" = c(0.32,0.68,0.76), "011O" = c(0.35,0.65,0.5), "021AA" = c(0.35,0.65,0.5), "033" = c(0.5,0.5,0.5), "010" = c(0.35,0.65,0.5), "110" = c(0.25,0.5,0.75), "112AA" = c(0.34,0.5,0.66), "011A" = c(0.35,0.65,0.65), "100" = c(0.31,0.68,0.5), "001" = c(0.4,0.6,0.5));

	centre.y <- list("012AA" = c(0.5,0.5,0.5), "111A" = c(0.5,0.5,0.5), "011O" = c(0.5,0.5,0.5), "021AA" = c(0.5,0.5,0.5), "033" = c(0.5,0.5,0.5), "010" = c(0.5,0.5,0.55), "110" = c(0.5,0.5,0.5), "112AA" = c(0.5,0.5,0.5), "011A" = c(0.5,0.5,0.5), "100" = c(0.66,0.66,0.333), "001" = c(0.5,0.5,0.55));

	radii <- list("012AA" = c(0.2,0.35,0.2), "111A" = c(0.27,0.27,0.14), "011O" = c(0.25,0.25,0), "021AA" = c(0.3,0.3,0.12), "033" = c(0.4,0.4,0.4), "010" = c(0.25,0.25,0.25), "110" = c(0.2,0.2,0.2), "112AA" = c(0.15,0.35,0.15), "011A" = c(0.3,0.3,0.15), "100" = c(0.216,0.216,0.216), "001" = c(0.25,0.25,0));


	######################Rotations
	
	#Need to check certain special rotations
	if(sp.case.name == "111A" || sp.case.name == "011A"){
		for(i in 1:3){
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), i, reverse);
			if (tmp$areas[7] == 0 & tmp$areas[4] == 0) { break; }
			if (tmp$areas[7] == 0 & tmp$areas[6] == 0) { 
				tmp <- VennDiagram::rotate.sp(tmp$areas, 1, reverse = TRUE,additional.rot = TRUE,additional.o7=tmp$o7,additional.o3=tmp$o3);
				break;
			}
		}
	} else if(sp.case.name == "121AO"){ #Need to check all possible rotations
		for (i in 1:6) {
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), (i-1) %% 3 + 1, reverse = (i>3));
			if (0 == tmp$areas[3] & 0 == tmp$areas[4] & 0 == tmp$areas[5] & 0 == tmp$areas[6]) { break; }
		}
	} else if(sp.case.name == "022AAAO"){#Need to make sure reverse = FALSE
		for (i in 1:3) {
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), i, reverse = FALSE);
			if (0 == tmp$areas[3] & 0 == tmp$areas[4] & 0 == tmp$areas[6] & 0 == tmp$areas[7]) { break; }
		}
	} else {#Normal rotations
		break.ind <- c(2,4,5,6);#Break if these tmp$areas are equal to zero
		#Get the break.ind by indexing a list by the sp.case.name
		break.ind <- area.zeroes[[sp.case.name]];
		for (i in 1:3) {
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), i, reverse);
			if (all(tmp$areas[break.ind]==0)) { break; }
		}
	}

	a1 <- tmp$areas[1];
	a2 <- tmp$areas[2];
	a3 <- tmp$areas[3];
	a4 <- tmp$areas[4];
	a5 <- tmp$areas[5];
	a6 <- tmp$areas[6];
	a7 <- tmp$areas[7];

	# 3-vector rotations
	fill <- fill[tmp$o3];
	cat.col <- cat.col[tmp$o3];
	category <- category[tmp$o3];
	lwd <- lwd[tmp$o3];
	lty <- lty[tmp$o3];
	col <- col[tmp$o3];
	alpha <- alpha[tmp$o3];
	cat.dist <- cat.dist[tmp$o3];
	cat.cex <- cat.cex[tmp$o3];
	cat.fontface <- cat.fontface[tmp$o3];
	cat.fontfamily <- cat.fontfamily[tmp$o3];
	cat.just <- cat.just[tmp$o3];

	# 7-vector rotations
	label.col <- label.col[tmp$o7];
	cex <- cex[tmp$o7];
	fontface <- fontface[tmp$o7];
	fontfamily <- fontfamily[tmp$o7];

	a1.x.pos <- 0;
	a1.y.pos <- 0;
	a2.x.pos <- 0;
	a2.y.pos <- 0;
	a3.x.pos <- 0;
	a3.y.pos <- 0;
	a4.x.pos <- 0;
	a4.y.pos <- 0;
	a5.x.pos <- 0;
	a5.y.pos <- 0;
	a6.x.pos <- 0;
	a6.y.pos <- 0;
	a7.x.pos <- 0;
	a7.y.pos <- 0;

	###########Calculations of [xy].centre[1-3] and a[1-7].[xy].pos

	########### Get the areas and positions by indexing into a vector using the sp.case.name. Don't have to do other caculations because these aren't the scaled cases. This is done directly in the return function call

	r1 <- radii[[sp.case.name]][1];
	r2 <- radii[[sp.case.name]][2];
	r3 <- radii[[sp.case.name]][3];

	a.list = c(r1, r2, r3);
	b.list = c(r1, r2, r3);

	if(sp.case.name == "001" || sp.case.name == "011O")
	{
		a.list = c(r1,r2,0.25);
		if(sp.case.name == "001") {
			b.list = c(r1,r2,0.18);
		}
		else {#sp.case.name == "011O"
			b.list = c(r1,r2,0.2);
		}
	}
	
	if(!sp.case.name %in% c("011A","022AAAO","023","032","033","111A","121AO"))
	{
		reverse = FALSE;
	}
	
	straight.reverse = FALSE;
	
	if(!sp.case.name %in% c("011A","022AAAO","022AAOO","111A","120","121AO","122AAOO"))
	{
		straight.reverse = TRUE;
	}

	return(
		VennDiagram::draw.sp.case(
			area.list = c(a1, a2, a3, a4, a5, a6, a7),
			enabled.areas = setdiff(1:7,unlist(area.zeroes[sp.case.name])),#The enabled areas are those with non-zero areas
			area.x = a.x[[sp.case.name]],
			area.y = a.y[[sp.case.name]],
			attach.label.to = area.labels[[sp.case.name]],
			x.centres = centre.x[[sp.case.name]],
			y.centres = centre.y[[sp.case.name]],
			a.list = a.list,
			b.list = b.list,
			straight.reverse = straight.reverse,
			reverse = reverse,
			category = category,
			cat.default.pos = cat.default.pos,
			lwd = lwd,
			lty = lty,
			col = col,
			label.col = label.col,
			cex = cex,
			fontface = fontface,
			fontfamily = fontfamily,
			cat.pos = cat.pos,
			cat.dist = cat.dist,
			cat.col = cat.col,
			cat.cex = cat.cex,
			cat.fontface = cat.fontface,
			cat.fontfamily = cat.fontfamily,
			cat.just = cat.just,
			fill = fill,
			alpha = alpha,
			print.mode=print.mode,
			sigdigs=sigdigs,
			...
			)
		);
	}
