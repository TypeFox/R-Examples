#Called if the sp.case.name is not one of these: 022AAAO, 022AAOO, 023, 032, 120, 121AO, 122AAOO, 130
#if(!sp.case.name %in% c("022AAAO", "022AAOO", "023", "032", "120", "121AO", "122AAOO", "130"))

draw.sp.case.scaled <- function(
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

	######################Rotations
	
	#Need to check certain special rotations
	if(sp.case.name == "111A" || sp.case.name == "011A"){
		for(i in 1:3){
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), i, reverse);
			if (tmp$areas[7] == 0 & tmp$areas[4] == 0) { break; }
			if (tmp$areas[7] == 0 & tmp$areas[6] == 0) { 
				tmp <- VennDiagram::rotate.sp(tmp$areas, 1, reverse = TRUE);
				break;
			}
		}
	} else if(sp.case.name == "121AO"){ #Need to check all possible rotations
		for (i in 1:6) {
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), (i-1) %% 3 + 1, reverse = (i>3));
			if (0 == tmp$areas[3] & 0 == tmp$areas[4] & 0 == tmp$areas[5] & 0 == tmp$areas[6]) { break; }
		}
	} else if(sp.case.name == "022AAAO"){#Need to make sure reverse = FALSE
		for (i in 1:6) {
			tmp <- VennDiagram::rotate.sp(c(a1, a2, a3, a4, a5, a6, a7), (i-1) %% 3 + 1, reverse = i>3);
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

	########### Calculations of [xy].centre[1-3] and a[1-7].[xy].pos
	########### Calculate the areas, radii and positions for each case seperately
	
	if (sp.case.name == "022AAAO"){
		if (scaled) {
			r1 <- sqrt((a1 + a2 + a5) / pi);
			r2 <- sqrt((a2 + a5) / pi);
			r3 <- sqrt(a5 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
		}
		else {
			r1 <- 0.4;
			r2 <- 0.25;
			r3 <- 0.1;
			}
		
		x.centre.1 <- 0.5;
		y.centre.1 <- 0.5;
		x.centre.2 <- 0.5 - offset * (r1 - r2);
		y.centre.2 <- 0.5;
		x.centre.3 <- 0.5 - offset * (r1 - r3);
		y.centre.3 <- 0.5;
		
		a1.x.pos <- (1 + r2 - offset * (r1 - r2) + r1) / 2;
		a1.y.pos <- 0.5;
		a2.x.pos <- (1 + r3 - offset * (r1 - r3) + r2 - offset * (r1 - r2)) / 2;
		a2.y.pos <- 0.5;
		a5.x.pos <- x.centre.3;
		a5.y.pos <- 0.5;
		
	} else if (sp.case.name == "022AAOO"){
		if (scaled) {
			if (a2 >= a7) { d <- find.dist((a2 + a5), (a7 + a5), a5) }
			else  { d <- find.dist((a7 + a5), (a2 + a5), a5) }
			r1 <- sqrt((a2 + a5) / pi);
			r2 <- sqrt((a2 + a5) / pi);
			r3 <- sqrt((a5 + a7) / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			d <- d * shrink.factor;
			}
		else {
			r1 <- 0.2;
			r2 <- 0.2;
			r3 <- 0.2;
			d <- 0.2;
			}

		x.centre.1 <- (1 + r1 - r2 - d) / 2;
		y.centre.1 <- 0.5;
		x.centre.2 <- x.centre.1;
		y.centre.2 <- 0.5;
		x.centre.3 <- x.centre.1 + d;
		y.centre.3 <- 0.5;
		
		a2.x.pos <- (x.centre.1 + x.centre.3 - r1 - r3) / 2;
		a2.y.pos <- 0.5;
		a5.x.pos <- (x.centre.1 + x.centre.3 + r1 - r3) / 2;
		a5.y.pos <- 0.5;
		a7.x.pos <- (x.centre.1 + x.centre.3 + r1 + r3) / 2;
		a7.y.pos <- 0.5;
		
	} else if (sp.case.name == "023"){
		if (scaled) {
			r1 <- sqrt((a2 + a5) / pi);
			r2 <- sqrt(a5 / pi);
			r3 <- sqrt((a2 + a5) / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			}
		else {
			r1 <- 0.4;
			r2 <- 0.2;
			r3 <- 0.4;
			}

		x.centre.1 <- 0.5;
		y.centre.1 <- 0.5;
		x.centre.2 <- 0.5 - offset * (r1 - r2);
		y.centre.2 <- 0.5;
		x.centre.3 <- 0.5;
		y.centre.3 <- 0.5;

		a2.x.pos <- (x.centre.1 + x.centre.2 + r1 + r2) / 2;
		a2.y.pos <- 0.5;
		a5.x.pos <- x.centre.2;
		a5.y.pos <- 0.5;
		
	} else if (sp.case.name == "032"){
		if (scaled) {
			r1 <- sqrt((a1 + a5) / pi);
			r2 <- sqrt(a5 / pi);
			r3 <- sqrt(a5 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			}
		else {
			r1 <- 0.4;
			r2 <- 0.2;
			r3 <- 0.2;
			}
		
		x.centre.1 <- 0.5;
		y.centre.1 <- 0.5;
		x.centre.2 <- 0.5 - offset * (r1 - r2);
		y.centre.2 <- 0.5;
		x.centre.3 <- 0.5 - offset * (r1 - r3);
		y.centre.3 <- 0.5;

		a1.x.pos <- (x.centre.1 + x.centre.2 + r1 + r2) / 2;
		a1.y.pos <- 0.5;
		a5.x.pos <- x.centre.2;
		a5.y.pos <- 0.5;
		
	} else if (sp.case.name == "120"){
		if (scaled) {
			if (a1 >= a3) {	d <- find.dist(a1 + a2, a3 + a2, a2); }
			if (a1 < a3)  { d <- find.dist(a3 + a2, a1 + a2, a2); }
			r1 <- sqrt((a1 + a2) / pi);
			r2 <- sqrt((a3 + a2) / pi);
			r3 <- sqrt(a7 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			d <- d * shrink.factor;
			}
		else {
			r1 <- 0.2;
			r2 <- 0.2;
			r3 <- 0.2;
			d <- 0.2;
			}

		upper.y <- 0.66;
		lower.x <- 0.5;

		x.centre.1 <- (1 + r1 - r2 - d) / 2;
		x.centre.2 <- x.centre.1 + d;
		y.centre.1 <- upper.y;
		y.centre.2 <- upper.y;
		x.centre.3 <- lower.x;

		if (scaled) {
			if (a1 >= a3) {
				y.centre.3 <- y.centre.1 - sqrt(((r1 + r3) * (1 + sep.dist))^ 2 - (x.centre.1 - x.centre.3) ^2);
				}
			if (a1 < a3) {
				y.centre.3 <- y.centre.2 - sqrt(((r2 + r3) * (1 + sep.dist)) ^ 2 - (x.centre.2 - x.centre.3) ^2);
				}
			}
		else {
			if (a1 >= a3) {
				y.centre.3 <- y.centre.1 - sqrt((r1 + r3 + 0.03) ^ 2 - (x.centre.1 - x.centre.3) ^2);
				}
			if (a1 < a3) {
				y.centre.3 <- y.centre.2 - sqrt((r2 + r3 + 0.03) ^ 2 - (x.centre.2 - x.centre.3) ^2);
				}
			}

		a1.x.pos <- (x.centre.1 + x.centre.2 - r1 - r2) / 2;
		a1.y.pos <- upper.y;
		a3.x.pos <- (x.centre.1 + x.centre.2 + r1 + r2) / 2;
		a3.y.pos <- upper.y;
		a2.x.pos <- (x.centre.1 + x.centre.2 + r1 - r2) / 2;
		a2.y.pos <- upper.y;
		a7.x.pos <- x.centre.3;
		a7.y.pos <- y.centre.3;
		
	} else if (sp.case.name == "121AO"){
		if (scaled) {
			r1 <- sqrt((a1 + a2) / pi);
			r2 <- sqrt(a2 / pi);
			r3 <- sqrt(a7 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			}
		else {
			r1 <- 0.2;
			r2 <- 0.1;
			r3 <- 0.2;
			}

		x.centre.1 <- r1;
		y.centre.1 <- 0.5;
		x.centre.2 <- x.centre.1 - offset * (r1 - r2);
		y.centre.2 <- 0.5;
		x.centre.3 <- x.centre.1 + (1 + sep.dist) * (r1 + r3);
		y.centre.3 <- 0.5;

		a1.x.pos <- ((x.centre.1 - r1) + (x.centre.2 - r2)) / 2;
		a1.y.pos <- 0.5;
		a2.x.pos <- x.centre.1;
		a2.y.pos <- 0.5;
		a7.x.pos <- x.centre.3;
		a7.y.pos <- 0.5;
		
	} else if (sp.case.name == "122AAOO"){
		if (scaled) {
			r1 <- sqrt(a2 / pi);
			r2 <- sqrt(a2 / pi);
			r3 <- sqrt(a7 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			}
		else {
			r1 <- 0.2;
			r2 <- 0.2;
			r3 <- 0.2;
			}
		
		x.centre.1 <- r1;
		y.centre.1 <- 0.5;
		x.centre.2 <- r1;
		y.centre.2 <- 0.5;
		x.centre.3 <- r1 + (1 + sep.dist) * (r1 + r3);
		y.centre.3 <- 0.5;

		a2.x.pos <- x.centre.1;
		a2.y.pos <- 0.5;
		a7.x.pos <- x.centre.3;
		a7.y.pos <- 0.5;
		
	} else if (sp.case.name == "130"){
		if (scaled) {
			r1 <- sqrt(a1 / pi);
			r2 <- sqrt(a3 / pi);
			r3 <- sqrt(a7 / pi);
			shrink.factor <- 0.2 / max(r1, r2, r3);
			r1 <- r1 * shrink.factor;
			r2 <- r2 * shrink.factor;
			r3 <- r3 * shrink.factor;
			}
		else {
			r1 <- 0.18;
			r2 <- 0.18;
			r3 <- 0.18;
			}

		a <- (r1 + r2) * (1 + sep.dist[1]);
		b <- (r2 + r3) * (1 + sep.dist[2]);
		c <- (r1 + r3) * (1 + sep.dist[3]);

		beta <- (a^2 + c^2 - b^2) / (2 * a * c);
		gamma <- sqrt(1 - beta^2);
		x.centre.1 <- (r1 - r2 - a + 1) / 2;
		x.centre.3 <- x.centre.1 + c * beta;
		y.centre.3 <- (r3 - r1 + 1 - c * gamma) / 2;
		y.centre.1 <- y.centre.3 + c * gamma;
		x.centre.2 <- x.centre.1 + a;
		y.centre.2 <- y.centre.1;

		a1.x.pos <- x.centre.1;
		a1.y.pos <- y.centre.1;
		a3.x.pos <- x.centre.2;
		a3.y.pos <- y.centre.2;
		a7.x.pos <- x.centre.3;
		a7.y.pos <- y.centre.3;
		
	} else {
		flog.info(paste0("The special case is not in the scaled cases: ",sp.case.name),name="VennDiagramLogger");
	}

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
			area.x = c(a1.x.pos, a2.x.pos, a3.x.pos, a4.x.pos, a5.x.pos, a6.x.pos, a7.x.pos),
			area.y = c(a1.y.pos, a2.y.pos, a3.y.pos, a4.y.pos, a5.y.pos, a6.y.pos, a7.y.pos),
			attach.label.to = area.labels[[sp.case.name]],
			x.centres = c(x.centre.1, x.centre.2, x.centre.3),
			y.centres = c(y.centre.1, y.centre.2, y.centre.3),
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
