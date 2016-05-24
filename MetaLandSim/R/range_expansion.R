range_expansion <-
function(rl, percI, param, b=1, tsteps, iter)
  {
    if (class(rl) != "landscape")
    {
        stop(paste(rl, " should be an object of class class 'landscape'.", sep = ""), 
            call. = FALSE)
    }
    mapsize <- rl$mapsize
    dist_m <- rl$minimum.distance
    areaM <- rl$mean.area
    areaSD <- rl$SD.area
    Npatch <- rl$number.patches
    disp <- rl$dispersal
node.expansion <- function(occ_landscape, param, b, node, tsteps) {
		output0 <- c()
        npatch <- occ_landscape$number.patches
        mapsize <- occ_landscape$mapsize
        ID_land <- max(occ_landscape$nodes.characteristics$ID)
        nrow_land <- nrow(occ_landscape$nodes.characteristics)
		areaM2 <- occ_landscape$mean.area
		areaSD2 <- occ_landscape$SD.area
		dispersal <- 1/param[1,1]
	    ocupp <- "N"
		for(j in 1:tsteps) {

     		ocupp <- "N"
			
			occ_landscape_new <- spom(sp = occ_landscape, kern = "op1", conn = "op1", colnz = "op1", ext = "op1", param_df = param, b, 
                c1 = NULL, c2 = NULL, z = NULL, R = NULL, succ="none")
				
			if(sum(occ_landscape_new$nodes.characteristics$species2)==0) {
				message(paste("Empty landscape. Next simulation... ","- time step ", j,sep=""))
				break		
		}
			
			perc_occup <- (sum(occ_landscape_new$nodes.characteristics$species2[1:nrow_land])*100)/nrow_land
			
			if(node == "North") v0 <- occ_landscape_new$nodes.characteristics[ which(occ_landscape_new$nodes.characteristics$y > (mapsize-dispersal) ), ]$species2
			if(node == "South") v0 <- occ_landscape_new$nodes.characteristics[ which(occ_landscape_new$nodes.characteristics$y < dispersal ), ]$species2
			if(node == "East") v0 <- occ_landscape_new$nodes.characteristics[ which(occ_landscape_new$nodes.characteristics$x > (mapsize-dispersal) ), ]$species2
			if(node == "West") v0 <- occ_landscape_new$nodes.characteristics[ which(occ_landscape_new$nodes.characteristics$x < dispersal ), ]$species2
			
			if(sum(v0)==0) ocupp <- "N"
			if(sum(v0)!=0) ocupp <- "Y"
			
			 if(ocupp == "N") {
		     message(paste("No transition between landscape units. ","- time step ", j,sep=""))
			 occ_landscape <- occ_landscape_new
             occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[, -c(9, 11)]
             names(occ_landscape$nodes.characteristics)[names(occ_landscape$nodes.characteristics) == "species2"] <- "species"
             occ_landscape <- occ_landscape[-9]
			 occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-(nrow_land + 1), ]
			 occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-(nrow_land + 1), ]
	    	 occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, -(nrow_land + 1)]
			 occ_landscape$number.patches <- occ_landscape$number.patches - 1
		     class(occ_landscape) <- "metapopulation"
		     }
     		if(ocupp == "Y"){
			message(paste("Transition between landscapes. New landscape created. ","- time step ", j,sep=""))
			
			if(node == "North") {
			v2 <- length(occ_landscape$nodes.characteristics[ which(occ_landscape$nodes.characteristics$y < dispersal ), ]$ID)
			entry <- "S"
			}
			
			if(node == "South") {
			v2 <- length(occ_landscape$nodes.characteristics[ which(occ_landscape$nodes.characteristics$y > (mapsize-dispersal)), ]$ID)
			entry <- "N"
			}
			
			if(node == "East") {
			v2 <- length(occ_landscape$nodes.characteristics[ which(occ_landscape$nodes.characteristics$x < dispersal ), ]$ID)
			entry <- "W"
			}
			
			if(node == "West") {
			v2 <- length(occ_landscape$nodes.characteristics[ which(occ_landscape$nodes.characteristics$x  > (mapsize-dispersal)), ]$ID)
			entry <- "E"
			}
			
			number_patches <- v2*(nrow_land/100)
			
     		rl1 <- rland.graph(mapsize = occ_landscape$mapsize, dist_m = occ_landscape$minimum.distance, 
                  areaM2, areaSD2, Npatch = npatch, disp = occ_landscape$dispersal, plotG = FALSE)
				  
			occ_landscape <- suppressWarnings(species.graph(rl = rl1, method = "number", parm = number_patches, nsew = entry, plotG = FALSE))

			message(paste("Saving in output file! ","- time step ", j,sep=""))
			output0 <- c(output0,j)
			
			}
			
}

if(sum(output0)==0) output0 <- 0
output1 <- cbind(mapsize,output0)
output1 <- cbind(cumsum(output1[,1]),output1[,2])
colnames(output1) <- c("DISTANCE", "TIME STEP")
output1 <- as.data.frame(output1)
return(output1)
}
    distance <- mapsize * 1:tsteps
    outputN <- as.data.frame(distance)
    outputS <- as.data.frame(distance)
    outputE <- as.data.frame(distance)
    outputW <- as.data.frame(distance)
    for (i in 1:iter) {
        sp1 <- species.graph(rl = rl, method = "percentage", parm = percI, nsew = "none", plotG = FALSE)
        nodeN <- node.expansion(occ_landscape = sp1, param, b, node = "North", tsteps)
        outputN <- suppressWarnings(cbind(outputN, c(nodeN[, 2], rep(NA, nrow(outputN) - length(nodeN[, 2])))))
		message(paste("###### North sub-model concluded for iteration ", i,sep=""))
        nodeS <- node.expansion(occ_landscape = sp1, param, b, node = "South", tsteps)
        outputS <- suppressWarnings(cbind(outputS, c(nodeS[, 2], rep(NA, nrow(outputS) - length(nodeS[, 2])))))
		message(paste("###### South sub-model concluded for iteration ", i,sep=""))
        nodeE <- node.expansion(occ_landscape = sp1, param, b, node = "East", tsteps)
        outputE <- suppressWarnings(cbind(outputE, c(nodeE[, 2], rep(NA, nrow(outputE) - length(nodeE[, 2])))))
		message(paste("###### East sub-model concluded for iteration ", i,sep=""))
        nodeW <- node.expansion(occ_landscape = sp1, param, b, node = "West", tsteps)
        outputW <- suppressWarnings(cbind(outputW, c(nodeW[, 2], rep(NA, nrow(outputW) - length(nodeW[, 2])))))
		message(paste("###### West sub-model concluded for iteration ", i,sep=""))
		message(paste("################### Completed iteration ",i,"! ###################",sep=""))
    }
    outputN[is.na(outputN)] <- 0
    outputS[is.na(outputS)] <- 0
    outputE[is.na(outputE)] <- 0
    outputW[is.na(outputW)] <- 0
	for (x in 2:(iter + 1)) {
        i_N <- which(outputN[, x] != 0)
        outputN[i_N, x] <- 1
        i_S <- which(outputS[, x] != 0)
        outputS[i_S, x] <- 1
        i_E <- which(outputE[, x] != 0)
        outputE[i_E, x] <- 1
        i_W <- which(outputW[, x] != 0)
        outputW[i_W, x] <- 1
    }
	outputN <- cbind(outputN[, 1], rowSums(as.data.frame(outputN[, 2:(ncol(outputN))])))
    outputS <- cbind(outputS[, 1], rowSums(as.data.frame(outputS[, 2:(ncol(outputS))])))
    outputE <- cbind(outputE[, 1], rowSums(as.data.frame(outputE[, 2:(ncol(outputE))])))
    outputW <- cbind(outputW[, 1], rowSums(as.data.frame(outputW[, 2:(ncol(outputW))])))
	outputN <- cbind(outputN[, 1:2], outputN[, 2]/iter)
    outputS <- cbind(outputS[, 1:2], outputS[, 2]/iter)
    outputE <- cbind(outputE[, 1:2], outputE[, 2]/iter)
    outputW <- cbind(outputW[, 1:2], outputW[, 2]/iter)
    outputN <- as.data.frame(outputN)
    outputS <- as.data.frame(outputS)
    outputE <- as.data.frame(outputE)
    outputW <- as.data.frame(outputW)
	names(outputN) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputS) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputE) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputW) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
	message("Preparing graphic output!")
	
    p1 <- gvisLineChart(outputN, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the North", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#006400'}]", 
        backgroundColor = "#D1EEEE"))
    p2 <- gvisLineChart(outputS, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the South", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#0000FF'}]", 
        backgroundColor = "#D1EEEE"))
    p3 <- gvisLineChart(outputE, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the East", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#8B0000'}]", 
        backgroundColor = "#D1EEEE"))
    p4 <- gvisLineChart(outputW, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the West", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#FF4500'}]", 
        backgroundColor = "#D1EEEE"))
    ln1 <- gvisMerge(p1, p2, horizontal = TRUE)
    ln2 <- gvisMerge(p3, p4, horizontal = TRUE)
    ln.final <- gvisMerge(ln1, ln2, horizontal = FALSE)
    ln.final$html$caption <- paste("<div><span>Range expansion in the four cardinal directions, considering that ", 
        percI, "% of the patches are occupied in the first landscape mosaic.</span><br />", 
        sep = "")
    ln.final$html$footer <- paste("\n<!-- htmlFooter -->\n<span> \n  ",R.Version()$version.string,"&#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-", packageVersion("googleVis"),"</a>\n  &#8226; MetaLandSim-",packageVersion("MetaLandSim"),"\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n", sep="")
    output <- list(NORTH = outputN, SOUTH = outputS, EAST = outputE, WEST = outputW)
    plot(ln.final)
    class(output) <- "expansion"
    return(output)
}