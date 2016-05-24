manage_expansion_sim <- function(mapsize, dist_m, areaM, areaSD, Npatch,percI, param, b=1, tsteps, iter, 
								 variable,var_min,var_max,by)
	{

range_expansion1 <-
function(rl, percI, param, b, tsteps, iter)
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
        output0 <- as.data.frame(matrix(nrow = tsteps, ncol = 2))  
        output0[, 1] <- 1:nrow(output0)
        npatch <- occ_landscape$number.patches
        mapsize <- occ_landscape$mapsize
        ID_land <- max(occ_landscape$nodes.characteristics$ID)
        mapsize <- occ_landscape$mapsize
        nrow_land <- nrow(occ_landscape$nodes.characteristics)
        clock <- 1
        if (node == "North") {
            xx1xx <- "S"
            xx2xx <- NA
            xx3xx <- mapsize
            xx4xx <- "North"
            xx7xx <- NA  
            xx8xx <- 0  
            xx9xx <- "South"  
        }
        if (node == "South") {
            xx1xx <- "N"
            xx2xx <- NA
            xx3xx <- 0
            xx4xx <- "South"
            xx7xx <- NA  
            xx8xx <- mapsize  
            xx9xx <- "North"  
        }
        if (node == "East") {
            xx1xx <- "W"
            xx2xx <- mapsize
            xx3xx <- NA
            xx4xx <- "East"
            xx7xx <- 0  
            xx8xx <- NA  
            xx9xx <- "West"  
        }
        if (node == "West") {
            xx1xx <- "E"
            xx2xx <- 0
            xx3xx <- NA
            xx4xx <- "West"
            xx7xx <- mapsize  
            xx8xx <- NA  
            xx9xx <- "East"  
        }
        repeat {
            if (node == "North") {
                xx5xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  2])[abs(mapsize - occ_landscape$nodes.characteristics[, 2]) != 
                  0])
                xx6xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 2])
                xx10xx <- min(occ_landscape$nodes.characteristics[, 2][occ_landscape$nodes.characteristics[, 
                  2] != 0])  
                xx11xx <- occ_landscape$nodes.characteristics[, 2]  
            }
            if (node == "South") {
                xx5xx <- min(occ_landscape$nodes.characteristics[, 2][occ_landscape$nodes.characteristics[, 
                  2] != 0])
                xx6xx <- occ_landscape$nodes.characteristics[, 2]
                xx10xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  2])[abs(mapsize - occ_landscape$nodes.characteristics[, 2]) != 
                  0])  
                xx11xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 2])  
            }
            if (node == "East") {
                xx5xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  1])[abs(mapsize - occ_landscape$nodes.characteristics[, 1]) != 
                  0])
                xx6xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 1])
                xx10xx <- min(occ_landscape$nodes.characteristics[, 1][occ_landscape$nodes.characteristics[, 
                  1] != 0])  
                xx11xx <- occ_landscape$nodes.characteristics[, 1]  
            }
            if (node == "West") {
                xx5xx <- min(occ_landscape$nodes.characteristics[, 1][occ_landscape$nodes.characteristics[, 
                  1] != 0])
                xx6xx <- occ_landscape$nodes.characteristics[, 1]
                xx10xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  1])[abs(mapsize - occ_landscape$nodes.characteristics[, 1]) != 
                  0])  
                xx11xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 1])  
            }
            ocupp <- "N"
            if (clock == 1 || (clock >= 2 & ocupp == "N")) {
                occ_landscape$nodes.characteristics[nrow_land + 1, ] <- NA
                occ_landscape$nodes.characteristics$ID[(nrow_land + 1)] <- (ID_land + 
                  1)
                area1 <- mean(occ_landscape$nodes.characteristics$areas[1:nrow_land])
                occ_landscape$nodes.characteristics$radius[(nrow_land + 1)] <- sqrt((area1 * 
                  10000)/pi)
                occ_landscape$nodes.characteristics$areas[(nrow_land + 1)] <- (occ_landscape$nodes.characteristics$radius[(nrow_land + 
                  1)]) * mapsize
                occ_landscape$nodes.characteristics$y[nrow_land + 1] <- xx3xx
                occ_landscape$nodes.characteristics$x[nrow_land + 1] <- xx2xx
                max_cluster <- max(occ_landscape$nodes.characteristics$cluster, na.rm = T)
                occ_landscape$nodes.characteristics$cluster[(nrow_land + 1)] <- c(max_cluster + 
                  1)
                occ_landscape$nodes.characteristics$colour <- as.character(occ_landscape$nodes.characteristics$colour)
                occ_landscape$nodes.characteristics$colour[(nrow_land + 1)] <- xx4xx
                occ_landscape$nodes.characteristics$species[(nrow_land + 1)] <- 0
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 1)] <- xx5xx
                occ_landscape$number.patches <- occ_landscape$number.patches + 1
                dist_nodos <- occ_landscape$distance.to.neighbours
                dist_nodos[, (nrow_land + 1)] <- c(xx6xx)
                colnames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                dist_nodos[(nrow_land + 1), ] <- c(xx6xx, 0)
                rownames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                occ_landscape$distance.to.neighbours <- dist_nodos
            }
            if (clock >= 2 & ocupp == "Y") {
                rl1 <- rland.graph(mapsize = occ_landscape$mapsize, dist_m = occ_landscape$minimum.distance, 
                  areaM = occ_landscape$mean.area, areaSD = occ_landscape$SD.area, 
                  Npatch = npatch, disp = occ_landscape$dispersal, plotG = FALSE)
                occ_landscape <- suppressWarnings(species.graph(rl = rl1, method = "percentage", 
                  parm = 0, plotG = FALSE))
                occ_landscape$nodes.characteristics[nrow_land + 2, ] <- NA
                occ_landscape$nodes.characteristics$ID[(nrow_land + 1):(nrow_land + 
                  2)] <- (ID_land + 1):(ID_land + 2)
                area1 <- mean(occ_landscape$nodes.characteristics$areas[1:nrow_land])
                occ_landscape$nodes.characteristics$radius[(nrow_land + 1):(nrow_land + 
                  2)] <- sqrt((area1 * 10000)/pi)
                occ_landscape$nodes.characteristics$areas[(nrow_land + 1):(nrow_land + 
                  2)] <- (occ_landscape$nodes.characteristics$radius[(nrow_land + 
                  1)]) * mapsize
                occ_landscape$nodes.characteristics$y[nrow_land + 1] <- xx3xx
                occ_landscape$nodes.characteristics$x[nrow_land + 1] <- xx2xx
                occ_landscape$nodes.characteristics$y[nrow_land + 2] <- xx8xx
                occ_landscape$nodes.characteristics$x[nrow_land + 2] <- xx7xx
                max_cluster <- max(occ_landscape$nodes.characteristics$cluster, na.rm = T)
                occ_landscape$nodes.characteristics$cluster[(nrow_land + 1):(nrow_land + 
                  2)] <- c(max_cluster + 1)
                occ_landscape$nodes.characteristics$colour <- as.character(occ_landscape$nodes.characteristics$colour)
                occ_landscape$nodes.characteristics$colour[(nrow_land + 1)] <- xx4xx
                occ_landscape$nodes.characteristics$colour[(nrow_land + 2)] <- xx9xx
                occ_landscape$nodes.characteristics$species[(nrow_land + 1)] <- 0
                occ_landscape$nodes.characteristics$species[(nrow_land + 2)] <- 1
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 1)] <- xx5xx
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 2)] <- xx10xx
                occ_landscape$number.patches <- occ_landscape$number.patches + 2
                dist_nodos <- occ_landscape$distance.to.neighbours
                dist_nodos[, (nrow_land + 1)] <- NA
                dist_nodos[, (nrow_land + 2)] <- NA
                dist_nodos[(nrow_land + 1), ] <- NA
                dist_nodos[(nrow_land + 2), ] <- NA
                dist_nodos[, (nrow_land + 1)] <- c(xx6xx, 0, mapsize)
                dist_nodos[, (nrow_land + 2)] <- c(xx11xx, mapsize, 0)
                dist_nodos[(nrow_land + 1), ] <- c(xx6xx, 0, mapsize)
                dist_nodos[(nrow_land + 2), ] <- c(xx11xx, mapsize, 0)
                colnames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                colnames(dist_nodos)[(nrow_land + 2)] <- xx9xx
                rownames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                rownames(dist_nodos)[(nrow_land + 2)] <- xx9xx
                occ_landscape$distance.to.neighbours <- dist_nodos
            }
            class(occ_landscape) <- "metapopulation"
			occ_landscape_new <- spom(sp = occ_landscape, kern = "op1", 
                conn = "op1", colnz = "op1", ext = "op1", param_df = param, b = b, 
                c1 = NULL, c2 = NULL, z = NULL, R = NULL)
            if (occ_landscape_new$nodes.characteristics$species2[nrow_land + 1] == 
                0) 
                ocupp <- "N"
            if (occ_landscape_new$nodes.characteristics$species2[nrow_land + 1] == 
                1) 
                ocupp <- "Y"
            if (ocupp == "N") {
                occ_landscape <- occ_landscape_new
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[, 
                  -c(9, 11)]
                names(occ_landscape$nodes.characteristics)[names(occ_landscape$nodes.characteristics) == 
                  "species2"] <- "species"
                occ_landscape <- occ_landscape[-9]
            }
            if (ocupp == "Y") 
                value <- clock
            if (ocupp == "N") 
                value <- 0
            output0[clock, 2] <- value
            if (clock == tsteps) 
                break
            if (nrow(occ_landscape$nodes.characteristics) == (npatch + 1)) {
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-(nrow_land + 
                  1), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-(nrow_land + 
                  1), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, 
                  -(nrow_land + 1)]
                occ_landscape$number.patches <- occ_landscape$number.patches - 1
            }
            if (nrow(occ_landscape$nodes.characteristics) == (npatch + 2)) {
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-c((nrow_land + 
                  1):(nrow_land + 2)), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-c((nrow_land + 
                  1):(nrow_land + 2)), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, 
                  -c((nrow_land + 1):(nrow_land + 2))]
                occ_landscape$number.patches <- occ_landscape$number.patches - 2
            }
            clock <- clock + 1
        }
        if (sum(output0[, 2]) == 0) {
            out_vec <- 0  
            output1 <- cbind(0, out_vec)
            colnames(output1) <- c("DISTANCE", "TIME STEP")
            output1 <- as.data.frame(output1)
			return(output1)

        }
        if (sum(output0[, 2]) != 0) {
            output1 <- output0[output0[, 2] != 0, ]  
            output1[, 1] <- output1[, 1] * mapsize
            colnames(output1) <- c("DISTANCE", "TIME STEP")
            output1 <- as.data.frame(output1)
            return(output1)
        }
    } 
    distance <- mapsize * 1:tsteps
    outputN <- distance
    outputS <- distance
    outputE <- distance
    outputW <- distance
    for (i in 1:iter) {
        sp1 <- species.graph(rl = rl, method = "percentage", parm = percI, nsew = "none", 
            plotG = FALSE)
        nodeN <- node.expansion(occ_landscape = sp1, param, b, node = "North", 
            tsteps)
        outputN <- suppressWarnings(cbind(outputN, c(nodeN[, 2], rep(NA, length(outputN) - 
            length(nodeN[, 2])))))
        nodeS <- node.expansion(occ_landscape = sp1, param, b, node = "South", 
            tsteps)
        outputS <- suppressWarnings(cbind(outputS, c(nodeS[, 2], rep(NA, length(outputS) - 
            length(nodeS[, 2])))))
        nodeE <- node.expansion(occ_landscape = sp1, param, b, node = "East", 
            tsteps)
        outputE <- suppressWarnings(cbind(outputE, c(nodeE[, 2], rep(NA, length(outputE) - 
            length(nodeE[, 2])))))
        nodeW <- node.expansion(occ_landscape = sp1, param, b, node = "West", 
            tsteps)
        outputW <- suppressWarnings(cbind(outputW, c(nodeW[, 2], rep(NA, length(outputW) - 
            length(nodeW[, 2])))))
    
	}
    outputN[is.na(outputN)] <- 0
    outputS[is.na(outputS)] <- 0
    outputE[is.na(outputE)] <- 0
    outputW[is.na(outputW)] <- 0
	
	tstepsN <- outputN[,2]
	tstepsS <- outputS[,2]
	tstepsE <- outputE[,2]
	tstepsW <- outputW[,2]
	
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
    outputN <- cbind(outputN[, 1], rowSums(outputN[, 2:(ncol(outputN))]))
    outputS <- cbind(outputS[, 1], rowSums(outputS[, 2:(ncol(outputS))]))
    outputE <- cbind(outputE[, 1], rowSums(outputE[, 2:(ncol(outputE))]))
    outputW <- cbind(outputW[, 1], rowSums(outputW[, 2:(ncol(outputW))]))
    outputN <- cbind(outputN[, 1:2], outputN[, 2]/iter, tstepsN)
    outputS <- cbind(outputS[, 1:2], outputS[, 2]/iter, tstepsS)
    outputE <- cbind(outputE[, 1:2], outputE[, 2]/iter, tstepsE)
    outputW <- cbind(outputW[, 1:2], outputW[, 2]/iter, tstepsW)
    outputN <- as.data.frame(outputN)
    outputS <- as.data.frame(outputS)
    outputE <- as.data.frame(outputE)
    outputW <- as.data.frame(outputW)
    names(outputN) <- c("DISTANCE", "OCCUPATION", "PROPORTION", "TIME STEPS")
    names(outputS) <- c("DISTANCE", "OCCUPATION", "PROPORTION", "TIME STEPS")
    names(outputE) <- c("DISTANCE", "OCCUPATION", "PROPORTION", "TIME STEPS")
    names(outputW) <- c("DISTANCE", "OCCUPATION", "PROPORTION", "TIME STEPS")
    output <- list(NORTH = outputN, SOUTH = outputS, EAST = outputE, WEST = outputW)
    class(output) <- "expansion"
    return(output)
  }
								 
var1 <- rep(seq(from=var_min,to=var_max,by=by)) 

outputN <- data.frame(matrix(nrow = tsteps, ncol = length(var1)))
outputS <- data.frame(matrix(nrow = tsteps, ncol = length(var1)))
outputE <- data.frame(matrix(nrow = tsteps, ncol = length(var1)))
outputW <- data.frame(matrix(nrow = tsteps, ncol = length(var1)))
clock <- 1

namesN <- as.character()
namesS <- as.character()
namesE <- as.character()
namesW <- as.character()

for (i in 1:length(var1)){

if (variable=="dist") dist_m <- var1[i]
if (variable=="area") areaM <- var1[i]
if (variable=="npatch") Npatch <- var1[i]
if (variable=="sizevar") areaSD <- var1[i]

rl <- rland.graph(mapsize, dist_m, areaM, areaSD, Npatch, disp=0, plotG=FALSE)

rg <- range_expansion1(rl, percI, param, b, tsteps, iter)

if(clock==1){
outputN[,clock] <- rg$NORTH$DISTANCE/1000
outputS[,clock] <- rg$SOUTH$DISTANCE/1000
outputE[,clock] <- rg$EAST$DISTANCE/1000
outputW[,clock] <- rg$WEST$DISTANCE/1000
}

outputN[,clock+1] <- rg$NORTH[,4]
outputS[,clock+1] <- rg$SOUTH[,4]
outputE[,clock+1] <- rg$EAST[,4]
outputW[,clock+1] <- rg$WEST[,4]

clock <- clock+1
}

if (variable=="area") varname <- "MEAN PATCH AREA"
if (variable=="dist") varname <- "MINIMUM DISTANCE BETWEEN PATCHES"
if (variable=="npatch") varname <- "NUMBER OF PATCHES"
if (variable=="sizevar") varname <- "STAND. DEV. OF PATCH SIZE"

names(outputN) <- c("DISTANCE(km)",as.character(var1))
names(outputS) <- c("DISTANCE(km)",as.character(var1))
names(outputE) <- c("DISTANCE(km)",as.character(var1))
names(outputW) <- c("DISTANCE(km)",as.character(var1))

	outN <- as.data.frame(matrix(ncol=3,nrow=length(var1)))
	names(outN) <- c(varname,"MEAN EXPANSION SPEED","MAXIMUM EXPANSION DISTANCE")
	outN[,1] <- var1
	for(i in 1:nrow(outN)){
	val1 <- outputN[,i+1]
	val2 <- length(val1[val1!=0])
	if (val2==0) {
	outN[i,2] <- NA
	outN[i,3] <- NA	
	}
	if (val2!=0) {
	val3 <- outputN[,1][val2]
	outN[i,2] <- val3/tsteps
	outN[i,3] <- val3
	}
	}
	
	outS <- as.data.frame(matrix(ncol=3,nrow=length(var1)))
	names(outS) <- c(varname,"MEAN EXPANSION SPEED","MAXIMUM EXPANSION DISTANCE")
	outS[,1] <- var1
	for(i in 1:nrow(outS)){
	val1 <- outputN[,i+1]
	val2 <- length(val1[val1!=0])
	if (val2==0) {
	outS[i,2] <- NA
	outS[i,3] <- NA	
	}
	if (val2!=0) {
	val3 <- outputN[,1][val2]
	outS[i,2] <- val3/tsteps
	outS[i,3] <- val3
	}
	}
	
	outE <- as.data.frame(matrix(ncol=3,nrow=length(var1)))
	names(outE) <- c(varname,"MEAN EXPANSION SPEED","MAXIMUM EXPANSION DISTANCE")
	outE[,1] <- var1
	for(i in 1:nrow(outE)){
	val1 <- outputN[,i+1]
	val2 <- length(val1[val1!=0])
	if (val2==0) {
	outE[i,2] <- NA
	outE[i,3] <- NA	
	}
	if (val2!=0) {
	val3 <- outputN[,1][val2]
	outE[i,2] <- val3/tsteps
	outE[i,3] <- val3
	}
	}	

	outW <- as.data.frame(matrix(ncol=3,nrow=length(var1)))
	names(outW) <- c(varname,"MEAN EXPANSION SPEED","MAXIMUM EXPANSION DISTANCE")
	outW[,1] <- var1
	for(i in 1:nrow(outW)){
	val1 <- outputN[,i+1]
	val2 <- length(val1[val1!=0])
	if (val2==0) {
	outW[i,2] <- NA
	outW[i,3] <- NA	
	}
	if (val2!=0) {
	val3 <- outputN[,1][val2]
	outW[i,2] <- val3/tsteps
	outW[i,3] <- val3
	}
	}


	pN2 <- gvisLineChart(outN, xvar = varname, yvar = "MAXIMUM EXPANSION DISTANCE", options = list(title = "North-Distance", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'Maximum Expansion Distance(km)'}", hAxis = "{title: 'Variable'}", series = "[{color: '#006400'}]", 
        backgroundColor = "#D1EEEE"))
		
	pS2 <- gvisLineChart(outS, xvar = varname, yvar = "MAXIMUM EXPANSION DISTANCE", options = list(title = "South-Distance", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'Maximum Expansion Distance(km)'}", hAxis = "{title: 'Variable'}", series = "[{color: '#0000FF'}]", 
        backgroundColor = "#D1EEEE"))
		
	pE2 <- gvisLineChart(outE, xvar = varname, yvar = "MAXIMUM EXPANSION DISTANCE", options = list(title = "East-Distance", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'Maximum Expansion Distance(km)'}", hAxis = "{title: 'Variable'}", series = "[{color: '#8B0000'}]", 
        backgroundColor = "#D1EEEE"))
		
	pW2 <- gvisLineChart(outW, xvar = varname, yvar = "MAXIMUM EXPANSION DISTANCE", options = list(title = "West-Distance", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'Maximum Expansion Distance(km)'}", hAxis = "{title: 'Variable'}", series = "[{color: '#FF4500'}]", 
        backgroundColor = "#D1EEEE"))		
		
	ln1 <- gvisMerge(pN2, pS2, horizontal = TRUE)
    ln2 <- gvisMerge(pE2, pW2, horizontal = TRUE)
    ln.final <- gvisMerge(ln1, ln2, horizontal = FALSE)
		
    ln.final$html$caption <- paste("<div><span> MetaLandSim range expansion simulation - Expansion scenarios produced 
	considering variation in ",varname,". </span><br />", sep = "")
    ln.final$html$footer <- paste("\n<!-- htmlFooter -->\n<span> \n  ",R.Version()$version.string,"&#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-", packageVersion("googleVis"),"</a>\n  &#8226; MetaLandSim-",packageVersion("MetaLandSim"),"\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n", sep="")
    
	output <- list(NORTH = outN, SOUTH = outS, EAST = outE, WEST = outW, SimN = outputN, SimS = outputS, SimE = outputE, SimW = outputW)
	
    plot(ln.final)
	
    return(output)
	}
	