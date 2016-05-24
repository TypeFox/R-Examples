########################################################################
# This section contains a couple of general GIS related functions
########################################################################

# Earth radius: WSG84 equaritorial radius (km)

earth.radius <- 6378.137

# Mathematical function to calculate the central angle of a great circle defined by it's endpoints given as spherical coordinates 
# excluding the radius which is constant by definition. 

central.angle <- function(theta1, phi1, theta2, phi2) {
	2 * asin((sin((theta1 - theta2) / 2)^2 + cos(theta1) * cos(theta2) * sin((phi1 - phi2) / 2)^2)^0.5)
}

# Difference in theta per difference of surface projected y on a spherical surface in degrees / km
dtheta.dy <- function(r) 1 / r * 180 / pi

# Difference in phi per difference of surface projected x on a spherical surface in degrees / km
dphi.dx <- function(r, theta) 1 / (cos(theta * pi / 180) * r) * 180 / pi

##############################
# GIS.Exposure
##########################
# Exposure computes exposure using a given concentration matrix and protected population data from Heande.
# Inputs: 
#	Concentration.matrix - A matrix containing spatially dependent concentratio data; 
#	LO & LA - coordinates of the center of the concentration matrix; 
#	distx & disty - maximum displacement from center of concentration matrix, assumed symmetrical, defaults to 10.5 km (PLTTI matrix);
#	resolution - resolution of concentration matrix, length of side of grid element which are assumed squares, defaults to 1 km (PILTTI matrix)
# Output:
#	An ovariable containing result of Population * Concentration. Output and marginal slots are defined. Spatial information lost in 
# 	summation (though there is an easy way around it). 
#########################################

GIS.Exposure <- function(
		Concentration.matrix, 
		dbug = FALSE,
		...
) {
	bounds <- unique(Concentration.matrix@output[c("LAbin", "LObin")])
	LAlower <- NA
	LAupper <- NA
	LOlower <- NA
	LOupper <- NA
	koord_lower <- list()
	koord_upper <- list()
	a <- 1
	
	Population <- data.frame()
	
	for (i in 1:nrow(bounds)) {
		if (dbug) {
			cat("Internal loop", a, "start time:", as.character(Sys.time()), ".\n")
			a <- a+1
		}
		
		tmp <- strsplit(as.character(bounds$LAbin[i]), ",")[[1]]
		LAlower[i] <- as.numeric(substring(tmp[1], 2, nchar(tmp[1])))
		LAupper[i] <- as.numeric(substring(tmp[2], 1, nchar(tmp[2])-1))
		tmp <- strsplit(as.character(bounds$LObin[i]), ",")[[1]]
		LOlower[i] <- as.numeric(substring(tmp[1], 2, nchar(tmp[1])))
		LOupper[i] <- as.numeric(substring(tmp[2], 1, nchar(tmp[2])-1))
		
		# Multiple database reads --> extremely slow
		
		#koord_lower[[i]] <- koordGT(LAlower[i], LOlower[i])
		#koord_upper[[i]] <- koordGT(LAupper[i], LOupper[i])
		
		# Use error handling in case no data found within bounds
		#pop <- tryCatch(
		#	tidy(
		#		opbase.data(
		#			"Op_en2949", 
		#			subset = "2012",
		#			range = list(
		#				XKOORD = c(koord_lower[[i]]$E, koord_upper[[i]]$E), #XKOORD = c(LOlower, LOupper),
		#				YKOORD = c(koord_lower[[i]]$N, koord_upper[[i]]$N)#YKOORD = c(LAlower, LAupper)
		#			)
		#		)
		#	), 
		#	error = function(...) return(NULL)
		#)
		#if (!is.null(pop)) {
		#	pop <- merge(bounds[i,], pop)
		#	
		#	Population <- rbind(Population, pop)
		#}
	}
	
	pop_format <- CRS("+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs")
	longlat_format <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
	
	koord_lower <- SpatialPoints(cbind(LOlower, LAlower), longlat_format)
	koord_upper <- SpatialPoints(cbind(LOupper, LAupper), longlat_format)
	
	koord_lower <- spTransform(koord_lower, pop_format)
	koord_upper <- spTransform(koord_upper, pop_format)
	
	XKOORD <- c(min(coordinates(koord_lower)[,1]), max(coordinates(koord_upper)[,1]))
	YKOORD <- c(min(coordinates(koord_lower)[,2]), max(coordinates(koord_upper)[,2]))
	
	if (dbug) cat(XKOORD, YKOORD, "\n")
	
	Population <- tryCatch(
		tidy(
			opbase.data(
				"Op_en2949", 
				subset = "2012",
				range = list(
					XKOORD = XKOORD,
					YKOORD = YKOORD
				)#,
				#...
			)
		), 
		error = function(...) return(NULL)
	)
	
	if (is.null(Population)) stop("Data download failed!")
	if (nrow(Population) == 0) stop("No population data at these coordinates.")
	if (dbug) cat(nrow(Population), "\n")
	
	Population$LObin <- NA
	Population$LAbin <- NA
	
	for (i in 1:nrow(bounds)) {
		cond <- (
			Population$XKOORD >=  coordinates(koord_lower)[i,1] & 
			Population$XKOORD < coordinates(koord_upper)[i,1] &
			Population$YKOORD >=  coordinates(koord_lower)[i,2] & 
			Population$YKOORD < coordinates(koord_upper)[i,2]
		)
		if (dbug) {
			cat("Bound", i, "matching rows:", sum(cond), "\n")
		}
		Population[cond, "LObin"] <- as.character(bounds[i, "LObin"])
		Population[cond, "LAbin"] <- as.character(bounds[i, "LAbin"])
	}
	# Remove rows that do not fall into bins. (The fact this happens is indicative of some problems)
	Population <- Population[!is.na(Population$LObin),] 
	colnames(Population)[colnames(Population)=="Result"] <- "PopulationResult"
	
	Population$PopulationResult <- ifelse(Population$PopulationResult < 0, 0, Population$PopulationResult)
	
	Population <- Ovariable("Population", output = Population, marginal = colnames(Population) %in% c("Iter", "LObin", "LAbin", "HAVAINTO"))
	
	if(dbug) {
		cat(colnames(Concentration.matrix@output), "\n")
		cat(colnames(Population@output), "\n")
	}
	
	# Calculating exposure. 
	out <- Population * Concentration.matrix
	
	return(out)
}

######################################################
# GIS.Concentration.matrix 
##############################################
# Computes a concentration matrix from given emission and coordinates, based on random sampling PILTTI source-receptor-matrices.
# Inputs:
#	Emission - emission of substance in Mga^-1; 
#	LO & LA - coordinates where emission occurs; 
#	distx & disty - maximum displacement in kilometers from center of desired matrix, assumed symmetrical, 
#		defaults to 10.5 km (PLTTI matrix);
#	resolution - resolution of desired matrix (length in kilometers of side of grid element which are assumed squares), 
#		defaults to 1 km (PILTTI matrix);
#	N - number of iterations to be run
# Output:
#	An ovariable containing spatially dependent concentration data. Output and marginal slots are defined. 
##################################

GIS.Concentration.matrix <- function(
		Emission, 
		LO, 
		LA, 
		distx = 10.5, 
		disty = 10.5, 
		resolution = 1, 
		N = 1000, 
		dbug = FALSE, 
		...
) {
	LaPerKm <- dtheta.dy(earth.radius)
	LoPerKm <- dphi.dx(earth.radius, LA)
	
	# PILTTI source-receptor-matrices
	
	PILTTI.matrix <- tidy(op_baseGetData("opasnet_base", "Op_en5797", ...), objname = "PILTTI.matrix") # unit: ugm^-3/Mga^-1
	
	if (N == 0) {
		PILTTI.matrix <- as.data.frame(as.table(tapply(PILTTI.matrix[["PILTTI.matrixResult"]], PILTTI.matrix[,c("dx", "dy")], mean)))
		colnames(PILTTI.matrix)[colnames(PILTTI.matrix) == "Freq"] <- "PILTTI.matrixResult"
		#PILTTI.matrix <- PILTTI.matrix[,c("dy", "dx", "PILTTI.matrixResult")]
	} else {
		# Sampling; first make lists containing row numbers of individual matrices defined in the data. 
		ID.list <- tapply(1:nrow(PILTTI.matrix), PILTTI.matrix[,c("Kaupunki", "Vuosi", "Tyyppi")], list)
		# Then randomly pick N elements of that list to a new list. 
		ID.list.samples <- sample(ID.list, N, replace = TRUE)
		# For which we find the length of each individual list.
		ID.sample.lengths <- sapply(ID.list.samples, length)
		# Take all the values in the list and make one big vector out of it. 
		ID.vec <- unlist(ID.list.samples)
		# Use that vector to select corresponding rows from the original data. 
		PILTTI.matrix <- PILTTI.matrix[ID.vec, c("dx", "dy", "PILTTI.matrixResult")]
		# Add iteration indicator by repeating the numbers 1 to N according to the lengths of the list elements. 
		PILTTI.matrix$Iter <- rep(1:N, times = ID.sample.lengths)
	}
	
	PILTTI.matrix$dy <- as.numeric(as.character(PILTTI.matrix$dy))
	PILTTI.matrix$dx <- as.numeric(as.character(PILTTI.matrix$dx))
	
	# dx and dy in PILTTI matrix is given in meters
	PILTTI.matrix$LObin <- cut(PILTTI.matrix$dx / 1000 * LoPerKm + LO, breaks = LO + seq(-distx, distx, resolution) * LoPerKm)
	PILTTI.matrix$LAbin <- cut(PILTTI.matrix$dy / 1000 * LaPerKm + LA, breaks = LA + seq(-disty, disty, resolution) * LaPerKm)
	
	PILTTI.matrix <- new(
			"ovariable",
			name = "PILTTI.matrix",
			output = PILTTI.matrix,
			marginal = colnames(PILTTI.matrix) %in% c("LObin", "LAbin", "Iter")
	)
	
	if(dbug) {
		cat(colnames(PILTTI.matrix@output), "\n")
		cat(colnames(Emission@output), "\n")
	}
	# Calculate concentratios based on emission and source-receptor-matrix
	out <- PILTTI.matrix * Emission
	return(out)
}

