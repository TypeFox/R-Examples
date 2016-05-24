# date: 08.08.2013

# S4 method

# x = a input inlet node of class SpatialPointsDataframe
# y = a riverIO of class SpatialLinesDataFrame
# z = a nodeIO of class SpatialPointsDataFrame
# zhyd = SpatialPolygonsDataFrame for searching the tributary

setGeneric("Watershed.Tributary", function(x,xo,y,z,zhyd,c1,...) standardGeneric("Watershed.Tributary"))

setMethod("Watershed.Tributary", signature=c("SpatialPointsDataFrame", "SpatialPointsDataFrame",
	"SpatialLinesDataFrame","SpatialPointsDataFrame","SpatialPolygonsDataFrame","SpatialPolygonsDataFrame"), 
	function(x, xo, y, z, zhyd, c1){
	# x = c1_inlet
	# xo =c1_outlet 
	# y = riverIO
	# z = nodeIO
	# zhyd = zhyd_sb1
	# c1 = c1 #first watershed

	# identifying the tributary in c1_inlet
		touch = SpDF_Touch(x, y); touch
		tributaryId = touch[[2]]; tributaryId
		if(length(tributaryId)==2){ # 2 tributary rivers
			id = lapply(tributaryId, function(x){list(slot(y, "data")["OBJECTID"]==x)}); id
			id = list(id[[1]][[1]]|id[[2]][[1]])	
		}else if(length(tributaryId)==3){ # 3 tributary rivers
			id = lapply(tributaryId, function(x){list(slot(y, "data")["OBJECTID"]==x)}); id
			id1 = list(id[[1]][[1]]|id[[2]][[1]]); id1
			id = list(id1[[1]]|id[[3]][[1]]); id
			tributary = SpDF_Subset(id, y); tributary
			noId = RiverStation(xo, tributary); noId # no tributary node
			noId = which(touch[[2]]==slot(noId, "data")["OBJECTID"])
			id = lapply(tributaryId[-noId], function(x){list(slot(y, "data")["OBJECTID"]==x)}); id
			id = list(id[[1]][[1]]|id[[2]][[1]])	
		}
		tributary = SpDF_Subset(id, y)
		# plot(c1)
		# plot(tributary, add=T)
		# plot(z, add=T)
		# plot.PointAttribute(x=z, y="OBJECTID", dist=350, cex=.6)
		# plot.PolyLineAttribute(x=c1, y="OBJECTID", dist=350, cex=.6)
		
	# identifying nodes of tributary
		a = SpDF_Touch(z, tributary)
		id = which(a[[1]][,2]==1)
		node_tribId = a[[1]][id,1]; node_tribId
		id1 = (slot(z, "data")["OBJECTID"]==node_tribId[1]); id1
		id2 = (slot(z, "data")["OBJECTID"]==node_tribId[2]); id2
		id = 	list(id1|id2); id
		node_trib= SpDF_Subset(id, z); node_trib
		
	# identifying tributary watershed
		id = list(gIntersects(node_trib, zhyd, byid=T)); id
		c2c3 = SpDF_Subset(id, zhyd)	# current tributary watershed
		#plot(basin)
		#plot(c1, col="green", add=T)
		#plot(c2c3)
		#plot(c2c3, col="green", add=T)
		#plot(station, col="red", add=T)
		#plot(node_trib, col="red", add=T)
		#plot(y, add=T)
				
	# identifying greater tributary watershed
		areas = slot(c2c3, "data")["Shape_Area"]
		which(areas == max(areas))
		maxAreas = max(areas); maxAreas
		id = list(areas == max(areas)); id
		c2 = SpDF_Subset(id, c2c3)

		if(length(colnames(slot(c1, "data")))==length((slot(c1, "data")=="order"))){
			# does not exist "order" slot		
			slot(c1,"data")["order"]=1
		}
		
		# assigning order to watershed
		# nchar(slot(c1, "data")["order"])
		# nchar("1") # = 1 => 1/10 = 1/1e1 = 1/10^1
		# nchar("1.1") # = 3 => 1/100 = 1/1e2 = 1/10^(3-1)
		# nchar("1.11") # = 4 => 1/1000 = 1/1e3 = 1/10^(4-1)
		# nchar("1.111") # = 5 => 1/10000 = 1/1e4 = 1/10^(5-1)
		if (nchar(slot(c1, "data")["order"]) == 1){
			slot(c2, "data")["order"] = as.character(as.numeric(slot(c1, "data")["order"]) + 1/10^1)
			} else{
				n = nchar(slot(c1, "data")["order"]); n
				slot(c2, "data")["order"] = as.character(as.numeric(slot(c1, "data")["order"]) + 1/10^(n-1))
			}
		

		# identifying the second tributary
			id = list(areas != max(areas)); id
			c3 = SpDF_Subset(id, c2c3)
	
			# assigning order to watershed
				# nchar(slot(c1, "data")["order"])
				# nchar("1") # = 1 => 2/10 = 2/1e1 = 2/10^1
				# nchar("1.1") # = 3 => 2/100 = 2/1e2 = 2/10^(3-1)
				# nchar("1.11") # = 4 => 2/1000 = 2/1e3 = 2/10^(4-1)
				if (nchar(slot(c1, "data")["order"]) == 1){
					slot(c3, "data")["order"] = as.character(as.numeric(slot(c1, "data")["order"]) + 2/10^1)
				}else{
					n = nchar(slot(c1, "data")["order"]); n
					slot(c3, "data")["order"] = as.character(as.numeric(slot(c1, "data")["order"]) + 2/10^(n-1))
				}
			
	
				# plot(sb1)			
				# plot(c2c3)		
				# plot(c1, add=T, col="white")
				# plot(c2, col="green", add=T)
				# plot(c3, col="coral", add=T)
				# plot.PolyLineAttribute(c1, "order", 450, 0.8)
				# plot.PolyLineAttribute(c2, "order", 450, 0.8)
				# plot.PolyLineAttribute(c3, "order", 450, 0.8)				
				# plot(x, add= T)	
				# plot(xo, add= T)	
				# plot(node_trib, col = "brown", add= T)

	return(list(c2c3, c2, c3, node_trib))
}
)