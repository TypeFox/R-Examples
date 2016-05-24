# By: J.A. Torres-Matallana
# 05.07.2012
# Institute for Geoinformatics (ifgi) 
# University of Muenster, Germany

# S4 method

# station
# subbasin
# zhyd
# river
# c1, 0 if does not exist on contrary a SpatialPolygons DataFrame that represents
#		the initial watershed
 
setGeneric("Watershed.Order", function(x,...) standardGeneric("Watershed.Order"))

setMethod("Watershed.Order", signature="Watershed", 
	function(x,...){
	# Watershed.Order = function(station,subbasin,zhyd,river,c1){
	# x = watershed
	
	station = slot(x,"station")
	subbasin = slot(x,"subbasin")
	zhyd = slot(x,"zhyd")
	river = slot(x,"river")
	c1 = slot(x,"c1")
	node = slot(x,"node")
		# identifying subbasin that contains "station"
			id = list(gIntersects(station, subbasin, byid=T))
			sb1 = SpDF_Subset(id, subbasin) # current subbasin
				#class(sb1) #SpPolysDF
				#str(sb1)
				#plot(basin)
				#plot(sb1, add=T)
				#plot(station, col="red", add=T)
			
		# subsetting zhyd that belongs to subbasin
			#class(zhyd) #SpPolysDF
			#class(sb1) #SpPolysDF
			# id =  SpDFInt_Poly(sb1, zhyd)
			# id =  intersect(zhyd, sb1)
			id =  list(gIntersects(sb1, zhyd, byid=T))

			#time1 = proc.time()
			zhyd_sb1 = SpDF_Subset(id, zhyd) # zhyd in subbasin. To optimize
				#time2 = proc.time()
				#time2 -time1
				#zhyd_sb1 = parallel(SpDF_Subset(id, zhyd)) # zhyd in subbasin. To optimize
				#zhyd_sb1 = collect(zhyd_sb1)[[1]]
				#plot(basin)
				#plot(zhyd_sb1, add=T)
				#plot(zhyd)
				#plot(station, col="red", add=T)

	 if(identical(c1, subbasin)){ 	
		# c1 does not exist
		# identifying current watershed that contains "station"
			#class(zhyd_sb1) #SpPolysDF
			#class(station) #SpPoints
			id =  list(gIntersects(zhyd_sb1, station, byid=T))
			c1 = SpDF_Subset(id, zhyd_sb1) # current watershed
				#class(c1) #SpPolysDF
				#str(c1)
				#plot(basin)
				#plot(zhyd_sb1, add=T)
				#plot(c1, col="green", add=T)
				#plot(station, col="red", add=T)
			
		# assigning order to watershed
			slot(c1, "data")["order"] = 1
	}				
	
	# Identifying inlet and outlet stretches
	#class(c1) # SpPolysDF
	#class(river) # SpLinesDF
	# id = SpDFCrosses_Poly(c1, river)
	# plot(c1)
	id = list(gCrosses(c1, river, byid=T))		
		#length(id[[1]])
		#length(river)
	if(length(id[[1]])==length(which(id[[1]]==FALSE))){# there are no riverIO
		# typical coastal watershed
		c1_outlet = 0
		c1_inlet = 0
		riverIO = 0
	}else{
	riverIO = SpDF_Subset(id, river)
		#class(riverIO) # SpLinesDF
		#str(riverIO, level=1)
		#plot(c1)
		#plot(riverIO, col="blue", add=T)
		#plot_tPolyLine(riverIO, "STRAHLER", 250, 0.8)
		#plot_tPolyLine(c1, "OBJECTID", 450, 0.8)
		#plot(station, col="green", add=T)
	
		#length(riverIO)
		# if  length(riverIO) == 1

		# Subsetting "node" that intersects "c1"
			#class(c1) #SpatialPolygonsDataFrame
			#class(node) #SpatialPoinsDataFrame
			id = list(gIntersects(c1, node,  byid=T))
			c1_node = SpDF_Subset(id, node) # all nodes inside watershed
		
		# Subsetting "river" that intersects "c1"
			id = list(gIntersects(c1, river, byid=T))
			c1_river = SpDF_Subset(id, river)
		
		# Identifying inlet and outlet nodes
			#class(riverIO) # SpLinesDF
			#class(node) # SpPointsDF
			# subsetting points of IO stretches by buffering the stretches	
				buffer = gBuffer(riverIO, width=100)
					#class(buffer)
					#plot(buffer)
	
		#id = SpDFInt_Poly(buffer, node)
			id = list(gIntersects(buffer, node, byid=T))
				#class(id)
				#id
				#length(id[[1]])
				#length(node)
		nodeIO = SpDF_Subset(id, node)
			#class(nodeIO)
	
		#plot(c1)
		#plot(nodeIO, col="red", add=T)
		#plot(station, col="green", add=T)
		#plot(buffer, lwd=1, add=T)
		#plot(rWeser, add=T)
		#plot(riverIO, col="blue", add=T)
		#plot_tPoint(nodeIO, "OBJECTID", 1200, .5)
		#plot_tPoint(nodeIO, "ELEV", 350, .8)
	
		# determining inlet and output nodes
			# determining distances of nodeIO to c1
				boundary = gBoundary(c1)
				dist = gDistance(nodeIO, boundary, byid =T); dist		
				dist
	}
	
		if (class(riverIO)=="SpatialLinesDataFrame"){
			if (length(riverIO) == 1){	
				a = Watershed.IOR1(x=nodeIO, dist=dist)
				print("length(riverIO) == 1")
			}	
				
			if	(length(riverIO) == 2){	
				a = Watershed.IOR2(x=nodeIO, dist=dist, node=c1_node)
				print("length(riverIO) == 2")
			}
			if(length(riverIO) == 3){
				a = Watershed.IOR3(x=nodeIO, y=riverIO, dist=dist)
				print("length(riverIO) == 3")
			}
			
			if(length(riverIO) == 4){
				a = Watershed.IOR4(x=nodeIO, y=riverIO, dist=dist)
				print("length(riverIO) == 4")
			}	
			c1_inlet = a[[1]]
			c1_outlet= a[[2]]
		}else{
			riverIO = 0
			c1_inlet = 0
			c1_outlet = 0
		}
						
# identifying tributary watershedatersheds
if(class(c1_inlet)=="SpatialPointsDataFrame"){
	a = Watershed.Tributary(x=c1_inlet,xo= c1_outlet,y=riverIO,z=nodeIO,zhyd=zhyd_sb1, c1=c1)
	c2c3 = a[[1]]
	c2 = a[[2]]
	c3 = a[[3]]
	node_trib = a[[4]]
	}else{
		c2c3 = 0
		c2 = 0
		c3= 0
		node_trib = 0
	}
return(list(c1=c1,c1_inlet=c1_inlet,c1_outlet=c1_outlet,c2=c2,c3=c3,node_trib=node_trib,sb1=sb1,riverIO=riverIO,nodeIO=nodeIO,c1_river=c1_river,c1_node=c1_node))	
}
)
