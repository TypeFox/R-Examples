### R code from vignette source 'Watersheds.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: c1
###################################################
	library(Watersheds)
	data(WatershedsData)
	#ls()
	#str(WatershedsData)
	#str(WatershedsData["basin"])


	# plotting river Weser basin
	plot(WatershedsData["ctry"][[1]], col="gray60")
	plot(WatershedsData["basin"][[1]], col="gray30", add=TRUE)
	title("River Weser basin, Germany")

	# plotting subbasins river Weser basin
	plot(WatershedsData["basin"][[1]])
	plot(WatershedsData["subbasin"][[1]], col="gray60",add=TRUE)
	plot(WatershedsData["rWeser"][[1]],col="blue",lwd=2,add=TRUE)
	plot(WatershedsData["rAller"][[1]],col="blue",lwd=1,add=TRUE)
	plot(WatershedsData["rDiemel"][[1]],col="blue",lwd=1,add=TRUE)
	plot(WatershedsData["rFulda"][[1]],col="blue",lwd=1,add=TRUE)
	plot(WatershedsData["rHunte"][[1]],col="blue",lwd=1,add=TRUE)
	plot(WatershedsData["rWerra"][[1]],col="blue",lwd=1,add=TRUE)
	plot(WatershedsData["rWiumme"][[1]],col="blue",lwd=1,add=TRUE)
	title("Subbasins River Weser")

	# plotting primary zhyd watersheds and drainage network inside river Werra subbasin
	# subsetting the river Werra subbasin
	id = list(gIntersects(WatershedsData["rWerra"][[1]], WatershedsData["subbasin"][[1]],byid=TRUE))
	subbasin_rWerra = SpDF_Subset(id,WatershedsData["subbasin"][[1]])
	
	# subsetting the river Werra zhyd watersheds
	id = list(gIntersects(WatershedsData["rWerra"][[1]], WatershedsData["zhyd"][[1]],byid=TRUE))
	zhyd_rWerra = SpDF_Subset(id,WatershedsData["zhyd"][[1]])
	plot(subbasin_rWerra, col="grey60")
	plot(zhyd_rWerra,col="grey50",add=TRUE)
	plot(WatershedsData["rWerra"][[1]],col="blue",lwd=1,add=TRUE)
	title("Subbasin River Weser and primary zhyd watersheds")

	# subsetting the river Werra river drainage watersheds
	id = list(gIntersects(subbasin_rWerra, WatershedsData["river"][[1]],byid=TRUE))
	river_rWerra = SpDF_Subset(id,WatershedsData["river"][[1]])
	plot(subbasin_rWerra,col="grey60")
	plot(WatershedsData["rWerra"][[1]],col="blue",lwd=3,add=TRUE)
	plot(river_rWerra,col="blue1",add=TRUE)
	title("Subbasin River Weser and drainage network")


###################################################
### code chunk number 2: c2
###################################################
	station1 = WatershedsData["station"][[1]]
	subbasin1 = WatershedsData["subbasin"][[1]]
	zhyd1 = WatershedsData["zhyd"][[1]]
	river1 = WatershedsData["river"][[1]]
	node1 = WatershedsData["node"][[1]]
	
	station1 = SpatialPoints(station1, 
		proj4string=slot(subbasin1,"proj4string"))
	watershed = new("Watershed",station=station1,subbasin=subbasin1,
		zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)
	class(watershed)


###################################################
### code chunk number 3: c3
###################################################
	station1 = WatershedsData["station"][[1]]
	subbasin1 = WatershedsData["subbasin"][[1]]
	zhyd1 = WatershedsData["zhyd"][[1]]
	river1 = WatershedsData["river"][[1]]
	node1 = WatershedsData["node"][[1]]	

	station1 = SpatialPoints(coords=cbind(4328448.74, 3118576.86), 
		proj4string=slot(subbasin1,"proj4string"))
	watershed = new("Watershed",station=station1,subbasin=subbasin1,
		zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)

	a = Watershed.Order(watershed)
	c1 = a[[1]]
	c1_inlet = a[[2]]
	c1_outlet = a[[3]]
	c2 = a[[4]]
	c3 = a[[5]]
	node_trib = a[[6]]
	sb1 = a[[7]]
	riverIO = a[[8]]
	nodeIO = a[[9]]			
	c1_river = a[[10]]
	c1_node = a[[11]]	

	bbox1 = slot(c1, "bbox")
	bbox = matrix(0,2,2)
	bbox[,1] = bbox1[,1]*.998
	bbox[,2] = bbox1[,2]*1.002
	
	plot(c1, xlim=bbox[1,], ylim=bbox[2,],col="gray50")			
	plot(c2, col="gray75", add=TRUE)
	plot(c3, col="gray85", add=TRUE)
	plot(slot(watershed,"station"),pch=24, bg="blue",add= TRUE)
	plot.PolyLineAttribute(c1, "order", 450, 0.8)
	plot.PolyLineAttribute(c2, "order", 450, 0.8)
	plot.PolyLineAttribute(c3, "order", 450, 0.8)				
	plot(c1_river, col="blue", add=TRUE)
	plot(c1_node,pch=21,bg="blue",cex=.5,add=TRUE)
	plot(nodeIO,pch=21,bg="blue",cex=.5,add=TRUE)
	plot(c1_inlet, pch=21, bg="green",add= TRUE)
	plot(c1_outlet,pch=21, bg="red",add= TRUE)
	plot.PointAttribute(nodeIO,"ELEV",600,0.7)
	title(main="Current zhyd watershed (1)",
		sub="First order tributary watersheds (1.1, 1.2)")


###################################################
### code chunk number 4: c4
###################################################
	station1 = SpatialPoints(coords=cbind(4328650,3174450), 
		proj4string=slot(subbasin1,"proj4string"))
	watershed = new("Watershed",station=station1,subbasin=subbasin1,
		zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)

	a = Watershed.Order(watershed)
	c1 = a[[1]]
	node_trib = a[[6]]
	c1_river = a[[10]]

	watershed2 = new("Watershed", station=node_trib, subbasin=subbasin1, zhyd=zhyd1, river=river1, c1=c1,node=node1)
	c23 = Watershed.Order2(watershed2)
	c2 = c23[[1]]
	c3 = c23[[2]]	
	
	c2.0 = c2[[1]]
	c2_inlet = c2[[2]]
	c2_outlet = c2[[3]]
	c2.1 = c2[[4]]
	c2.2 = c2[[5]]
	c2_node_trib = c2[[6]]
	c2_sb1 = c2[[7]]
	c2_riverIO = c2[[8]]
	c2_nodeIO = c2[[9]]			
	c2_river = c2[[10]]
	c2_node = c2[[11]]	
			
	c3.0 = c3[[1]]
	c3_inlet = c3[[2]]
	c3_outlet = c3[[3]]
	c3.1 = c3[[4]]
	c3.2 = c3[[5]]
	c3_node_trib = c3[[6]]
	c3_sb1 = c3[[7]]
	c3_riverIO = c3[[8]]
	c3_nodeIO = c3[[9]]			
	c3_river = c3[[10]]
	c3_node = c3[[11]]	
	
	# subsetting river networks
	id = list(gIntersects(c2.1, WatershedsData$river,byid=TRUE))
	c21_river = SpDF_Subset(id,WatershedsData$river)

	id = list(gIntersects(c2.2, WatershedsData$river,byid=TRUE))
	c22_river = SpDF_Subset(id,WatershedsData$river)
	
	id = list(gIntersects(c3.1, WatershedsData$river,byid=TRUE))
	c31_river = SpDF_Subset(id,WatershedsData$river)

	id = list(gIntersects(c3.2, WatershedsData$river,byid=TRUE))
	c32_river = SpDF_Subset(id,WatershedsData$river)
	
	# plots
	bbox1 = slot(c3.2, "bbox")
	bbox = matrix(0,2,2)
	bbox[,1] = bbox1[,1]*.995
	bbox[,2] = bbox1[,2]*1.005
							
	plot(c1, col="gray50", xlim=bbox[1,], ylim=bbox[2,])	
	plot(c2.0, col = "gray95", add=TRUE)
	plot(c3.0, col="gray79", add=TRUE)
	plot(c2.1, col="gray78", add=TRUE)
	plot(c2.2, col="gray85", add=TRUE)
	plot(c3.1, col="gray53", add=TRUE)
	plot(c3.2, col="gray63", add=TRUE)

	plot(c1_river, col="blue",add=TRUE)		
	plot(c2_river, col="blue",add=TRUE)		
	plot(c3_river, col="blue",add=TRUE)			
	plot(c21_river, col="blue",add=TRUE)			
	plot(c22_river, col="blue",add=TRUE)			
	plot(c31_river, col="blue",add=TRUE)			
	plot(c32_river, col="blue",add=TRUE)	
	
	title(main="Current zhyd watershed and \n 1st and 2nd order tributary watersheds")


###################################################
### code chunk number 5: c5
###################################################
station1 = SpatialPoints(coords=cbind(4232972,3327634), 
		proj4string=slot(subbasin1,"proj4string"))
	watershed = new("Watershed",station=station1,subbasin=subbasin1,
		zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)

	a = Watershed.Order(watershed)
	c1 = a[[1]]
	nodeIO = a[[9]]			
	c1_river = a[[10]]
			
	# determining inlet and outlet watershed nodes
		# determining distances of nodeIO to c1
		boundary = gBoundary(c1)
		dist = gDistance(nodeIO, boundary, byid =TRUE)
		a = Watershed.IOR1(x=nodeIO, dist=dist)
		c1_inlet = a["inlet"][[1]]; c1_inlet
		c1_outlet = a["outlet"][[1]]; c1_outlet
		
	plot(c1,col="gray50")			
	plot(station1,pch=24, bg="blue",add= TRUE)
	plot(c1_river, col="blue", add=TRUE)
	plot(c1_outlet,pch=21, bg="red",add= TRUE)
	plot.PointAttribute(c1_outlet,"ELEV",700,0.8)
	title(main="Watershed outlet, case I")


###################################################
### code chunk number 6: c6
###################################################
	station1 = SpatialPoints(coords=cbind(4330341.36,3284797.06), 
			proj4string=slot(subbasin1,"proj4string"))
		watershed = new("Watershed",station=station1,subbasin=subbasin1,
			zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)
	
		a = Watershed.Order(watershed)
		c1 = a[[1]]
		nodeIO = a[[9]]			
		c1_river = a[[10]]
		c1_node = a[[11]]	
				
		# determining inlet and outlet watershed nodes
			# determining distances of nodeIO to c1
			boundary = gBoundary(c1)
			dist = gDistance(nodeIO, boundary, byid =TRUE)
			a = Watershed.IOR2(x=nodeIO, dist=dist, node=c1_node)	
		c1_inlet = a["inlet"][[1]]; c1_inlet
		c1_outlet = a["outlet"][[1]]; c1_outlet
			
		plot(c1,col="gray60")			
		plot(station1,pch=24, bg="blue",add= TRUE)
		plot(c1_river, col="blue", add=TRUE)
		plot(c1_outlet,pch=21, bg="red",add= TRUE)
		plot.PointAttribute(c1_outlet,"ELEV",700,0.8)
		title(main="Watershed outlet, case II")


###################################################
### code chunk number 7: c7
###################################################
	station1 = SpatialPoints(coords=cbind(4217199.42,3353511.83), 
		proj4string=slot(subbasin1,"proj4string"))
	watershed = new("Watershed",station=station1,subbasin=subbasin1,
		zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)

	a = Watershed.Order(watershed)
	c1 = a[[1]]
	riverIO = a[[8]]
	nodeIO = a[[9]]			
	c1_river = a[[10]]
			
	# determining inlet and outlet watershed nodes
		# determining distances of nodeIO to c1
		boundary = gBoundary(c1)
		dist = gDistance(nodeIO, boundary, byid =TRUE)
		a = Watershed.IOR3(x=nodeIO, y=riverIO, dist=dist)	
		c1_inlet = a["inlet"][[1]]; c1_inlet
		c1_outlet = a["outlet"][[1]]; c1_outlet
		
	plot(c1,col="gray60")			
	plot(station1,pch=24, bg="blue",add= TRUE)
	plot(c1_river, col="blue", add=TRUE)
	plot(c1_outlet,pch=21, bg="red",add= TRUE)
	plot(c1_inlet,pch=21, bg="green",add= TRUE)
	plot.PointAttribute(c1_outlet,"ELEV",1000,0.8)
	plot.PointAttribute(c1_inlet,"ELEV",1000,0.8)
	title(main="Watershed outlet and inlet, case III")	


###################################################
### code chunk number 8: c8
###################################################
		station1 = SpatialPoints(coords=cbind(4357947,3284525), 
			proj4string=slot(subbasin1,"proj4string"))
		watershed = new("Watershed",station=station1,subbasin=subbasin1,
			zhyd=zhyd1,river=river1,c1=subbasin1,node=node1)
	
		a = Watershed.Order(watershed)
		c1 = a[[1]]
		riverIO = a[[8]]
		nodeIO = a[[9]]			
		c1_river = a[[10]]
				
		# determining inlet and outlet watershed nodes
			# determining distances of nodeIO to c1
			boundary = gBoundary(c1)
			dist = gDistance(nodeIO, boundary, byid =TRUE)
			a = Watershed.IOR4(x=nodeIO, y=riverIO, dist=dist)
		c1_inlet = a["inlet"][[1]]; c1_inlet
		c1_outlet = a["outlet"][[1]]; c1_outlet
			
		plot(c1,col="gray60")			
		plot(station1,pch=24, bg="blue",add= TRUE)
		plot(c1_river, col="blue", add=TRUE)
		plot(c1_outlet,pch=21, bg="red",add= TRUE)
		plot(c1_inlet,pch=21, bg="green",add= TRUE)
		plot.PointAttribute(c1_outlet,"ELEV",1000,0.8)
		plot.PointAttribute(c1_inlet,"ELEV",1000,0.8)
		title(main="Watershed outlet and inlet, case IV")


