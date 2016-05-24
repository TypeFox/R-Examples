# determining inlet and output nodes
Watershed.IOR2 = function(x, dist, node){	
	nodeIO = x
	c1_node = node
	id = list(dist == min(dist)); id
	outinlet = SpDF_Subset(id, nodeIO)
	outinletElev = slot(outinlet, "data")["ELEV"]; outinletElev
	#plot.PointAttribute(outinlet, "ELEV")
	dist = gDistance(outinlet, c1_node, byid=T); dist
		#plot(c1)
		#plot(c1_node, add=T)
		#plot.PointAttributte(c1_node, "ELEV")
	id = list(dist == max(dist)); id
	farther = SpDF_Subset(id, c1_node); farther
	fartherElev = slot(farther, "data")["ELEV"]; fartherElev
		#averageElev = slot(c1_node, "data")["ELEV"]; averageElev
		#averageElev = mean(averageElev[[1]]); averageElev
	if(fartherElev > outinletElev & max(slot(nodeIO, "data")["ELEV"])>5){
		c1_inlet = 0
		c1_outlet = outinlet	
	}else if(max(slot(nodeIO, "data")["ELEV"])<5){# probable coastal basin
		c1_inlet = 0
		c1_outlet = 0
	}else{
		c1_outlet = 0	
		c1_inlet = outinlet	
	}
	return(list(inlet=c1_inlet, outlet=c1_outlet))
}		
	