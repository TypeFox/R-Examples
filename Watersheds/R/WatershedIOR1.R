	# determining inlet and output nodes
		# determining distances of nodeIO to c1
Watershed.IOR1 = function(x, dist){
	nodeIO=x
	id = list(dist == min(dist)); id
	distorder = sort(dist); distorder
	if(distorder[1]==distorder[2]){
		dist
		outinlet = SpDF_Subset(id, nodeIO); outinlet
		id =list(slot(outinlet, "data")["ELEV"]==min(slot(outinlet, "data")["ELEV"]))
		c1_outlet = SpDF_Subset(id, nodeIO); c1_outlet
		c1_inlet = 0
	}else{
	c1_outlet = SpDF_Subset(id, nodeIO)
	c1_inlet = 0
	}
	return(list(inlet=c1_inlet, outlet=c1_outlet))
}