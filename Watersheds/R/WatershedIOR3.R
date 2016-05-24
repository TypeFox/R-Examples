
# if(length(riverIO) == 3)
# x= nodeIO; y = riverIO
#	nodeIO
#	riverIO
#	length(riverIO)
# check code for eliminating closer and closerelev, maybe them are not necessary
Watershed.IOR3 = function(x, y, dist){
	nodeIO = x; nodeIO
	riverIO = y; riverIO
	dist = dist
	
	touch = SpDF_Touch(nodeIO, riverIO)
	touch = touch[[1]]
	id = which(touch[,2] == max(touch[,2])); id
	touch1 = vector()
	touch1 = rbind(touch1,touch[id,]);  length(touch1) # node with maximum rivers touching
	touchId = touch1[,1]; touchId # inlet node
	
	# identifying touching nodeIO and closer to boundary
		id = (dist == min(dist)); id
		id1 = (dist == sort(dist, F)[2]); id1
		id2 = id | id1; id2
		closerId = touch[id2,1]; closerId #probable inlet and outlet nodes
		id = (slot(nodeIO, "data")["OBJECTID"] == closerId[1]); id
		id1 = (slot(nodeIO, "data")["OBJECTID"] == closerId[2]); id1
		id2 = list(id | id1); id2
		c1_outinlet = SpDF_Subset(id2, nodeIO); c1_outinlet #probable inlet and outlet nodes

	# identyfing the distance to the second closer nodeIO
		c1_outinlet
		mini = vector() #closer node to inlet and outlet nodes
		for(i in 1:length(c1_outinlet)){
			objectId = slot(c1_outinlet, "data")["OBJECTID"][[1]][i]; objectId
			id2 = list(slot(c1_outinlet, "data")["OBJECTID"] == objectId); id2
			closeri = SpDF_Subset(id2, c1_outinlet); closeri
			neari = gDistance(closeri, nodeIO, byid=T); neari
			sorti = sort(neari); sorti
			id2 = list(neari == sorti[2]); id2
			closeri = SpDF_Subset(id2, nodeIO); closeri
			mini = rbind(mini, c(objectId, sorti[2], slot(closeri,"data")["OBJECTID"][[1]][1]))
		}

	# defining inlet or outlet node
		objectId = touchId; objectId
		object2Id = mini[which(mini[,1]!=objectId),1]; object2Id
		c1_inletId = objectId; c1_inletId
		c1_outletId = object2Id; c1_outletId								
		id = list(slot(nodeIO, "data")["OBJECTID"] == c1_inletId)
		c1_inlet = SpDF_Subset(id, nodeIO); c1_inlet
		id = list(slot(nodeIO, "data")["OBJECTID"] == c1_outletId)
		c1_outlet = SpDF_Subset(id, nodeIO); c1_outlet
		print("end WatershedIO_R3")
		return(list(inlet=c1_inlet, outlet=c1_outlet))					
}# end function