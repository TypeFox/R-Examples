
# if(length(riverIO) == 4)
# x= nodeIO; y = riverIO
#	nodeIO
#	riverIO
#	length(riverIO)
Watershed.IOR4 = function(x, y, dist){
	nodeIO = x
	riverIO = y
	
	touch = SpDF_Touch(nodeIO, riverIO)
	touch = touch[[1]] ##review
	id = which(touch[,2] == max(touch[,2])); id
	maxtouch = max(touch[,2]); maxtouch
	touch1 = vector()
	touch1 = rbind(touch1,touch[id,]);  length(touch1) # node with maximum rivers touching
	touchId = touch1[,1]; touchId

	# identifying touching nodeIO and closer to boundary
		dist
		nodeIO
		id = (dist == min(dist)); id
		id1 = (dist == sort(dist, F)[2]); id1
		id2 = id | id1; id2
		closerId = touch[id2,1]; closerId #probable inlet and outlet nodes
		id = (slot(nodeIO, "data")["OBJECTID"] == closerId[1]); id
		id1 = (slot(nodeIO, "data")["OBJECTID"] == closerId[2]); id1
		id2 = list(id | id1); id2
		c1_outinlet = SpDF_Subset(id2, nodeIO); c1_outinlet #probable inlet and outlet nodes
	
	## identyfing the distance to the second closer nodeIO
		# c1_outinlet
		# mini = vector() #closer node to inlet and outlet nodes
		# for(i in 1:length(c1_outinlet)){
			# objectId = slot(c1_outinlet, "data")["OBJECTID"][[1]][i]; objectId
			# id2 = list(slot(c1_outinlet, "data")["OBJECTID"] == objectId); id2
			# closeri = SpDF_Subset(id2, c1_outinlet); closeri
			# disti = gDistance(closeri, nodeIO, byid=T); disti
			# sorti = sort(disti); sorti
			# id2 = list(disti == sorti[2]); id2
			# closeri = SpDF_Subset(id2, nodeIO); closeri
			# mini = rbind(mini, c(objectId, sorti[2], slot(closeri,"data")["OBJECTID"][[1]][1]))
		# }
		# mini
	# identifying inlet and outlet nodes if maxtouch==3
		if(maxtouch > 2){# just a inlet node could have 3 or more rivers touching
			objectId = touch[id,1]; objectId						
			id = list(slot(c1_outinlet, "data")["OBJECTID"]!=objectId); id						
			c1_outlet = SpDF_Subset(id, c1_outinlet); 
			c1_inletId = objectId; c1_inletId
			id = list(slot(nodeIO, "data")["OBJECTID"] == c1_inletId)
			c1_inlet = SpDF_Subset(id, nodeIO); c1_inlet
			return(list(inlet=c1_inlet, outlet=c1_outlet))					
		}
	# identyfing the distance to the farther nodeIO
		c1_outinlet
		maxi = vector() #farther nodeIO to inlet and outlet nodes
		for(i in 1:length(c1_outinlet)){
			objectId = slot(c1_outinlet, "data")["OBJECTID"][[1]][i]; objectId
			id2 = list(slot(c1_outinlet, "data")["OBJECTID"] == objectId); id2
			closeri = SpDF_Subset(id2, c1_outinlet); closeri
			disti = gDistance(closeri, nodeIO, byid=T); disti
			sorti = sort(disti, T); sorti
			id2 = list(disti == sorti[1]); id2
			closeri = SpDF_Subset(id2, nodeIO); closeri
			maxi = rbind(maxi, c(objectId, sorti[1], slot(closeri,"data")["OBJECTID"][[1]][1]))
		}
		maxi
	# identifying the elevation of the second closer nodeIO
	# for defining inlet or outlet node
		maxi
		objectId = maxi[1,1]; objectId						
		fartherId = maxi[1,3]; fartherId
		object2Id = maxi[2,1]; object2Id						
		farther2Id = maxi[2,3]; farther2Id
		objectIdElev = touch[which(touch[,1]==objectId),3]; objectIdElev
		fartherIdElev = touch[which(touch[,1]==fartherId),3]; fartherIdElev
		object2IdElev = touch[which(touch[,1]==object2Id),3]; object2IdElev	
		farther2IdElev = touch[which(touch[,1]==farther2Id),3]; farther2IdElev												
		if (fartherIdElev < farther2IdElev){
			c1_inletId = objectId; c1_inletId
			c1_outletId = object2Id; c1_outletId								
			id = list(slot(nodeIO, "data")["OBJECTID"] == c1_inletId)
			c1_inlet = SpDF_Subset(id, nodeIO)
			id = list(slot(nodeIO, "data")["OBJECTID"] == c1_outletId)
			c1_outlet = SpDF_Subset(id, nodeIO)
			return(list(inlet=c1_inlet, outlet=c1_outlet))					
		}else{
			c1_inletId = object2Id; c1_inletId
			c1_outletId = objectId; c1_outletId								
			id = list(slot(nodeIO, "data")["OBJECTID"] == c1_inletId)
			c1_inlet = SpDF_Subset(id, nodeIO)
			id = list(slot(nodeIO, "data")["OBJECTID"] == c1_outletId)
			c1_outlet = SpDF_Subset(id, nodeIO)
			return(list(inlet=c1_inlet, outlet=c1_outlet))					
		}
} # end function