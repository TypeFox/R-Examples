
# Date: 08.08.2013
# Touch function. Identifies which nodes has touching lines
# and retrives a list with two elements. The first one a matrix 
# with the OBJECTID of the node (column 1), the maximum number of 
# lines that are touching the node (col2), and the elevation of 
# that node (col 3). The second one a matrix with the OBJECTID 
# of the lines that are touching the node

SpDF_Touch = function(x, y){
	# x = nodeIO*
	# x = c1_inlet
	# y = riverIO
	# y = tributary*
	# i = 1


	n = as.matrix(cbind(slot(x, "data")["OBJECTID"], slot(x, "data")["ELEV"])); n
	touch = vector()	# necessary review the length asignment
	riverTouch = vector()	# necessary review the length asignment
	#touch = matrix(0,1,3); dim(touch)
	#riverTouch = matrix(0,1,3); dim(riverTouch)	
	for(i in 1:(length(n)/2)){
		id = list(slot(x, "data")["OBJECTID"] == slot(x, "data")["OBJECTID"][[1]][i]); id
		node1 = SpDF_Subset(id, x); dim(node1)
		river = RiverStation(x=node1, y=y); dim(river)
		# if(i==1){
			# j = length(c(slot(x, "data")["OBJECTID"][[1]][i], c(dim(river)[1]), 
		# slot(x, "data")["ELEV"][[1]][i])); j
			# k = length(slot(river, "data")["OBJECTID"][[1]]); k
			# touch = matrix(0,1,j); dim(touch)
			# riverTouch = matrix(0,1,k); dim(riverTouch)
		# }
		if(class(river)!="SpatialLinesDataFrame"){#skipping no intersecting nodes
			next
		}	
		touch = rbind(touch, c(slot(x, "data")["OBJECTID"][[1]][i], c(dim(river)[1]), 
		slot(x, "data")["ELEV"][[1]][i]))
		riverTouch = rbind(riverTouch, slot(river, "data")["OBJECTID"][[1]])
	}
	#touch = touch[2:dim(riverTouch)[1],]; touch
	#riverTouch = riverTouch[2:dim(riverTouch)[1],]; riverTouch
	return(list(touch, riverTouch))
}
