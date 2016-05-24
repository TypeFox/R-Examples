# Date: 25.07.2013

# input x=id, y=zhyd. Identiying the index of true intersection
	SpDF_Subset = function(x, y){
									#x = list(gIntersects(zhyd, sb1, byid=T))
									#y = zhyd
									id1 = lapply(x, function(x) {
										id1 = which(x == "TRUE")
										})
										length(id1)[[1]]
										id1
									# subsetting zhyd by id1. zhyd that contains riverStation
 									 	n	 = length(id1[[1]])#; n
										nn = dim(y)[1]#; nn
										m  = dim(x[[1]])[2]#; m
										id2 = vector()

										for (j in 1:m-1){
											for (i in 1:n){
													if (id1[[1]][i] <= nn){id2 = cbind(id2, id1[[1]][i])} else if
														(id1[[1]][i] > j*nn & id1[[1]][i] <= (j+1)*nn ){id2 = cbind(id2, id1[[1]][i] - j*nn)}
												}
										}
										id2 = unique(as.vector(id2))
										y_id = y[id2,]
									}