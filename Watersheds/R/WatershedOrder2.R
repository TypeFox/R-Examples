
# Date: 09.08.2013

# Identifying third level tributaries Watershed

setGeneric("Watershed.Order2", function(watershed,...) standardGeneric("Watershed.Order2"))

setMethod("Watershed.Order2", signature=c("Watershed"),
	function(watershed,...){
		station = slot(watershed,"station")
		subbasin = slot(watershed,"subbasin")
		zhyd = slot(watershed,"zhyd")
		river = slot(watershed,"river")
		c1 = slot(watershed,"c1")
		node = slot(watershed,"node")
			
			#station=node_trib
			
			if(length(station)==2){
				tribId = slot(station,"data")["OBJECTID"]; tribId
				id = list(slot(station,"data")["OBJECTID"]==tribId[[1]][1])
				station2 = SpDF_Subset(id, station)
				station2 = SpatialPoints(coords=slot(station2,"coords"), proj4string=slot(station2,"proj4string"))
				id = list(slot(station,"data")["OBJECTID"]==tribId[[1]][2])
				station3 = SpDF_Subset(id, station)
				station3 = SpatialPoints(coords=slot(station3,"coords"), proj4string=slot(station3,"proj4string"))


				# identifying current watershed that contains "station" with multicore
					##begin parallel process
					##time1 = proc.time()
					#id =  list(gIntersects(zhyd, station2, byid=T))
					#c2 = parallel(SpDF_Subset(id, zhyd)) # current watershed
					#id2 =  list(gIntersects(zhyd, station3, byid=T))
					#c3 = parallel(SpDF_Subset(id2, zhyd)) # current watershed
					#c23 = collect(list(c2, c3))
					#time2 = proc.time()
					#time2 - time1									
				
				# identifying current watershed that contains "station", without multicore					
					id  = list(gIntersects(zhyd, station2, byid=T))
					c2  = (SpDF_Subset(id, zhyd)) # current watershed
					id2 = list(gIntersects(zhyd, station3, byid=T))
					c3  = (SpDF_Subset(id2, zhyd)) # current watershed
					c23 = (list(c2, c3))


				watershed1 = new("Watershed",station=station2,subbasin=subbasin,zhyd=zhyd,river=river,c1=c23[[1]],node=node)
				watershed2 = new("Watershed",station=station3,subbasin=subbasin,zhyd=zhyd,river=river,c1=c23[[2]],node=node)


			## begin parallel process with multicore
				##time1 = proc.time()
				##c2 = Watershed.Order(watershed1)		
				##c3 = Watershed.Order(watershed2)
				 #c2 = parallel(Watershed.Order(watershed1))		
				 #c3 = parallel(Watershed.Order(watershed2))
				 #c23 = collect(list(c2, c3))
				## str(c23)
				##time2 = proc.time()
				##time2 - time1									

			## without parallel process with multicore
				 c2  = (Watershed.Order(watershed1))		
				 c3  = (Watershed.Order(watershed2))
				 c23 = (list(c2, c3))
				}
		return(c23)
		}
)