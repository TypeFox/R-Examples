### Classdefinition inherited_stfdf
inherited_stfdf <- setClass("inherited_stfdf", 		# space-time full data frame
	slots =c(ValueIDs = "STFDF", DerivedFromIDs = "STFDF", MetadataRel= "STFDF", Metadata = "data.frame"),		
  contains = "STFDF",									# superclass
  
  validity = function(object) {  
	sp = object@sp@coords
	spValueID = object@ValueIDs@sp@coords
	spDerID = object@DerivedFromIDs@sp@coords
	
	time = object@time
	timeValueID = object@ValueIDs@time
	timeDerID = object@DerivedFromIDs@time
	nrowtime <- nrow(object@time)
	
	data = nrow(object@data)
	dataValueID = nrow(object@ValueIDs@data)
	dataDerID = nrow(object@DerivedFromIDs@data)
  
  stopifnot(
	# content of SP and Time must be the same. The data must have the same length. Otherwise the inherited_stfdf is wrong generated.
	(sp == spValueID) &&
	(spValueID == spDerID) &&
	(time == timeValueID) &&
	(timeValueID == timeDerID) &&
	(data == dataValueID) &&
	( dataValueID == dataDerID) #&&
	#(data == nrowtime )
	)
    return(TRUE)
  }
)

### Constructor class inherited_stfdf
inherited_stfdf <- function(sp, time, data, endtime, ValueIDs, DerivedFromIDs, MetadataRel, Metadata)
	new("inherited_stfdf",
		sp = sp, time = time, data = data, endTime = endtime,
									ValueIDs=ValueIDs,
									DerivedFromIDs=	DerivedFromIDs,
									MetadataRel = MetadataRel,
									Metadata=Metadata
		)


# Sub-Setting
setMethod("[", "inherited_stfdf",
function(x, i, j, ... ){
	# callNextMethod call the method defined for the parent class
	 object <- callNextMethod(x,i,j,...)						# Overwrite-method. get Slots from the main-data-stfdf

# wenn weitere Parameter angegeben sind (ausser Zeit und Ort) wird auch nach allen Metadaten gesucht.
# Wenn hingegen nur eine bestimmte Variable selectiert ist, sollen natuerlich auch die Metadaten reduziert werden.
	
	dots = list(...)
	#check if variables are choosed
	if (length(dots) > 0) {
		missing.k = FALSE
		k = dots[[1]]
	} else{
		missing.k = TRUE
	}
	if(missing.k){
		meta.new <- x@Metadata
		}
	else {
		if(class(k) == 'character'){
			selAtChr = which(x@Metadata$variable == k)
			selection = x@Metadata[selAtChr,]
		}
		else{
		getColName = colnames(x@data[k])
			selAtChr = which(x@Metadata$variable == getColName)
			selection = x@Metadata[selAtChr,]
		}
		meta.new <- selection
	}

	if("xts" %in% class(object)){
		warning("If only one spatial point is selected a conversion to a spacetime-object is impossible. So a xts is returned.")
		to.return = object
	}
	else if(class(object) == "numeric"){
		warning("If only one spatial point and one time is selected a conversion to a spacetime-object is impossible. So numeric is returned.")
		to.return = object
	}
	else if(class(object) == "SpatialPointsDataFrame"){
		warning("If only one time point is selected a conversion to a spacetime-object is impossible. So a SpatialPointsDataFrame is returned.")
		to.return = object
	}
	else if(class(object) == "data.frame"){
		warning("If only one time point and one location is selected a conversion to a spacetime-object is impossible. So a data frame is returned.")
		to.return = object
	}
	else{
		stfdf = object
	 to.return <- inherited_stfdf(sp=stfdf@sp, time=stfdf@time, data=stfdf@data, endtime = stfdf@endTime,
	                    ValueIDs=x@ValueIDs[i,j,...],
						DerivedFromIDs=x@DerivedFromIDs[i,j,...],
						MetadataRel = x@MetadataRel[i,j,...],
						Metadata=meta.new)
						} 
	 return(to.return)

 }
 )
 
 
 # compare two datasets
 setMethod("==", signature(e1="inherited_stfdf",e2="inherited_stfdf"),
	function(e1,e2){		

		sp1 = e1@sp@coords
		sp2 = e2@sp@coords
		spValueID1 = e1@ValueIDs@sp@coords
		spValueID2 = e2@ValueIDs@sp@coords
		spDerID1 = e1@DerivedFromIDs@sp@coords
		spDerID2 = e2@DerivedFromIDs@sp@coords
		
		time1 = e1@time
		time2 = e2@time
		timeValueID1 = e1@ValueIDs@time
		timeValueID2 = e2@ValueIDs@time
		timeDerID1 = e1@DerivedFromIDs@time
		timeDerID2 = e2@DerivedFromIDs@time
		
		# test the content of sp and time. Thus current only updates of data works.

		if(!(											
			(sp1 == sp2) &&
			(spValueID1 == spValueID2) &&
			(spDerID1 == spDerID2) &&
			
			(time1 == time2) &&
			(timeValueID1 == timeValueID2) &&
			(timeDerID1 == timeDerID2) 
		    )
		) return(FALSE)
		
		#comparing values	<- to-do: do the same for e1@sp, e1@time, e1@DerivedFromIDs@data, e1@ValueIDs@data!
		data1 <- e1@data
		data1[is.na(data1)]<-"NA"
			
		data2 <- e2@data
		data2[is.na(data2)]<-"NA"
		
		# compare dataframe
		# compare both objects. different values get the boolean false. The same get true.
		if(any(dim(data1)!=dim(data2))) return(FALSE)

		different <- e1@data != e2@data
				
		return(different)
	}
)
