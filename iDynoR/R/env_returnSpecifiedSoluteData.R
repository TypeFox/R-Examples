env_returnSpecifiedSoluteData <-
function(xmlResultFile,soluteRequested)
{
	# The solutes begin in the XML file at position xmlResultData[[1]][[2]]
	# Check this is valid input
	if(soluteRequested<=xmlSize(xmlResultFile[[1]])-1)
	{
		# Get the result tags for the specified solute
		# +1 as the first tag is not a solute tag
		soluteData<-xmlResultFile[[1]][[soluteRequested+1]][[1]]

		# Important here to grab the sizes of the grid, as this is how the data will be processed
		# This can be found in the header of each solute
		res<-env_returnSoluteGridRes(xmlResultFile,soluteRequested)
		nI<-env_returnSoluteGridIVoxels(xmlResultFile,soluteRequested)
		nJ<-env_returnSoluteGridJVoxels(xmlResultFile,soluteRequested)
		nK<-env_returnSoluteGridKVoxels(xmlResultFile,soluteRequested)

		# Now split the data by the semi-colon that separates each box
		result<-toString(soluteData)
		result<- strsplit(result,";\n")
		
		# Now remove the semicolon from the last entry
		result[[1]][length(result[[1]])]<- substr(result[[1]][length(result[[1]])], 1, nchar(result[[1]][length(result[[1]])])-1)

		# Counter to the position in the list of solute grid.
		gridPosCounter<-1

		soluteTable <- matrix(nrow=(nI*nJ), ncol=5)
		k<-0
		for(i in 0:(nI-1))
		{
			for(j in 0:(nJ-1))
			{
				soluteTable[gridPosCounter,1]<-soluteRequested
				soluteTable[gridPosCounter,2]<-i
				soluteTable[gridPosCounter,3]<-j
				soluteTable[gridPosCounter,4]<-k
				soluteTable[gridPosCounter,5]<-as.numeric(result[[1]][gridPosCounter])
			
				gridPosCounter<-gridPosCounter+1
			}
		}

		colnames(soluteTable)<-c("solute_id","location_x","location_y","location_z","concentration")

		return(soluteTable)
	}
	else
	{
		print("Error Reading Iterate")
	}
	
}
