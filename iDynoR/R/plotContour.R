plotContour <-
function(resultFileFolder,timepoint,soluteReqd,folderForGraphOut)
{
	if(is.numeric(timepoint))
	{
		# Set up the graph file location
		GRAPHFILE = paste(folderForGraphOut,"/Contour_Solute_",soluteReqd,"_Iteration_",timepoint,".pdf",sep="")
		pdf(GRAPHFILE,width=10,height=7)

		xmlResultFile<-readSimResultFile(resultFileFolder,"env_State",timepoint)

		soluteData<-env_returnSpecifiedSoluteData(xmlResultFile,soluteReqd)
	
		nums<-c(1:nrow(soluteData))
		minix <- min(soluteData[,2])
		maxix <- max(soluteData[,2])
		x <- minix:maxix 
		miniy <- min(soluteData[,3])
		maxiy <- max(soluteData[,3])
		y <- miniy:maxiy
		temp <- matrix(nrow=length(y), ncol=length(x))
		for(i in x)
		{
			for(j in y)
			{
				if(minix==0)
				{
					iv <- i+1
			 	}
			  	else
				{
					iv <- i
			  	}
	  
				if(miniy==0)
				{
					jv <- j+1
	  			}
	  			else
				{
					jv <- j
			  	}
	  		
				temp[iv,jv] <- soluteData[,5][which(soluteData[,2]==i)][which(soluteData[,3][which(soluteData[,2]==i)]==j)]
			}
		}

		filled.contour(1:dim(temp)[2],1:dim(temp)[1],t(temp))

		dev.off()

	}
	else
	{
		print("Error: Iteration Not Recognised")
	}



}
