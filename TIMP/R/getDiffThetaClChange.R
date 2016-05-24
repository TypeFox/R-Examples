"getDiffThetaClChange" <-
function (th, parorder, modellist, thetaClasslist, diffchange) 
{
      pcnt <- 1
      for(diffs in diffchange) {
	  d <- diffs$dataset[1]
	  parchange <- modellist[[d]]@parnames
	  
	  if(diffs$what %in% parchange){
		 parslot <- getPar(modellist[[d]],
		 parorder[[pcnt]], th, thetaClasslist[[d]])
		 pcnt <- pcnt + 1  
		 if(is.list( slot(modellist[[d]], diffs$what))) {
			tlist <- list()
			cnt <- 1
			
			for(i in 1:length(slot(modellist[[d]],diffs$what))){ 
			      	tlist[[i]] <- parslot[cnt:(cnt + 
					    length(slot(modellist[[d]], 
					    diffs$what)[[i]]) -1)]
				cnt <- cnt + length(slot(modellist[[d]], 
				diffs$what)[[i]])
			}
			parslot <- tlist
		 }
		 for (datai in diffs$dataset)
			  slot(thetaClasslist[[datai]], diffs$what) <- parslot
           }
	   if(diffs$what == "prelspec"){
 	  	   thL <- addPrelCl(thetaClasslist[[d]], modellist[[d]], 
		   th, parorder[[pcnt]])
		   for (datai in diffs$dataset)
			  thetaClasslist[[datai]] <- thL
		   pcnt <- pcnt + 1  
          }
       }
      
       thetaClasslist
}

