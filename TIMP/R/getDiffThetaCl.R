"getDiffThetaCl" <-
function(th, thetaClasslist, mod) {
	modeldiff <- mod@modeldiffs 
	modellist <- mod@modellist 
	parorder <- mod@parorderdiff 
	parorderchange <- mod@parorderchange
	pcnt <- 1  
	if(length(modeldiff$change)!=0)   
	   thetaClasslist <- getDiffThetaClChange(th, parorderchange,
	   modellist, thetaClasslist, modeldiff$change)
	if(length(modeldiff$free)!=0 || length(modeldiff$add) != 0)
		for(diff in append(modeldiff$free, modeldiff$add)){ 
			  
			  partmp <- th[parorder[[pcnt]]$ind]	
			  removepar <- parorder[[pcnt]]$rm 
			  pcnt <- pcnt + 1
			  if(diff$what %in% modellist[[diff$dataset[1]]]@positivepar)
			     partmp <- exp(partmp) 
			  for(fx in removepar){
				 if(fx %in% modellist[[diff$dataset]]@fvecind[[diff$what]])
				       partmp <- append(partmp,
				       unlist(slot(modellist[[diff$dataset]], 
				       diff$what))[fx], after=(fx-1))
				  else  ## add prel code, this won't work 
					partmp <- append(partmp, 0, after=(fx-1))
			  }
			  if(length(diff$ind) == 2)
			     for (d in diff$dataset) 
				slot(thetaClasslist[[d]], 
				diff$what)[[diff$ind[1]]][diff$ind[2]] <- partmp
			  if(length(diff$ind) == 1)
			     for (d in diff$dataset) 
				slot(thetaClasslist[[d]], 
				diff$what)[diff$ind]  <- partmp
          
		}
	if(length(modeldiff$remove)!=0)
		for(diff in modeldiff$remove){ 
			  if(length(diff$ind) == 2) 
			     for (d in diff$dataset) 
				slot(thetaClasslist[[d]], 
				diff$what)[[diff$ind[1]]][- diff$ind[2]] 
			  if(length(diff$ind) == 1) 
			     for (d in diff$dataset) 
				slot(thetaClasslist[[d]], 
				diff$what)[-diff$ind]  
		}
	if(length(modeldiff$rel)!=0)
	  for(diff in modeldiff$rel){
	   ds1 <- diff$dataset1
	   ds2 <- diff$dataset2
           if(length(diff$rel) == 0 || diff$rel == "lin"){
	      if( length(diff$fix) !=0)
		   thscal <- diff$start
	      else 
		   thscal <- th[parorder[[pcnt]]$ind]
	      pcnt <- pcnt+1
	      if(length(diff$ind1)==1 && length(diff$ind2)==1){
		for(i in 1:length(ds1)){
                    slot(thetaClasslist[[ds1[i]]], 
	            diff$what1)[diff$ind1] <- slot(thetaClasslist[[ds2[i]]],
	diff$what2)[diff$ind2] * thscal[1] + thscal[2]

			       }
	      } 
	      if(length(diff$ind1)==1 && length(diff$ind2)==2){
		   for(i in 1:length(ds1))
		    slot(thetaClasslist[[ds1[i]]], 
	            diff$what1)[diff$ind1] <- 
	            slot(thetaClasslist[[ds2[i]]], 
	            diff$what2)[[diff$ind2[1]]][diff$ind2[2]] * thscal[1] + thscal[2]	
	      }
	      if(length(diff$ind1)==2 && length(diff$ind2)==1){
		  for(i in 1:length(ds1))
                    slot(thetaClasslist[[ds1[i]]], 
	            diff$what1)[[diff$ind1[1]]][diff$ind1[2]] <- 
	            slot(thetaClasslist[[ds2[i]]], 
	            diff$what2)[diff$ind2] * thscal[1] + thscal[2] 
	      }
	      if(length(diff$ind1)==2 && length(diff$ind2)==2){
	            for(i in 1:length(ds1))
		    slot(thetaClasslist[[ds1[i]]], 
	            diff$what1)[[diff$ind1[1]]][diff$ind1[2]] <- 
	            slot(thetaClasslist[[ds2[i]]], 
	            diff$what2)[[diff$ind2[1]]][diff$ind2[2]] * 
		    thscal[1] + thscal[2] 
              }
           }       
	}
	if(length(modeldiff$dscal) != 0) {
		for(i in 1:length(modeldiff$dscal)) {
		       parvec <- getPar(modellist[[modeldiff$dscal[[i]]$to]], 
					parorder[[pcnt]], th, thetaClasslist[[i]])
		       pcnt <- pcnt + 1
		       slot(thetaClasslist[[modeldiff$dscal[[i]]$to]], 
			    "drel") <- parvec  
		}
        }
	thetaClasslist
}

