"addPrelCl" <- function (thetaClass, model, th, po, addM=TRUE) 
{
     	## prelspec has structure 
	## list(list(what1, ind1, 
	##           what2, ind2, 
	##           rel, start), ...) 
     prelspec <- model@prelspec
     prel <- thetaClass@prel
     parvec <- getPar(model, po, th, thetaClass, addM) 
     if(!addM) 
	       parvec[model@mvecind[["prel"]]] <- prel[model@mvecind[["prel"]]]
     cnt <- 1
     thetaClass@prel <- parvec
     for(diffs in prelspec){
           if(length(diffs$rel) == 0 || diffs$rel == "lin"){
              constant <- 0
              if (length(diffs$start) > 1) {
		 constant<- thetaClass@prel[cnt+1]
	      }
	      newpar <- multiLin(thetaClass, diffs, thetaClass@prel[cnt]) + constant
	      if(length(diffs$ind1)==1)
		    slot(thetaClass, diffs$what1)[diffs$ind1] <- newpar 
	      if(length(diffs$ind1)==2) {
		    slot(thetaClass, diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- newpar 

                  }
	      cnt <- cnt + length(diffs$start)       
            }
	   else { 
             if(diffs$rel == "multilin"){
               newpar <- thetaClass@prel[cnt] + multiLin(thetaClass, diffs,thetaClass@prel[(cnt+1):(cnt+length(diffs$start)-1)] )
               if(length(diffs$ind1)==1)
                 slot(thetaClass, diffs$what1)[diffs$ind1] <- newpar 
               if(length(diffs$ind1)==2) 
                 slot(thetaClass, diffs$what1)[[diffs$ind1[1]]][diffs$ind1[2]] <- newpar 
	       cnt <- cnt + length(diffs$start) 
               
             }
           }
         }
     thetaClass
}

