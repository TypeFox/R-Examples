`PPT.Init` <-function(visible = TRUE,method=c("rcom","RDCOMClient"),addPres=TRUE){

    ppt <- list()
    ppt$method=match.arg(method)

    if(ppt$method=="rcom"){
    
    	if (!require(rcom)) {
       	warning("The package rcom is unavailable.")
	
		if(require(RDCOMClient)){
		warning("Using RDOMClient package instead of rcom.")
		}else{
		stop("Neither rcom or RDCOMClient packages are unavailable.")
		} 

    	}else{


	if("package:RDCOMClient" %in% search()){

		 warning("\nUsing rcom package. Detaching RDCOMClient package to avoid conflicts.")		
		 try(detach("package:RDCOMClient"))
	}


	}

    }


    if(ppt$method=="RDCOMClient"){
    
    	if (!require(RDCOMClient)) {
        	stop("The package RDCOMClient is unavailable. \n 
		To install RDCOMClient use:\n 
		install.packages('RDCOMClient' repos = 'http://www.omegahat.org/R')")
    	}

	if("package:rcom" %in% search()){

		 warning("Using RDCOMClient package. Detaching rcom package to avoid conflicts.")		
		 try(detach("package:rcom"))
	}

    }
	

    

    if(ppt$method=="rcom"){


	    ppt$ppt <- comCreateObject("PowerPoint.Application")
		
    		if(!comIsValidHandle(ppt$ppt)){
    
    			cat("Error: Unable to create valid handle for powerpoint presentation.\n")
    			cat("\n")
    			cat("Please ensure that you have installed the rcom package correctly.\n")
    			cat("\n")
    			cat("Run the following code in R:\n")
    
    			cat(">library(rcom)\n")
    			cat(">installstatconnDCOM()\n")
    			cat(">comRegisterRegistry()\n")
    			cat("\n")
    			cat("See 'http://learnserver.csd.univie.ac.at/rcomwiki/doku.php?id=wiki:how_to_install' for more info.\n")
    			cat("\n")

    			stop("Invalid handle for powerpoint presentation")
    		}
    

    }else{
    
	    ppt$ppt <- COMCreate("PowerPoint.Application")
	    ## Could do with an RDCOMClient validity check.

    }
    
    if(addPres){
    ppt$pres<-ppt$ppt[["Presentations"]]$add()
    }

    if (visible) {
    	ppt$ppt[["Visible"]]<-TRUE
    }

  
    return(invisible(ppt))
}

