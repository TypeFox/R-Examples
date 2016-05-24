classMODEL<-function(org,method,nameTF){
require("MEET")
model<-switch(org,
	 "Rattus norvegicus"=	switch(method,          	
								
								"Entropy"=EntropyRATTUS(nameTF),
								
         						"Divergence"=DivergenceRATTUS(nameTF),
                                
                                "Qresiduals"=QresidualsRATTUS(nameTF),
    
          						stop("Method not included")
        							 ),

	
	 "Mus musculus"= 	switch(method,   
    
                            "Entropy"=EntropyMUS(nameTF),
    
                            "Divergence"=DivergenceMUS(nameTF),
    
                            "Qresiduals"=QresidualsMUS(nameTF),
    
                            stop("Method not included")
    ),
    
          						
          						
          						
	 "Drosophila melanogaster"= switch(method,          	
								
                            "Entropy"=EntropyDROSOPHILA(nameTF),
    
                            "Divergence"=DivergenceDROSOPHILA(nameTF),
    
                            "Qresiduals"=QresidualsDROSOPHILA(nameTF),
    
                            stop("Method not included")
    ),
    
	 "Homo sapiens"= 	switch(method,          	
		
                            "Entropy"=EntropyHOMO(nameTF),
    
                            "Divergence"=DivergenceHOMO(nameTF),
    
                            "Qresiduals"=QresidualsHOMO(nameTF),
    
                            stop("Method not included")),
    
     stop("Organism not included"))
           
return(model)
}

