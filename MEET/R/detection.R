        detection <-function (iicc){
	require("MEET")
	
	if(length(iicc$model)==0){
        
        if(is.null(iicc$nameTF)==TRUE) { 
        
                iicc$parameterIdeal<-iicc$parameters
	 
                iicc$model<-Models(iicc)
	 
                results<-Prediction(iicc)
        }else{
            
            p<-classMODEL(org=iicc$organism,method=iicc$method,nameTF=iicc$nameTF)
           
            if (is.null(p)=="TRUE"){
                        results<-NULL
                        print("Not model")
           
            }else{
    
                        iicc$model<-slot(p,"model")
        
                        if (iicc$method=="Entropy" | iicc$method=="Divergence"){
                    
                            iicc$Transcriptionfactor<-iicc$model$Transcriptionfactor
            }
            
                       
            results<-Prediction(iicc)
            }
        }
	}else{
      results<-Prediction(iicc)
	}
	
	return(results)
}

