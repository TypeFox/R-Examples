DivergenceRATTUS<-function(nameTF){
    require("MEET")
    rattusdivergence<-new.env()
    data("RattusDivergence",envir=rattusdivergence)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=rattusdivergence)[1],envir=rattusdivergence))==nameTF)["TRUE"])=="TRUE"
){

                    print("Transcription Factor not included")
                    model<-NULL
            }else{
                        
                    r<-slot(model,"model")<-get(ls(envir=rattusdivergence)[1],envir=rattusdivergence)[[nameTF]]
                    r$model$Transcriptionfactor<-r$Transcriptionfactor
                    slot(model,"model")<-r$model
    }


    return(model)
    
}
