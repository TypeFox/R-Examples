DivergenceMUS<-function(nameTF){
    require("MEET")
    musdivergence<-new.env()
    data("MusDivergence",envir=musdivergence)
	model<-new("Model")
   
    if (is.na(summary(names(get(ls(envir=musdivergence)[1],envir=musdivergence))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")
        model<-NULL

    }else{

        r<-slot(model,"model")<-get(ls(envir=musdivergence)[1],envir=musdivergence)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor
        slot(model,"model")<-r$model
}


    
    return(model)
    
}