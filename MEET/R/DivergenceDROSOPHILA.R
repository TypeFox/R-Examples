DivergenceDROSOPHILA<-function(nameTF){
    require("MEET")
    drosophiladivergence<-new.env()
    data("DrosophilaDivergence",envir=drosophiladivergence)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=drosophiladivergence)[1],envir=drosophiladivergence))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")
        model<-NULL
    }else{
    
        r<-slot(model,"model")<-get(ls(envir=drosophiladivergence)[1],envir=drosophiladivergence)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor
        slot(model,"model")<-r$model

}

    return(model)
    
}