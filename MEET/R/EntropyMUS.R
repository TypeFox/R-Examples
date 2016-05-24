EntropyMUS<-function(nameTF){
    require("MEET")
    musentropy<-new.env()
    data("MusEntropy",envir=musentropy)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=musentropy)[1],envir=musentropy))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")
        model<-NULL

    }else{
        r<-slot(model,"model")<-get(ls(envir=musentropy)[1],envir=musentropy)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor
        slot(model,"model")<-r$model

}

    return(model)
    
}