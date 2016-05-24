EntropyDROSOPHILA<-function(nameTF){
    require("MEET")
    drosophilaentropy<-new.env()
    data("DrosophilaEntropy",envir=drosophilaentropy)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=drosophilaentropy)[1],envir=drosophilaentropy))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")
        model<-NULL
    }else{

        r<-slot(model,"model")<-get(ls(envir=drosophilaentropy)[1],envir=drosophilaentropy)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor

        slot(model,"model")<-r$model

    }

    return(model)
    
}