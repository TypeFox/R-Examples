EntropyHOMO<-function(nameTF){
    require("MEET")
    homoentropy<-new.env()
    data("HomoEntropy",envir=homoentropy)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=homoentropy)[1],envir=homoentropy))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")
        model<-NULL
    
    }else{

        r<-slot(model,"model")<-get(ls(envir=homoentropy)[1],envir=homoentropy)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor
        slot(model,"model")<-r$model

}


    return(model)
    
}