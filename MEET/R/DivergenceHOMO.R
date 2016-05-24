DivergenceHOMO<-function(nameTF){
    require("MEET")
    homodivergence<-new.env()
    data("HomoDivergence",envir=homodivergence)
	model<-new("Model")
    if (is.na(summary(names(get(ls(envir=homodivergence)[1],envir=homodivergence))==nameTF)["TRUE"])=="TRUE"){
        print("Transcription Factor not included")
        model<-NULL
    }else{

        r<-slot(model,"model")<-get(ls(envir=homodivergence)[1],envir=homodivergence)[[nameTF]]
        r$model$Transcriptionfactor<-r$Transcriptionfactor    
        slot(model,"model")<-r$model
    }


    return(model)
    
}