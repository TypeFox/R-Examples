QresidualsMUS<-function(nameTF){
    require("MEET")
    musQresiduals<-new.env()
    data("MusQresiduals",envir=musQresiduals)
	model<-new("Model")
    
    if (is.na(summary(names(get(ls(envir=musQresiduals)[1],envir=musQresiduals))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")

    }else{

        slot(model,"model")<-get(ls(envir=musQresiduals)[1],envir=musQresiduals)[[nameTF]]

    }

    return(model)
    
}