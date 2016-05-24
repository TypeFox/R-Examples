QresidualsDROSOPHILA<-function(nameTF){
    require("MEET")
    drosophilaQresiduals<-new.env()
    data("DrosophilaQresiduals",envir=drosophilaQresiduals)
	model<-new("Model")
    if (is.na(summary(names(get(ls(envir=drosophilaQresiduals)[1],envir=drosophilaQresiduals))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")

    }else{

        slot(model,"model")<-get(ls(envir=drosophilaQresiduals)[1],envir=drosophilaQresiduals)[[nameTF]]


    }


    return(model)
    
}