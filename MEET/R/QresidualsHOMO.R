QresidualsHOMO<-function(nameTF){
    require("MEET")
    homoQresiduals<-new.env()
    data("HomoQresiduals",envir=homoQresiduals)
	model<-new("Model")
    if (is.na(summary(names(get(ls(envir=homoQresiduals)[1],envir=homoQresiduals))==nameTF)["TRUE"])=="TRUE"){

        print("Transcription Factor not included")

    }else{

        slot(model,"model")<-get(ls(envir=homoQresiduals)[1],envir=homoQresiduals)[[nameTF]]


    }
return(model)
    
}