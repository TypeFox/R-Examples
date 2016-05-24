QresidualsRATTUS<-function(nameTF){
    require("MEET")
    rattusQresiduals<-new.env()
    data("RattusQresiduals",envir=rattusQresiduals)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=rattusQresiduals)[1],envir=rattusQresiduals))==nameTF)["TRUE"])=="TRUE"){

                print("Transcription Factor not included")

            }else{

                slot(model,"model")<-get(ls(envir=rattusQresiduals)[1],envir=rattusQresiduals)[[nameTF]]

    }
   
    return(model)
    
}
	


