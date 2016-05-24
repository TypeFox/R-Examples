EntropyRATTUS<-function(nameTF){
    require("MEET")
    rattusentropy<-new.env()
    data("RattusEntropy",envir=rattusentropy)
	model<-new("Model")

    if (is.na(summary(names(get(ls(envir=rattusentropy)[1],envir=rattusentropy))==nameTF)["TRUE"])=="TRUE"){
             
                print("Transcription Factor not included")
                model<-NULL
       
        }else{
                r<-slot(model,"model")<-get(ls(envir=rattusentropy)[1],envir=rattusentropy)[[nameTF]]
                r$model$Transcriptionfactor<-r$Transcriptionfactor
                slot(model,"model")<-r$model
    }

  return(model)

}
