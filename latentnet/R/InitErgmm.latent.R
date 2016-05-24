InitErgmm.euclidean<-function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  if(nargs()<2)
    stop(paste("euclidean() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  if(d<=0) stop("Invalid latent space dimensionality given", call.=FALSE)
  if(!is.null(model[["d"]]) && model[["d"]]>0){
    stop("Only one latent position term can be added!")
  }

  model[["latent"]] <- "negative.Euclidean"
  
  model[["d"]] <- d
  model[["G"]] <- G

  model[["prior"]][["Z.var.mul"]]<-var.mul
  model[["prior"]][["Z.var"]]<-var
  model[["prior"]][["Z.var.df.mul"]]<-var.df.mul
  model[["prior"]][["Z.var.df"]]<-var.df
  model[["prior"]][["Z.mean.var.mul"]]<-mean.var.mul
  model[["prior"]][["Z.mean.var"]]<-mean.var
  model[["prior"]][["Z.pK.mul"]]<-pK.mul
  model[["prior"]][["Z.pK"]]<-pK
 
  if(!("Z.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.var"]]<-model[["prior"]][["Z.var.mul"]]*(network.size(model[["Yg"]])/max(1,model[["G"]]))^(2/model[["d"]])
  if(!("Z.mean.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.mean.var"]]<-model[["prior"]][["Z.mean.var.mul"]]*model[["prior"]][["Z.var"]]*max(1,model[["G"]])^(2/model[["d"]])
  if(!("Z.var.df" %in% names(model[["prior"]]))) model[["prior"]][["Z.var.df"]]<-model[["prior"]][["Z.var.df.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  if(!("Z.pK" %in% names(model[["prior"]]))) model[["prior"]][["Z.pK"]]<-model[["prior"]][["Z.pK.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  
  
  model
}


InitErgmm.bilinear<-function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  if(nargs()<2)
    stop(paste("bilinear() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  if(d<=0) stop("Invalid latent space dimensionality given", call.=FALSE)
  if(!is.null(model[["d"]]) && model[["d"]]>0){
    stop("Only one latent position term can be added!")
  }

  model[["latent"]] <- "bilinear"
  
  model[["d"]] <- d
  model[["G"]] <- G

  if(G>0) warning("For a bilinear (inner-product) latent position effect, two actors being closer in a clustering sense does not necessarily mean a higher expected value of their relationship. Thus, clustering might not be interpretable.")
  
  model[["prior"]][["Z.var.mul"]]<-var.mul
  model[["prior"]][["Z.var"]]<-var
  model[["prior"]][["Z.var.df.mul"]]<-var.df.mul
  model[["prior"]][["Z.var.df"]]<-var.df
  model[["prior"]][["Z.mean.var.mul"]]<-mean.var.mul
  model[["prior"]][["Z.mean.var"]]<-mean.var
  model[["prior"]][["Z.pK.mul"]]<-pK.mul
  model[["prior"]][["Z.pK"]]<-pK
  
  if(!("Z.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.var"]]<-model[["prior"]][["Z.var.mul"]]*(network.size(model[["Yg"]])/max(1,model[["G"]]))^(2/model[["d"]])
  if(!("Z.mean.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.mean.var"]]<-model[["prior"]][["Z.mean.var.mul"]]*model[["prior"]][["Z.var"]]*max(1,model[["G"]])^(2/model[["d"]])
  if(!("Z.var.df" %in% names(model[["prior"]]))) model[["prior"]][["Z.var.df"]]<-model[["prior"]][["Z.var.df.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  if(!("Z.pK" %in% names(model[["prior"]]))) model[["prior"]][["Z.pK"]]<-model[["prior"]][["Z.pK.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))  
  
  model
}
