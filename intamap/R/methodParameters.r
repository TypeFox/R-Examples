
methodParameters.copula = function(object) {
methodParameters = " "
  
copulaParams = object$copulaParams

methodParameters = paste(methodParameters,"copulaPar = list() \n")

if ("margin" %in% names(copulaParams)) {
  mparams = copulaParams$margin$params
  mmparams = paste(mparams,collapse = ",")
  mmparams = paste("c(",mmparams,")")
  mmparamsName = paste(names(mparams),collapse="\",\"")
  mmparamsName = paste("c(\"",mmparamsName,"\")")
  if (!is.null(names(mparams))) mmparams = paste(mmparams,
      "\n names(margin$params) = ",mmparamsName)
  mname = copulaParams$margin$name
  mlower = paste(copulaParams$margin$lower,collapse=",")
  mlower = paste("c(",mlower,")")
  mupper = paste(copulaParams$margin$upper,collapse=",")
  mupper = paste("c(",mupper,")")
  marginParameters = paste("margin = list()\n  margin$params =", mmparams,
                  "\n margin$name = \"",mname,
                  "\"\n margin$lower = ",mlower,
                  "\n margin$upper = ",mupper,
                  "\n copulaPar$margin = margin \n", sep="")
  methodParameters = paste(methodParameters,marginParameters)
}

if ("correlation" %in% names(copulaParams)) {
  cmodel = copulaParams$correlation$model
  cparams = copulaParams$correlation$params
  ccparams = paste(cparams,collapse=",")
  ccparams = paste("c(",ccparams,")")
  ccparamsName = paste(names(cparams),collapse="\",\"")
  ccparamsName = paste("\"",ccparamsName,"\"")
  if (!is.null(names(cparams))) ccparams = paste(cparams,
      "\n names(correlation$params) = ",ccparamsName)
  clower = paste(copulaParams$correlation$lower,collapse=",")
  clower = paste("c(",clower,")")
  cupper = paste(copulaParams$correlation$upper,collapse=",")
  cupper = paste("c(",cupper,")")
  correlationParameters = paste("correlation = list()\n correlation$model = \"",cmodel,
                  "\"\n correlation$params = ", ccparams,
                  "\n correlation$lower = ",clower,
                  "\n correlation$upper = ",cupper,
                  "\n copulaPar$correlation = correlation \n", sep="")
  methodParameters = paste(methodParameters,correlationParameters)
}

if ("anisotropy" %in% names(copulaParams)) {
  aparams = paste(copulaParams$anisotropy$params,collapse=",")
  aparams = paste("c(",aparams,")")
  alower = paste(copulaParams$anisotropy$lower,collapse=",")
  alower = paste("c(",alower,")")
  aupper = paste(copulaParams$anisotropy$upper,collapse=",")
  aupper = paste("c(",aupper,")")
  anisotropyParameters = paste("anisotropy = list()\n anisotropy$lower = ",alower,
                  "\n anisotropy$upper = ", aupper,
                  "\n anisotropy$params =", aparams,
                  "\n copulaPar$anisotropy = anisotropy \n", sep="")
  methodParameters = paste(methodParameters,anisotropyParameters)
}

if ("trend" %in% names(copulaParams)) {
  tF = paste(copulaParams$trend$F,collapse=",")
  tF = paste("as.matrix(c(",tF,"))")
  tFd = paste(dimnames(copulaParams$trend$F)[[1]],collapse="\",\"")
  tFd = paste("as.list( c(\"",tFd,"\"))",sep="")
  tparams = copulaParams$trend$params
  ttparams = paste(tparams,collapse=",")
  ttparams = paste("c(",ttparams,")")
  ttparamsName = paste(names(tparams),collapse="\",\"")
  ttparamsName = paste("\"",ttparamsName,"\"")
  if (!is.null(names(tparams))) ttparams = paste(tparams,
      "\n names(correlation$params) = ",ttparamsName)
  trendParameters = paste("trend = list() \n trend$F = ",tF,
                  "\n rownames(trend$F) = ",tFd,
                  "\n trend$params = ",ttparams,
                  "\n trend$lower = ",copulaParams$trend$lower,
                  "\n trend$upper = ",copulaParams$trend$upper,
                  "\n copulaPar$trend = trend \n", sep="")
  methodParameters = paste(methodParameters,trendParameters)
}

copulaP = paste("copula = list(method = \"",copulaParams$copula,"\") \n",
                 "copulaPar$copula = copula \n", sep="")
methodParameters = paste(methodParameters,copulaP)
methodParameters = paste(methodParameters,"object$copulaParams = copulaPar \n")

object$methodParameters = methodParameters
object = NextMethod(object)
object
}




methodParameters.default = function(object) {
  if ("methodParameters" %in% names(object)) {
    methodParameters = object$methodParameters
  } else methodParameters = " "
  if ("variogramModel" %in% names(object)) {
    vmodel = object$variogramModel
    if (is(vmodel,"variogramModel")) {
      vmodel = object$variogramModel
      nvar = dim(vmodel)[1]
      ncols = dim(vmodel)[2]
      mp = paste("vmodel = data.frame(matrix(0,nrow = ",nvar,",ncol = ",ncols,"))\n")
      mpName = paste(names(vmodel),collapse="\",\"")
      mpName = paste("\"",mpName,"\"",sep="")
      mp = paste(mp,"names(vmodel) = c(",mpName,")\n",sep="")
      for (nv in 1:nvar) {
        mp = paste(mp,"vmodel[",nv,",1] = \"",vmodel$model[nv],"\"\n",sep="")
        for (nc in 2:ncols) {
          mp = paste(mp,"vmodel[",nv,",",nc,"] = ",vmodel[nv,nc],"\n",sep="")
        }
      }
      mp = paste(mp,"class(vmodel) =  c(\"variogramModel\",\"data.frame\") \n")
      methodParameters = paste(methodParameters,mp,"\n object$variogramModel = vmodel \n")
    } else warning(paste("Not able to create methodParameters for variogram of class",class(vmodel)))
    if ("lambda" %in% names(object)) methodParameters =    
      paste(methodParameters,"object$lambda = ", object$lambda, "\n",sep="")
  }
  if ("anisPar" %in% names(object)) {
    anisPar = object$anisPar
    QQt = paste(anisPar$Q,collapse=",")
    QQt = paste("matrix(c(",QQt,"),ncol=3)")
    QQd = dimnames(anisPar$Q)[[2]]
    if (is.null(QQd)) {
      QQd = paste("c(1:",length(anisPar$Q),")" )
    } else {
      QQd = paste(QQd,collapse="\",\"")
      QQd = paste("list(\"",QQd,"\")",sep="")
    }
    am = paste("QQ = ",QQt,"\n colnames(QQ) =", QQd, "\n",
               "anisPar = list(ratio = ",anisPar$ratio,
               ", direction = ",anisPar$direction,
               ", Q= QQ,doRotation = ",anisPar$doRotation,")\n")
    methodParameters = paste(methodParameters,am)
    methodParameters = paste(methodParameters,"object$anisPar = anisPar \n")
  }
#  if ("TGcorrection" %in% names(object)) {
#    methodParameters = paste(methodParameters,"object$TGcorrection = ",object$TGcorrection,"\n")
#  }
  methodParameters = paste(methodParameters,"class(object) = \"",class(object), "\"\n",sep="")
  if (!is.null(object$params$set.seed)) methodParameters = 
     paste(methodParameters,"object$params$set.seed = ", object$params$set.seed, "\n",sep="")
  if (!is.null(object$params$testMean)) methodParameters = 
     paste(methodParameters,"object$params$testMean = ", object$params$testMean, "\n",sep="")
  object$methodParameters = methodParameters
  object
}

 
    
    

methodParameters.idw = function(object) {
  methodParameters = paste("object$inverseDistancePower = ",object$inverseDistancePower,"\n")
  object$methodParameters = methodParameters
  object = NextMethod(object)
  object
}

          
          
          
          
methodParameters.default.old = function(object) {
  if ("methodParameters" %in% names(object)) {
    methodParameters = object$methodParameters
  } else methodParameters = " "
  if ("variogramModel" %in% names(object)) {
    vmodel = object$variogramModel
    for (i in dim(vmodel)[1] : 1) {
      psill = vmodel$psill[i]
      vmod = vmodel$model[i]
      range = vmodel$range[i]
      anis = paste(vmodel$ang1[i],vmodel$anis1[i],sep=",")
      
      mp = paste("vgm(",psill,",\"",vmod,"\",",range,",anis =c(",anis,")",sep="")
      if (i == dim(vmodel)[1]) {
        methodParameters = paste(methodParameters,"\n vvmod = ",mp)
      } else {
        methodParameters = paste(methodParameters,",add.to = ",mp)
      }
    }
    for (i in 1:dim(vmodel)[1]) methodParameters = paste(methodParameters,")")
    if ("kappa" %in% names(vmodel)) {
      kappaPar = paste(vmodel$kappa,collapse=",")
      kappaPar = paste("c(",kappaPar,")")      
      methodParameters = paste(methodParameters,
          "\n vvmod$kappa = ", kappaPar,"\n")
    }
    if ("beta" %in% names(vmodel)) methodParameters = paste(methodParameters,
          "\n vvmod$beta = ",vmodel$beta[1],"\n",sep="")
    methodParameters = paste(methodParameters,"\n object$variogramModel = vvmod \n")
  }
  if ("anisPar" %in% names(object)) {
    anisPar = object$anisPar
    QQt = paste(anisPar$Q,collapse=",")
    QQt = paste("matrix(c(",QQt,"),ncol=3)")
    QQd = dimnames(anisPar$Q)[[2]]
    if (is.null(QQd)) {
      QQd = paste("c(1:",length(anisPar$Q),")" )
    } else {
      QQd = paste(QQd,collapse="\",\"")
      QQd = paste("list(\"",QQd,"\")",sep="")
    }
    am = paste("QQ = ",QQt,"\n colnames(QQ) =", QQd, "\n",
               "anisPar = list(ratio = ",anisPar$ratio,
               ", direction = ",anisPar$direction,
               ", Q= QQ,doRotation = ",anisPar$doRotation,")\n")
    methodParameters = paste(methodParameters,am)
    methodParameters = paste(methodParameters,"object$anisPar = anisPar \n")
  }
  methodParameters = paste(methodParameters,"class(object) = \"",class(object),"\"\n",sep="")
  object$methodParameters = methodParameters
  object
}
