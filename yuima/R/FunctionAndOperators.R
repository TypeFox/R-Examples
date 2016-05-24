# Method for construction of function and operator of yuima
# object
setMaps <- function(func, yuima, out.var = "", nrow =1 ,ncol=1){
  # A function has three kind of inputs
  # parameters that is a scalar
  # Process that is an object of class yuima
  # Time that is an object of sample grid
  res <- aux.setMaps(func, yuima, out.var = out.var, nrow =1 ,
              ncol=1, type="Maps")
  return(res)
#   if(missing(yuima)){
#     yuima.stop("yuima object is missing.")
#   }
#
#   if(missing(func)){
#     yuima.stop("function is missing.")
#     return(NULL)
#   }
#
# #   if(is.array(func)){
# #     dimens<-dim(func)
# #   }else{
# #     if(length(func)!=(nrow*ncol)){
# #       yuima.warn("nrow*ncol is different from the dim of image. f becomes a vector function")
# #       func<-as.matrix(func)
# #       dimens<-dim(func)
# #     }else{
# #       func<-matrix(func,nrow = nrow, ncol = ncol)
# #       dimens<-dim(func)
# #     }
# #   }
#
#   resFunc<-constFunc(func, nrow, ncol)
#
#   func <- resFunc$func
#   dimens <- resFunc$dimens
#
# #   if(is(yuima, "yuima.model")){
# #     mod<-yuima
# #     yuima<-setYuima(model = mod)
# #   }else{
# #     if(is(yuima, "yuima")){
# #       mod<-yuima@model
# #     }else{
# #       yuima.stop("yuima must be an object of class yuima or yuima.model")
# #     }
# #   }
#
#
#   modDum <- ExtYuimaMod(yuima)
#   mod <- modDum$mod
#   yuima <- modDum$yuima
#
#   paramfunc<-NULL
#   ddd<-prod(dimens)
#   funcList<-as.list(character(length=ddd))
#   func<-as.character(func)
#   for(i in c(1:ddd)){
#     funcList[[i]]<-parse(text=func[i])
#     paramfunc<-c(paramfunc,all.vars(funcList[[i]]))
#   }
# #  funcList<-array(funcList,dim=dimens)
# #   for(j in c(1:ncol)){
# #     for(i in c(1:nrow)){
# #       funcList[[i+(j-1)*nrow]]<-parse(text = func[i,j])
# #       paramfunc<-c(paramfunc,all.vars(funcList[[i+(j-1)*nrow]]))
# #     }
# #   }
#   paramfunc<-unique(paramfunc)
#   common<-mod@parameter@common
#
#   Cond<-(mod@parameter@all %in% paramfunc)
#   common <- c(common,mod@parameter@all[Cond])
#   Cond <- (paramfunc %in% mod@solve.variable)
#   if(sum(Cond)==0){
#     yuima.warn("function does not depend on solve.variable")
#   }
#   paramfunc<-paramfunc[!Cond]
#
#   Cond <- (paramfunc %in% mod@time.variable)
#   paramfunc <- paramfunc[!Cond]
#   if(length(out.var)==1){
#     out.var<-rep(out.var,ddd)
#   }
#   param <- new("param.Output",
#                out.var = out.var,
#                allparam = unique(c(paramfunc,mod@parameter@all)),
#                allparamMap = paramfunc,
#                common = common,
#                Input.var = mod@solve.variable,
#                time.var=mod@time.variable)
#
#   objFunc <- new("info.Output", formula = funcList,
#                  dimension=dimens, type ="Maps")
#
#   res<-new("yuima.Output",
#            param = param,
#            Output = objFunc,
#            yuima=yuima )
#
#   return(res)
}

aux.setMaps <- function(func, yuima, out.var = "",
                        nrow =1 ,ncol=1, type="Maps"){
  if(missing(yuima)){
    yuima.stop("yuima object is missing.")
  }

  if(missing(func)){
    yuima.stop("function is missing.")
    return(NULL)
  }

  #   if(is.array(func)){
  #     dimens<-dim(func)
  #   }else{
  #     if(length(func)!=(nrow*ncol)){
  #       yuima.warn("nrow*ncol is different from the dim of image. f becomes a vector function")
  #       func<-as.matrix(func)
  #       dimens<-dim(func)
  #     }else{
  #       func<-matrix(func,nrow = nrow, ncol = ncol)
  #       dimens<-dim(func)
  #     }
  #   }

  resFunc<-constFunc(func, nrow, ncol)

  func <- resFunc$func
  dimens <- resFunc$dimens

  #   if(is(yuima, "yuima.model")){
  #     mod<-yuima
  #     yuima<-setYuima(model = mod)
  #   }else{
  #     if(is(yuima, "yuima")){
  #       mod<-yuima@model
  #     }else{
  #       yuima.stop("yuima must be an object of class yuima or yuima.model")
  #     }
  #   }


  modDum <- ExtYuimaMod(yuima)
  mod <- modDum$mod
  yuima <- modDum$yuima

  paramfunc<-NULL
  ddd<-prod(dimens)
  funcList<-as.list(character(length=ddd))
  func<-as.character(func)
  for(i in c(1:ddd)){
    funcList[[i]]<-parse(text=func[i])
    paramfunc<-c(paramfunc,all.vars(funcList[[i]]))
  }
  #  funcList<-array(funcList,dim=dimens)
  #   for(j in c(1:ncol)){
  #     for(i in c(1:nrow)){
  #       funcList[[i+(j-1)*nrow]]<-parse(text = func[i,j])
  #       paramfunc<-c(paramfunc,all.vars(funcList[[i+(j-1)*nrow]]))
  #     }
  #   }
  paramfunc<-unique(paramfunc)
  common<-mod@parameter@common

  Cond<-(mod@parameter@all %in% paramfunc)
  common <- c(common,mod@parameter@all[Cond])
  Cond <- (paramfunc %in% mod@solve.variable)
  if(sum(Cond)==0){
    yuima.warn("function does not depend on solve.variable")
  }
  paramfunc<-paramfunc[!Cond]

  Cond <- (paramfunc %in% mod@time.variable)
  paramfunc <- paramfunc[!Cond]
  if(length(out.var)==1){
    out.var<-rep(out.var,ddd)
  }
  param <- new("param.Output",
               out.var = out.var,
               allparam = unique(c(paramfunc,mod@parameter@all)),
               allparamMap = paramfunc,
               common = common,
               Input.var = mod@solve.variable,
               time.var=mod@time.variable)

  objFunc <- new("info.Output", formula = funcList,
                 dimension=dimens, type = type,
                 param=param)

  res<-new("yuima.Output",
           Output = objFunc,
           yuima=yuima )

  return(res)
}


# setIntegral <- function(yuima, integrand, var.dx,
#   lower.var, upper.var, out.var = "", nrow =1 ,ncol=1,
#   type = "Integral"){
#   if(missing(yuima)){
#     yuima.stop("yuima object is missing.")
#   }
#   if(missing(integrand)){
#     yuima.stop("Integrand function is missing")
#   }
#   if(missing(var.dx)){
#     yuima.stop("dx object is missing.")
#   }
#
#   resFunc<-constFunc(func=integrand, nrow, ncol)
#   Integrand <- resFunc$func
#   dimension <- resFunc$dimens
#
#   modDum <- ExtYuimaMod(yuima)
#   mod <- modDum$mod
#   yuima <- modDum$yuima
#   paramIntegrand <- NULL
#   ddd <- prod(dimension)
#   IntegrandList <- as.list(character(length=ddd))
#   Integrand <- as.character(Integrand)
#
#   for(i in c(1:ddd)){
#     IntegrandList[[i]]<-parse(text=Integrand[i])
#     paramIntegrand<-c(paramIntegrand,all.vars(IntegrandList[[i]]))
#   }
#
#   paramIntegrand<-unique(paramIntegrand)
#   common<-mod@parameter@common
#
#   Cond<-(mod@parameter@all %in% paramIntegrand)
#   common <- c(common,mod@parameter@all[Cond])
#   # solve variable
#   Cond <- (paramIntegrand %in% mod@solve.variable)
#   if(sum(Cond)==0){
#     yuima.warn("Integrand fuction does not depend on solve.variable")
#   }
#
#   paramIntegrand <- paramIntegrand[!Cond]
#   # time variable
#   Cond <- (paramIntegrand %in% mod@time.variable)
#   paramIntegrand <- paramIntegrand[!Cond]
#   # upper.var
#   if((upper.var == mod@time.variable)||(lower.var == mod@time.variable)){
#     yuima.stop("upper.var or lower.var must be different from time.variable")
#   }
#
#   Cond <- (paramIntegrand %in% upper.var)
#   paramIntegrand <- paramIntegrand[!Cond]
#
#   Cond <- (paramIntegrand %in% lower.var)
#   paramIntegrand <- paramIntegrand[!Cond]
#
#   allparam <- c(mod@parameter@all, unique(paramIntegrand))
#
#   if(type == "Integral"){
#     cond1 <-c(var.dx %in% c(mod@solve.variable, mod@time.variable))
#     if(sum(cond1)!=dimension[2]){
#       yuima.stop("var.dx must be contains only components of solve variable or time variable")
#     }
#   }
#   my.param.Integral <- new("param.Integral",
#                            allparam = unique(allparam),
#                            common = common,
#                            Integrandparam = paramIntegrand)
#   my.variable.Integral <- new("variable.Integral",
#                               var.dx = var.dx,
#                               lower.var = lower.var,
#                               upper.var = upper.var,
#                               out.var = out.var,
#                               var.time = yuima@model@time.variable)
#   my.integrand <- new("Integrand",
#                       IntegrandList=IntegrandList,
#                       dimIntegrand = dimension)
#
#   my.Integral<-new("Integral.sde",
#                    param.Integral = my.param.Integral,
#                    variable.Integral = my.variable.Integral,
#                    Integrand = my.integrand)
#   res<-new("yuima.Integral",Integral=my.Integral, yuima=yuima)
#   return(res)
#
# #   param <- list(allparam=unique(allparam), common=common,
# #     IntegrandParam = paramIntegrand)
# #
# #   return(list(param = param, IntegrandList=IntegrandList,
# #     var.dx=var.dx, lower.var=lower.var, upper.var=upper.var,
# #     out.var=out.var, dimIntegrand = dimension))
# }
#
# setOperator <- function(operator, X, Y,
#   out.var = "", nrow =1 ,ncol=1){
#   if(is(X, "yuima.model")&& is(Y, "yuima.model")){
#     modtot <- rbind(X,Y)
#   }
#   #assign("mod1",mod1)
#   Oper<- strsplit(operator,split="")[[1]]
#   if(mod1@equation.number!=mod2@equation.number){
#     yuima.stop("the models must have the same dimension")
#   }
#   func <- matrix(character(),mod1@equation.number,1)
#   condX <- (Oper %in% "X")
#   condY <- (Oper %in% "Y")
#   for(i in c(1:mod1@equation.number)){
#     dummyCond <- Oper
#     dummyCond[condX] <- X@solve.variable[i]
#     dummyCond[condY] <- Y@solve.variable[i]
#     func[i,] <- paste0(dummyCond,collapse ="")
#   }
# #   res <- setMaps(func = func, yuima = modtot,
# #     out.var = out.var, nrow = nrow , ncol = ncol)
#    res <- aux.setMaps(func = func, yuima = modtot,
#     out.var = out.var, nrow = nrow ,
#     ncol=ncol, type="Operator")
#   return(res)
# }
#
# setIntensity <- function(...){
#   return(NULL)
# }

constFunc<-function(func, nrow, ncol){
  if(is.array(func)){
    dimens<-dim(func)
  }else{
    if(length(func)!=(nrow*ncol)){
      yuima.warn("nrow*ncol is different from the dim of image. f becomes a vector function")
      func<-as.matrix(func)
      dimens<-dim(func)
    }else{
      func<-matrix(func,nrow = nrow, ncol = ncol)
      dimens<-dim(func)
    }
  }
  return(list(func=func, dimens = dimens))
}

ExtYuimaMod <- function(yuima){
  if(is(yuima, "yuima.model")){
    mod<-yuima
    yuima<-setYuima(model = mod)
  }else{
    if(is(yuima, "yuima")){
      mod<-yuima@model
    }else{
      yuima.stop("yuima must be an object of class yuima or yuima.model")
    }
  }
  return(list(mod=mod, yuima=yuima))
}

