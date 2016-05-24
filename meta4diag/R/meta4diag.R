meta4diag = function(data=NULL, model.type = 1, 
                     var.prior = "invgamma", var2.prior = "invgamma", cor.prior = "normal",
                     var.par = c(0.25, 0.025), var2.par, cor.par = c(0,5),
                     wishart.par = c(4, 1, 1, 0),
                     init = c(0.01,0.01,0), link="logit", quantiles=c(0.025,0.5,0.975),
                     modality=NULL, covariates = NULL,
                     verbose = FALSE, nsample=FALSE){
  
  var.prior = tolower(var.prior)
  var2.prior = tolower(var2.prior)
  cor.prior = tolower(cor.prior)
   
  ################ check data
  if(!is.data.frame(data)){
    stop("Data MUST be a data frame!!!")
  }
  ################ check model.type
  if(length(model.type)!=1){
    stop("Argument \"model.type\" can ONLY be ONE integer of c(1,2,3,4)!!!")
  }
  if(!is.numeric(model.type)){
    stop("Argument \"model.type\" can ONLY be ONE integer of c(1,2,3,4)!!!")
  }else{ # model.type is numerical
    if(!(model.type %in% c(1,2,3,4))){
      stop("Argument \"model.type\" can ONLY be ONE integer of c(1,2,3,4)!!!")
    }
  }
  
  I = dim(data)[1]
  ################
  variables.names = colnames(data)
  if(!("studynames" %in% variables.names)){
    data$studynames = paste("study[",c(1:I),"]",sep="")
  }
  
  ################ check link
  if(!is.character(link)){
    stop("Argument \"link\" can ONLY be character. The options are \"logit\", \"probit\" and \"cloglog\"!!!")
  }
  if(length(link)!=1){
    stop("Argument \"link\" can ONLY be character. The options are \"logit\", \"probit\" and \"cloglog\"!!!")
  }
  ################ check verbose
  if(!is.logical(verbose)){
    stop("Argument \"verbose\" can ONLY be logic, either \"TRUE\" or \"FALSE\"!!!")
  }
  
  if(any(c(var.prior, var2.prior, cor.prior)=="invwishart")){
    var.prior = var2.prior = cor.prior = "invwishart"
  }else{
    if(missing(var2.par)){
      if(var2.prior==var.prior){
        var2.par = var.par
      }else{
        stop("Please give the parameters of the prior for second variance component!")
      }
    }
  }

  ################ Make prior, and in the makePrior function, check var.prior, var.par, cor.prior, cor.par and init
  outpriors = makePriors(var.prior=var.prior, var2.prior=var2.prior, cor.prior=cor.prior, 
                         var.par=var.par, var2.par=var2.par, cor.par=cor.par, init=init)
  
  ################ Make data, and in the makedata function, check covariates and compare
  outdata = makeData(data = data, model.type = model.type, modality = modality, covariates = covariates)
  
  ################ Run model in INLA
  model = runModel(outdata=outdata, outpriors=outpriors, link=link, quantiles=quantiles, verbose = verbose)
  
  ##########################  construct the result
  res = makeObject(outdata, outpriors, model, nsample=nsample)
  
  return(res)
}








# .ROC = function(x, add=FALSE){
#   modelnames = x$names.model
#   lm = length(modelnames)
#   is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#   if(length(model)!=1){stop("Argument \"model\" should be a inteeger.")}
#   if(!is.numeric(model)){
#     stop(paste("Argument \"model\" should be in c(1:",lm,").",sep=""))
#   }
#   if(!is.wholenumber(model)){
#     stop(paste("Argument \"model\" should be in c(1:",lm,").",sep=""))
#   }
#   fitname = x$names.fitted
#   fullname = paste("summary.predict.",fitname,sep="")
#   fitfullname = paste("summary.fitted.",fitname,sep="")
#   t = seq(0, 2*pi, by = 2*pi/100)
#   par(mfrow=c(2,2),mar=c(5.1,5.1,2.1,1.1))
#   for(i in 1:lm){
#     mean.A = x[[fullname[1]]][[i]][,1]
#     sd.A = x[[fullname[1]]][[i]][,2]
#     mean.B = x[[fullname[2]]][[i]][,1]
#     sd.B = x[[fullname[2]]][[i]][,2]
#     r = x$mean.correlation[i]
#     
#     I = length(mean.A)
#     f = qf(0.95, 2, I-2)
#     c = sqrt(2*f)
#     
#     A = mean.A[1] + sd.A[1]*c*cos(t)
#     B = mean.B[1] + sd.B[1]*c*cos(t + acos(r))
#     confidence.A = .invlogit(A)
#     confidence.B = .invlogit(B)
#     plot(confidence.B, confidence.A, type="l",xlim=c(0,1),ylim=c(0,1),xlab=fitname[2],ylab=fitname[1])
#     points(x[[fitfullname[2]]][[i]][,1],x[[fitfullname[1]]][[i]][,1],pch=1)
#     for(j in 2:I){
#       A = mean.A[j] + sd.A[j]*c*cos(t)
#       B = mean.B[j] + sd.B[j]*c*cos(t + acos(r))
#       confidence.A = .invlogit(A)
#       confidence.B = .invlogit(B)
#       lines(confidence.B, confidence.A)
#     }
#   }
# }