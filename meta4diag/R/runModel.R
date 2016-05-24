runModel <- function(outdata, outpriors, link="logit", quantiles = c(0.025, 0.5, 0.975), verbose=FALSE){
  if(requireNamespace("INLA", quietly = TRUE)){
    model.type = outdata$model.type
    N = dim(outdata$internaldata)[1]
    varnames = names(outdata$internaldata)
    nnv = c("studynames", "Y", "id", "Ntrials") 
    nnvn = unlist(lapply(nnv, function(x) which(varnames==x)))
    varnames = varnames[-nnvn]
    
    
    if(outpriors$wishart.flag){
      data_temp = outdata$internaldata
      data <- rbind(data_temp[seq(1,N,by=2),], data_temp[seq(2,N,by=2),])
      data$id <- 1:N
    }else{
      data = outdata$internaldata
    }
    #data = outdata$internaldata
    
    if(model.type==1){names.fitted = c("Se","Sp")}
    if(model.type==2){names.fitted = c("Se","1-Sp")}
    if(model.type==3){names.fitted = c("1-Se","Sp")}
    if(model.type==4){names.fitted = c("1-Se","1-Sp")}
    
    
    if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        
        names.summarized.fitted = paste("mean(",link,".",names.fitted,")",sep="")
        
        lc1.ind  = agrep("mu", varnames, max.distance=0)
        lc2.ind  = agrep("nu", varnames, max.distance=0)
        
        lc1text = paste("INLA::inla.make.lincomb(",paste(varnames[lc1.ind], "=1",sep="", collapse = ", "),")",sep="")
        lc2text = paste("INLA::inla.make.lincomb(",paste(varnames[lc2.ind], "=1",sep="", collapse = ", "),")",sep="")
        
        lc1 = eval(parse(text=lc1text))
        names(lc1) = names.summarized.fitted[1]
        lc2 = eval(parse(text=lc2text))
        names(lc2) = names.summarized.fitted[2]
        lc = c(lc1, lc2) 
      } else{
        studynames = outdata$originaldata$studynames
        
        mu.ind = which(data[,"mu"]==1)
        nu.ind = which(data[,"nu"]==1)
        
        lc1.ind  = c(agrep("mu", varnames, max.distance=0), agrep("alpha", varnames, max.distance=0))
        lc2.ind  = c(agrep("nu", varnames, max.distance=0), agrep("beta", varnames, max.distance=0))
        lc1text = paste("INLA::inla.make.lincombs(",paste(varnames[lc1.ind], "=", data[mu.ind,varnames[lc1.ind]],sep="", collapse = ", "),")",sep="")
        lc2text = paste("INLA::inla.make.lincombs(",paste(varnames[lc2.ind], "=", data[nu.ind,varnames[lc2.ind]],sep="", collapse = ", "),")",sep="")
        lc1 = eval(parse(text=lc1text))
        names(lc1) = paste("mean(",link,".",names.fitted[1],".",1:(0.5*N),")",sep="")
        lc2 = eval(parse(text=lc2text))
        names(lc2) = paste("mean(",link,".",names.fitted[2],".",1:(0.5*N),")",sep="")
        lc = c(lc1, lc2)   
      }
    }else{
      um = as.character(unique(outdata$originaldata[,outdata$modality.setting]))
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        lc1 = c()
        lc2 = c()
        for(i in 1:length(um)){
          umname = um[i]
          umname = paste(unlist(strsplit(basename(umname), "[-]")), collapse=".")
          names.summarized.fitted = paste("mean(",link,".",names.fitted,".",umname,")",sep="")
          lc1.ind = agrep(paste("mu.",umname,sep=""), varnames, max.distance=0)
          lc2.ind = agrep(paste("nu.",umname,sep=""), varnames, max.distance=0)
          lc1text_um = paste("INLA::inla.make.lincomb(",paste(varnames[lc1.ind], "=1", sep="", collapse = ", "),")",sep="")
          lc2text_um = paste("INLA::inla.make.lincomb(",paste(varnames[lc2.ind], "=1", sep="", collapse = ", "),")",sep="")
          lc1_um = eval(parse(text=lc1text_um))
          names(lc1_um) = names.summarized.fitted[1]
          lc1 = append(lc1, lc1_um)
          lc2_um = eval(parse(text=lc2text_um))
          names(lc2_um) = names.summarized.fitted[2]
          lc2 = append(lc2, lc2_um)
        }
        lc = c(lc1, lc2) 
      }else{
        lc1 = c()
        lc2 = c()
        for(i in 1:length(um)){
          umname = um[i]
          umname = paste(unlist(strsplit(basename(umname), "[-]")), collapse=".")
          mu.ind = which(data[,paste("mu.",umname,sep="")]==1)
          nu.ind = which(data[,paste("nu.",umname,sep="")]==1)
          
          lc1.ind  = c(agrep(paste("mu.",umname,sep=""), varnames, max.distance=0), agrep("alpha", varnames, max.distance=0))
          lc2.ind  = c(agrep(paste("nu.",umname,sep=""), varnames, max.distance=0), agrep("beta", varnames, max.distance=0))
          lc1text = paste("INLA::inla.make.lincombs(",paste(varnames[lc1.ind], "=", data[mu.ind,varnames[lc1.ind]],sep="", collapse = ", "),")",sep="")
          lc2text = paste("INLA::inla.make.lincombs(",paste(varnames[lc2.ind], "=", data[nu.ind,varnames[lc2.ind]],sep="", collapse = ", "),")",sep="")
          lc1_um = eval(parse(text=lc1text))
          names(lc1_um) = paste("mean(",link,".",names.fitted[1],".",umname,".",1:length(mu.ind),")",sep="")
          lc1 = append(lc1, lc1_um)
          lc2_um = eval(parse(text=lc2text))
          names(lc2_um) = paste("mean(",link,".",names.fitted[2],".",umname,".",1:length(nu.ind),")",sep="")
          lc2 = append(lc2, lc2_um)
        }
        lc = c(lc1, lc2) 
      }
    }
    
    quantiles = c(quantiles, 0.025, 0.5, 0.975)
    quantiles = sort(unique(quantiles))
    
    Ntrials = data$Ntrials
    
    if(!outpriors$wishart.flag){
      prec1 = outpriors$prec1
      prec2 = outpriors$prec2
      cor = outpriors$cor
      
      fm <- as.formula(paste("Y ~ f(id, model=\"2diid\", hyper=list(prec1 = prec1,prec2 = prec2,cor = cor), n=N) + ", 
                             paste(varnames, collapse= " + "), "-1"))
      
      model = INLA::inla(fm, family="binomial", data=data, Ntrials=Ntrials, quantiles=quantiles,
                         verbose=verbose, lincomb = lc,  
                         control.inla=list(strategy="gaussian", 
                                           lincomb.derived.correlation.matrix = TRUE),
                         control.predictor=list(compute=TRUE),
                         control.family = list(link = link),
                         control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE, mlik = TRUE,config=TRUE))
      
      if(!model$ok){
        stop("Something wrong while running model with data! Please set verbose=TRUE to check!!!!")
      }
    }else{
      prec1 = outpriors$prec1
      
      fm <- as.formula(paste("Y ~ f(id, model=\"iid2d\", hyper=list(prec1 = prec1), n=N) + ", 
                             paste(varnames, collapse= " + "), "-1"))
      
      
      
      model = INLA::inla(fm, family="binomial", data=data, Ntrials=Ntrials, quantiles=quantiles,
                         verbose=verbose, lincomb = lc,  
                         control.inla=list(strategy="gaussian", 
                                           lincomb.derived.correlation.matrix = TRUE),
                         control.predictor=list(compute=TRUE),
                         control.family = list(link = link),
                         control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE, mlik = TRUE,config=TRUE))
      
      if(!model$ok){
        stop("Something wrong while running model with data! Please set verbose=TRUE to check!!!!")
      } 
    }
    
    model$model.type = model.type
    model$link = link
    model$quantiles = quantiles
    model$verbose = verbose
    model$outdata = data
    model$wishart.flag = outpriors$wishart.flag
    
    if(link=="logit"){
      model$g = function(x){return(log(x/(1-x)))}
      model$inv.g = function(x){return(exp(x)/(1+exp(x)))}
    }
    if(link=="probit"){
      model$g = function(x){return(qnorm(x))}
      model$inv.g = function(x){return(pnorm(x))}
    }
    if(link=="cloglog"){
      model$g = function(x){return(log(-log(1-x)))}
      model$inv.g = function(x){return(1-exp(-exp(x)))}
    }
    
    return(model)
  }else{
    stop("INLA need to be installed and loaded!\n
       Please use the following commants to install and load INLA,\n
       install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/testing\") \n
       library(INLA) \n")
  }
}