makeObject <- function(outdata, outpriors, model, nsample=FALSE){
  if(requireNamespace("INLA", quietly = TRUE)){
    quantiles = sort(unique(c(model$quantiles,0.025,0.5,0.975)))
    effect.length = length(quantiles) + 2
    est = list()
    est$data = outdata$originaldata
    est$outdata = model$outdata
    colnames(est$data) = tolower(colnames(est$data))
    est$priors.density = outpriors$density
    
    
    ########## model type
    if(model$model.type==1){
      est$names.fitted = c("Se","Sp")
      est$names.transf.fitted = c("1-Se","1-Sp")
    }
    if(model$model.type==2){
      est$names.fitted = c("Se","1-Sp")
      est$names.transf.fitted = c("1-Se","Sp")
    }
    if(model$model.type==3){
      est$names.fitted = c("1-Se","Sp")
      est$names.transf.fitted = c("Se","1-Sp")
    }
    if(model$model.type==4){
      est$names.fitted = c("1-Se","1-Sp")
      est$names.transf.fitted = c("Se","Sp")
    }
    
    n.marginal.point = dim(model$marginals.fixed[[1]])[1]
    ########## model cpu & call
    est$cpu.used = model[["cpu.used"]]
    est$call = model[["call"]]
    
    ######## fixed
    est$summary.fixed = model[["summary.fixed"]][,c(1:effect.length)]
    est$marginals.fixed = model[["marginals.fixed"]]
    
    ############# hyperpar: need transformation
    names.var = c("var_phi","var_psi")
    names.cor = "cor"
    #names.var = paste("var(",model$link,"s.",est$names.fitted,")",sep="")
    #names.cor = paste("cor(",model$link,"s)",sep="")
    tau1.marginals = model[["marginals.hyperpar"]][[1]]
    var1.marginals = INLA::inla.tmarginal(function(x) 1/x, tau1.marginals, n=1024)
    tau2.marginals = model[["marginals.hyperpar"]][[2]]
    var2.marginals = INLA::inla.tmarginal(function(x) 1/x, tau2.marginals, n=1024)  
    marginals.hyperpar.temp = list()
    marginals.hyperpar.temp[[names.var[1]]] = var1.marginals
    marginals.hyperpar.temp[[names.var[2]]] = var2.marginals
    marginals.hyperpar.temp[[names.cor]] = model[["marginals.hyperpar"]][[3]]
    est$marginals.hyperpar = marginals.hyperpar.temp
    
    var1.marginals = est$marginals.hyperpar[[1]]
    var2.marginals = est$marginals.hyperpar[[2]]
    suminfo1 = .summary.marginal(var1.marginals,level=quantiles)
    suminfo2 = .summary.marginal(var2.marginals,level=quantiles)
    
    summary.hyperpar.temp = rbind(suminfo1, suminfo2, model[["summary.hyperpar"]][3,c(1:effect.length)])
    summary.hyperpar.temp = as.matrix(summary.hyperpar.temp)
    rownames(summary.hyperpar.temp) = c(names.var,names.cor)
    est$summary.hyperpar = summary.hyperpar.temp
    
    
    
    ######## study names
    studynames = as.character(outdata$originaldata$studynames)
    ######## Expected g=="logit","probit","cloglog" E(g(..)) and E(g(..))
    if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        est$summary.expected.gtransformed.accuracy = model[["summary.lincomb.derived"]][,c(2:(effect.length+1))]
        est$marginals.expected.gtransformed.accuracy = model[["marginals.lincomb.derived"]]
        est$correlation.matrix.linear.comb = model[["misc"]]$lincomb.derived.correlation.matrix
        est$covariance.matrix.linear.comb = model[["misc"]]$lincomb.derived.covariance.matrix
        est$correlation.linear.comb = model[["misc"]]$lincomb.derived.correlation.matrix[1,2]
        names(est$correlation.linear.comb) = "Cor(mu, nu)"
        est$covariance.linear.comb = model[["misc"]]$lincomb.derived.covariance.matrix[1,2]
        names(est$covariance.linear.comb) = "Cov(mu, nu)"
      } else{
        names.summarized.fixed1 = paste("mean(",model$link,".",est$names.fitted[1],".",studynames,")",sep="")
        names.summarized.fixed2 = paste("mean(",model$link,".",est$names.fitted[2],".",studynames,")",sep="")
        names.summarized.fixed = c(names.summarized.fixed1, names.summarized.fixed2)
        est$summary.expected.gtransformed.accuracy = as.matrix(model[["summary.lincomb.derived"]][,c(2:(effect.length+1))])
        est$marginals.expected.gtransformed.accuracy = model[["marginals.lincomb.derived"]]
        rownames(est$summary.expected.gtransformed.accuracy) = names.summarized.fixed
        names(est$marginals.expected.gtransformed.accuracy) = names.summarized.fixed
        
        est$correlation.linear.comb = do.call(rbind, lapply(1:length(studynames), 
            function(i) model[["misc"]]$lincomb.derived.correlation.matrix[i, i+length(studynames)]))
        colnames(est$correlation.linear.comb) = paste("Cor(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$correlation.linear.comb = as.matrix(est$correlation.linear.comb)
        rownames(est$correlation.linear.comb) = studynames
        
        est$covariance.linear.comb = as.matrix(do.call(rbind, lapply(1:length(studynames), 
            function(i) model[["misc"]]$lincomb.derived.covariance.matrix[i, i+length(studynames)])))
        colnames(est$covariance.linear.comb) = paste("Cov(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$covariance.linear.comb = as.matrix(est$covariance.linear.comb)
        rownames(est$covariance.linear.comb) = studynames
      }
    }else{
      um = as.character(unique(outdata$originaldata[,outdata$modality.setting]))
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        est$summary.expected.gtransformed.accuracy = model[["summary.lincomb.derived"]][,c(2:(effect.length+1))]
        est$marginals.expected.gtransformed.accuracy = model[["marginals.lincomb.derived"]]
        est$correlation.matrix.linear.comb = model[["misc"]]$lincomb.derived.correlation.matrix
        est$covariance.matrix.linear.comb = model[["misc"]]$lincomb.derived.covariance.matrix
        
        est$correlation.linear.comb = do.call(rbind, lapply(1:length(um), 
                     function(i) model[["misc"]]$lincomb.derived.correlation.matrix[i, i+length(um)]))
        colnames(est$correlation.linear.comb) = paste("Cor(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$correlation.linear.comb = as.matrix(est$correlation.linear.comb)
        rownames(est$correlation.linear.comb) = um
        
        est$covariance.linear.comb = do.call(rbind, lapply(1:length(um), 
                     function(i) model[["misc"]]$lincomb.derived.covariance.matrix[i, i+length(um)]))
        colnames(est$covariance.linear.comb) = paste("Cov(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$covariance.linear.comb = as.matrix(est$covariance.linear.comb)
        rownames(est$covariance.linear.comb) = um
      } else{
        names.summarized.fixed1 = paste("mean(",model$link,".",est$names.fitted[1],".",studynames,")",sep="")
        names.summarized.fixed2 = paste("mean(",model$link,".",est$names.fitted[2],".",studynames,")",sep="")
        names.summarized.fixed = c(names.summarized.fixed1, names.summarized.fixed2)
        est$summary.expected.gtransformed.accuracy = as.matrix(model[["summary.lincomb.derived"]][,c(2:(effect.length+1))])
        est$marginals.expected.gtransformed.accuracy = model[["marginals.lincomb.derived"]]
        rownames(est$summary.expected.gtransformed.accuracy) = names.summarized.fixed
        names(est$marginals.expected.gtransformed.accuracy) = names.summarized.fixed
        
        est$correlation.linear.comb = do.call(rbind, lapply(1:length(studynames), 
            function(i) model[["misc"]]$lincomb.derived.correlation.matrix[i, i+length(studynames)]))
        colnames(est$correlation.linear.comb) = paste("Cor(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$correlation.linear.comb = as.matrix(est$correlation.linear.comb)
        rownames(est$correlation.linear.comb) = studynames
        
        est$covariance.linear.comb = as.matrix(do.call(rbind, lapply(1:length(studynames), 
            function(i) model[["misc"]]$lincomb.derived.covariance.matrix[i, i+length(studynames)])))
        colnames(est$covariance.linear.comb) = paste("Cov(",paste("E[g(",est$names.fitted,")]",sep="", collapse=", "),")",sep="")
        est$covariance.linear.comb = as.matrix(est$covariance.linear.comb)
        rownames(est$covariance.linear.comb) = studynames
      }
    }
    
    ########## Expected accuracy E(..) and E(..)
    if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        ####### original----from model
        names.expected.accuracy.original = paste("mean(",est$names.fitted,")",sep="")
        
        names.temp = names(model[["marginals.lincomb.derived"]])
        marginals.expected.accuracy.original = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.original) = names.expected.accuracy.original
        
        se.marginals = marginals.expected.accuracy.original[[1]]
        sp.marginals = marginals.expected.accuracy.original[[2]]
        suminfo1 = .summary.marginal(se.marginals,level=quantiles)
        suminfo2 = .summary.marginal(sp.marginals,level=quantiles)
        summary.expected.accuracy.original = rbind(suminfo1, suminfo2)
        summary.expected.accuracy.original = as.matrix(summary.expected.accuracy.original)
        rownames(summary.expected.accuracy.original) = names.expected.accuracy.original
        
        ####### transformed----from inla.tmarginal
        names.expected.accuracy.transform = paste("mean(",est$names.transf.fitted,")",sep="")
        
        marginals.expected.accuracy.transform = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) 1-model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.transform) = names.expected.accuracy.transform
        
        se.marginals = marginals.expected.accuracy.transform[[1]]
        sp.marginals = marginals.expected.accuracy.transform[[2]]
        suminfo1 = .summary.marginal(se.marginals,level=quantiles)
        suminfo2 = .summary.marginal(sp.marginals,level=quantiles)
        summary.expected.accuracy.transform = rbind(suminfo1, suminfo2)
        summary.expected.accuracy.transform = as.matrix(summary.expected.accuracy.transform)
        rownames(summary.expected.accuracy.transform) = names.expected.accuracy.transform
        
        est$summary.expected.accuracy = rbind(summary.expected.accuracy.original, summary.expected.accuracy.transform) 
        est$marginals.expected.accuracy = append(marginals.expected.accuracy.original, marginals.expected.accuracy.transform)
      } else{
        ####### original----from model
        names.expected.accuracy.original1 = paste("mean(",est$names.fitted[1],".",studynames,")",sep="")
        names.expected.accuracy.original2 = paste("mean(",est$names.fitted[2],".",studynames,")",sep="")
        names.expected.accuracy.original = c(names.expected.accuracy.original1, names.expected.accuracy.original2)
        
        names.temp = names(model[["marginals.lincomb.derived"]])
        marginals.expected.accuracy.original = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.original) = names.expected.accuracy.original
        
        suminfo = lapply(1:length(marginals.expected.accuracy.original), function(x) .summary.marginal(marginals.expected.accuracy.original[[x]],level=quantiles))
        summary.expected.accuracy.original = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.original) = names.expected.accuracy.original
        ####### transformed----from inla.tmarginal
        names.expected.accuracy.transform1 = paste("mean(",est$names.transf.fitted[1],".",studynames,")",sep="")
        names.expected.accuracy.transform2 = paste("mean(",est$names.transf.fitted[2],".",studynames,")",sep="")
        names.expected.accuracy.transform = c(names.expected.accuracy.transform1, names.expected.accuracy.transform2)
        
        marginals.expected.accuracy.transform = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) 1-model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.transform) = names.expected.accuracy.transform
        
        suminfo = lapply(1:length(marginals.expected.accuracy.transform), function(x) .summary.marginal(marginals.expected.accuracy.transform[[x]],level=quantiles))
        summary.expected.accuracy.transform = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.transform) = names.expected.accuracy.transform
        
        est$summary.expected.accuracy = rbind(summary.expected.accuracy.original, summary.expected.accuracy.transform) 
        est$marginals.expected.accuracy = append(marginals.expected.accuracy.original, marginals.expected.accuracy.transform)
      }
    }else{
      um = as.character(unique(outdata$originaldata[,outdata$modality.setting]))
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        ####### original----from model
        names.expected.accuracy.original = unlist(lapply(est$names.fitted,function(x) paste("mean(",x,".",um,")",sep="")))
        
        names.temp = names(model[["marginals.lincomb.derived"]])
        marginals.expected.accuracy.original = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.original) = names.expected.accuracy.original
        
        suminfo = lapply(1:length(marginals.expected.accuracy.original), function(x) .summary.marginal(marginals.expected.accuracy.original[[x]],level=quantiles))
        summary.expected.accuracy.original = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.original) = names.expected.accuracy.original
        
        ####### transformed----from inla.tmarginal
        names.expected.accuracy.transform = unlist(lapply(est$names.transf.fitted,function(x) paste("mean(",x,".",um,")",sep="")))
        
        marginals.expected.accuracy.transform = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) 1-model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.transform) = names.expected.accuracy.transform
        
        suminfo = lapply(1:length(marginals.expected.accuracy.transform), function(x) .summary.marginal(marginals.expected.accuracy.transform[[x]],level=quantiles))
        summary.expected.accuracy.transform = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.transform) = names.expected.accuracy.transform
        
        est$summary.expected.accuracy = rbind(summary.expected.accuracy.original, summary.expected.accuracy.transform) 
        est$marginals.expected.accuracy = append(marginals.expected.accuracy.original, marginals.expected.accuracy.transform)
      } else{
        ####### original----from model
        names.expected.accuracy.original1 = paste("mean(",est$names.fitted[1],".",studynames,")",sep="")
        names.expected.accuracy.original2 = paste("mean(",est$names.fitted[2],".",studynames,")",sep="")
        names.expected.accuracy.original = c(names.expected.accuracy.original1, names.expected.accuracy.original2)
        
        names.temp = names(model[["marginals.lincomb.derived"]])
        marginals.expected.accuracy.original = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.original) = names.expected.accuracy.original
        
        suminfo = lapply(1:length(marginals.expected.accuracy.original), function(x) .summary.marginal(marginals.expected.accuracy.original[[x]],level=quantiles))
        summary.expected.accuracy.original = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.original) = names.expected.accuracy.original
        ####### transformed----from inla.tmarginal
        names.expected.accuracy.transform1 = paste("mean(",est$names.transf.fitted[1],".",studynames,")",sep="")
        names.expected.accuracy.transform2 = paste("mean(",est$names.transf.fitted[2],".",studynames,")",sep="")
        names.expected.accuracy.transform = c(names.expected.accuracy.transform1, names.expected.accuracy.transform2)
        
        marginals.expected.accuracy.transform = lapply(names.temp, function(y){
          INLA::inla.tmarginal(function(x) 1-model$inv.g(x), model[["marginals.lincomb.derived"]][[y]], n=n.marginal.point)
        })
        names(marginals.expected.accuracy.transform) = names.expected.accuracy.transform
        
        suminfo = lapply(1:length(marginals.expected.accuracy.transform), function(x) .summary.marginal(marginals.expected.accuracy.transform[[x]],level=quantiles))
        summary.expected.accuracy.transform = do.call(rbind, suminfo)
        rownames(summary.expected.accuracy.transform) = names.expected.accuracy.transform
        
        est$summary.expected.accuracy = rbind(summary.expected.accuracy.original, summary.expected.accuracy.transform) 
        est$marginals.expected.accuracy = append(marginals.expected.accuracy.original, marginals.expected.accuracy.transform)
      }
    }
    
    
    
    
    if(outpriors$wishart.flag){
      fitted1.ind = seq(1,dim(est$data)[1],by=1)
      fitted2.ind = seq((dim(est$data)[1]+1),dim(est$outdata)[1],by=1)
    }else{
      fitted1.ind = seq(1,dim(est$outdata)[1],by=2)
      fitted2.ind = seq(2,dim(est$outdata)[1],by=2)
    }
    
    ############# predict
    studynames = as.character(est$outdata$studynames[fitted1.ind])
    summary.predict.temp1 = as.matrix(model[["summary.linear.predictor"]][fitted1.ind,c(1:effect.length)])
    rownames(summary.predict.temp1) = studynames
    summary.predict.temp2 = as.matrix(model[["summary.linear.predictor"]][fitted2.ind,c(1:effect.length)])
    rownames(summary.predict.temp2) = studynames
    est[[paste("summary.predict.(",est$names.fitted[1],")",sep="")]] = summary.predict.temp1
    est[[paste("summary.predict.(",est$names.fitted[2],")",sep="")]] = summary.predict.temp2
    
    marginals.predict.temp1 = lapply(fitted1.ind, function(y){model[["marginals.linear.predictor"]][[y]]})
    names(marginals.predict.temp1) = studynames
    marginals.predict.temp2 = lapply(fitted2.ind, function(y){model[["marginals.linear.predictor"]][[y]]})
    names(marginals.predict.temp2) = studynames
    est[[paste("marginals.predict.(",est$names.fitted[1],")",sep="")]] = marginals.predict.temp1
    est[[paste("marginals.predict.(",est$names.fitted[2],")",sep="")]] = marginals.predict.temp2
    
    ############# transform to other accurary that not fitted here directly
    transfunc = function(x) 1-model$inv.g(x)
    marginals.transf.fitted1 = lapply(marginals.predict.temp1, function(x){INLA::inla.tmarginal(transfunc, x, n=n.marginal.point)})
    marginals.transf.fitted2 = lapply(marginals.predict.temp2, function(x){INLA::inla.tmarginal(transfunc, x, n=n.marginal.point)})
    
    suminfo1 = do.call(rbind,lapply(marginals.transf.fitted1, function(x) .summary.marginal(x,level=quantiles)))
    suminfo2 = do.call(rbind,lapply(marginals.transf.fitted2, function(x) .summary.marginal(x,level=quantiles)))
    
    est[[paste("summary.fitted.(",est$names.transf.fitted[1],")",sep="")]] = suminfo1
    est[[paste("summary.fitted.(",est$names.transf.fitted[2],")",sep="")]] = suminfo2
    
    ############# fitted
    summary.fitted.temp1 = as.matrix(model[["summary.fitted.values"]][fitted1.ind,c(1:effect.length)])
    rownames(summary.fitted.temp1) = studynames
    summary.fitted.temp2 = as.matrix(model[["summary.fitted.values"]][fitted2.ind,c(1:effect.length)])
    rownames(summary.fitted.temp2) = studynames
    est[[paste("summary.fitted.(",est$names.fitted[1],")",sep="")]] = summary.fitted.temp1
    est[[paste("summary.fitted.(",est$names.fitted[2],")",sep="")]] = summary.fitted.temp2
    
    marginals.fitted.temp1 = lapply(fitted1.ind, function(y){model[["marginals.fitted.values"]][[y]]})
    names(marginals.fitted.temp1) = studynames
    marginals.fitted.temp2 = lapply(fitted2.ind, function(y){model[["marginals.fitted.values"]][[y]]})
    names(marginals.fitted.temp2) = studynames
    est[[paste("marginals.fitted.(",est$names.fitted[1],")",sep="")]] = marginals.fitted.temp1
    est[[paste("marginals.fitted.(",est$names.fitted[2],")",sep="")]] = marginals.fitted.temp2
    
    
    ############# samples
    if(is.logical(nsample)){
      if(nsample==FALSE){
        est$misc$sample.flag = FALSE
      }else{
        est$misc$sample.flag = TRUE
        message("Argument \"nsample\" set to TRUE, we will give 5000 samples!")
        est$misc$nsample = 5000
      }
    }else if(is.numeric(nsample)){
      est$misc$sample.flag = TRUE
      if(abs(nsample - round(nsample)) < .Machine$double.eps^0.5){
        est$misc$nsample = nsample
      }else{
        # message("Argument \"nsample\" should be a integer, we round the given number to interger.")
        nsample = round(nsample, digits = 0)
        est$misc$nsample = nsample
      }
    }else{
      message("Argument \"nsample\" should be TRUE, FALSE or a integer, we set it to FALSE!")
    }
    if(est$misc$sample.flag){
      options(warn=-1)
      samples = INLA::inla.posterior.sample(n = nsample, model)
      options(warn=0)
      predictors.samples = do.call(cbind,lapply(c(1:nsample), function(x) samples[[x]]$latent[1:dim(outdata$internaldata)[1]]))
      fixed.samples = do.call(cbind,lapply(c(1:nsample), function(x){
        length.latent = length(samples[[x]]$latent)
        a = samples[[x]]$latent[(length.latent-dim(model$summary.fixed)[1]+1):length.latent]
        return(a)
      }))
      fixed.names = rownames(model$summary.fixed)
      rownames(fixed.samples) = fixed.names
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        ind.fitted1  = agrep("mu", fixed.names, max.distance=0)
        ind.fitted2  = agrep("nu", fixed.names, max.distance=0)
        if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
          mean.logit.fitted1.samples = fixed.samples[ind.fitted1,]
          mean.logit.fitted2.samples = fixed.samples[ind.fitted2,]
          
          mean.fitted1.samples = model$inv.g(mean.logit.fitted1.samples)
          mean.fitted2.samples = model$inv.g(mean.logit.fitted2.samples)
          
          if(model$model.type==1){
            mean.LRpos.samples = mean.fitted1.samples/(1-mean.fitted2.samples)
            mean.LRneg.samples = (1-mean.fitted1.samples)/mean.fitted2.samples
            mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
            mean.RD.samples = mean.fitted1.samples - (1-mean.fitted2.samples)
            mean.LLRpos.samples = log(mean.LRpos.samples)
            mean.LLRneg.samples = log(mean.LRneg.samples)
            mean.LDOR.samples = log(mean.DOR.samples)
          }
          if(model$model.type==2){
            mean.LRpos.samples = mean.fitted1.samples/mean.fitted2.samples
            mean.LRneg.samples = (1-mean.fitted1.samples)/(1-mean.fitted2.samples)
            mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
            mean.RD.samples = mean.fitted1.samples - mean.fitted2.samples
            mean.LLRpos.samples = log(mean.LRpos.samples)
            mean.LLRneg.samples = log(mean.LRneg.samples)
            mean.LDOR.samples = log(mean.DOR.samples)
          }
          if(model$model.type==3){
            mean.LRpos.samples = (1-mean.fitted1.samples)/(1-mean.fitted2.samples)
            mean.LRneg.samples = mean.fitted1.samples/mean.fitted2.samples
            mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
            mean.RD.samples = (1-mean.fitted1.samples) - (1-mean.fitted2.samples)
            mean.LLRpos.samples = log(mean.LRpos.samples)
            mean.LLRneg.samples = log(mean.LRneg.samples)
            mean.LDOR.samples = log(mean.DOR.samples)
          }
          if(model$model.type==4){
            mean.LRpos.samples = (1-mean.fitted1.samples)/mean.fitted2.samples
            mean.LRneg.samples = mean.fitted1.samples/(1-mean.fitted2.samples)
            mean.DOR.samples = mean.LRpos.samples/mean.LRneg.samples
            mean.RD.samples = (1-mean.fitted1.samples) - mean.fitted2.samples
            mean.LLRpos.samples = log(mean.LRpos.samples)
            mean.LLRneg.samples = log(mean.LRneg.samples)
            mean.LDOR.samples = log(mean.DOR.samples)
          }
          
          summary.LRpos.temp = .summary.samples(mean.LRpos.samples, level=quantiles)
          summary.LRneg.temp = .summary.samples(mean.LRneg.samples, level=quantiles)
          summary.DOR.temp = .summary.samples(mean.DOR.samples, level=quantiles)
          summary.RD.temp = .summary.samples(mean.RD.samples, level=quantiles)
          summary.LDOR.temp = .summary.samples(mean.LDOR.samples, level=quantiles)
          summary.LLRpos.temp = .summary.samples(mean.LLRpos.samples, level=quantiles)
          summary.LLRneg.temp = .summary.samples(mean.LLRneg.samples, level=quantiles)
          
          summary.summarized.statistics = rbind(summary.LRpos.temp, summary.LRneg.temp, summary.DOR.temp, summary.RD.temp, summary.LDOR.temp, summary.LLRpos.temp, summary.LLRneg.temp)
          rownames(summary.summarized.statistics) = c("mean(LRpos)","mean(LRneg)","mean(DOR)", "mean(RD)", "mean(LDOR)", "mean(LLRpos)", "mean(LLRneg)")
          est$summary.summarized.statistics = summary.summarized.statistics
        }else{ # no covariates, but modality
          mod.level = length(unique(outdata$originaldata[,outdata$modality.setting]))
          
          mean.g.fitted1.samples = lapply(1:mod.level, function(i) fixed.samples[ind.fitted1[i],])
          mean.g.fitted2.samples = lapply(1:mod.level, function(i) fixed.samples[ind.fitted2[i],])
          
          mean.fitted1.samples = lapply(1:mod.level, function(i)  model$inv.g(mean.g.fitted1.samples[[i]]))
          mean.fitted2.samples = lapply(1:mod.level, function(i)  model$inv.g(mean.g.fitted2.samples[[i]]))
          
          if(model$model.type==1){
            mean.LRpos.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]]/(1-mean.fitted2.samples[[i]]))
            mean.LRneg.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]])/mean.fitted2.samples[[i]])
            mean.DOR.samples = lapply(1:mod.level, function(i) mean.LRpos.samples[[i]]/mean.LRneg.samples[[i]])
            mean.RD.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]] - (1-mean.fitted2.samples[[i]]))
            mean.LLRpos.samples = lapply(1:mod.level, function(i) log(mean.LRpos.samples[[i]]))
            mean.LLRneg.samples = lapply(1:mod.level, function(i) log(mean.LRneg.samples[[i]]))
            mean.LDOR.samples = lapply(1:mod.level, function(i) log(mean.DOR.samples[[i]]))
          }
          if(model$model.type==2){
            mean.LRpos.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]]/mean.fitted2.samples[[i]])
            mean.LRneg.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]])/(1-mean.fitted2.samples[[i]]))
            mean.DOR.samples = lapply(1:mod.level, function(i) mean.LRpos.samples[[i]]/mean.LRneg.samples[[i]])
            mean.RD.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]] - mean.fitted2.samples[[i]])
            mean.LLRpos.samples = lapply(1:mod.level, function(i) log(mean.LRpos.samples[[i]]))
            mean.LLRneg.samples = lapply(1:mod.level, function(i) log(mean.LRneg.samples[[i]]))
            mean.LDOR.samples = lapply(1:mod.level, function(i) log(mean.DOR.samples[[i]]))
          }
          if(model$model.type==3){
            mean.LRpos.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]])/(1-mean.fitted2.samples[[i]]))
            mean.LRneg.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]]/mean.fitted2.samples[[i]])
            mean.DOR.samples = lapply(1:mod.level, function(i) mean.LRpos.samples[[i]]/mean.LRneg.samples[[i]])
            mean.RD.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]]) - (1-mean.fitted2.samples[[i]]))
            mean.LLRpos.samples = lapply(1:mod.level, function(i) log(mean.LRpos.samples[[i]]))
            mean.LLRneg.samples = lapply(1:mod.level, function(i) log(mean.LRneg.samples[[i]]))
            mean.LDOR.samples = lapply(1:mod.level, function(i) log(mean.DOR.samples[[i]]))
          }
          if(model$model.type==4){
            mean.LRpos.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]])/mean.fitted2.samples[[i]])
            mean.LRneg.samples = lapply(1:mod.level, function(i) mean.fitted1.samples[[i]]/(1-mean.fitted2.samples[[i]]))
            mean.DOR.samples = lapply(1:mod.level, function(i) mean.LRpos.samples[[i]]/mean.LRneg.samples[[i]])
            mean.RD.samples = lapply(1:mod.level, function(i) (1-mean.fitted1.samples[[i]]) - mean.fitted2.samples[[i]])
            mean.LLRpos.samples = lapply(1:mod.level, function(i) log(mean.LRpos.samples[[i]]))
            mean.LLRneg.samples = lapply(1:mod.level, function(i) log(mean.LRneg.samples[[i]]))
            mean.LDOR.samples = lapply(1:mod.level, function(i) log(mean.DOR.samples[[i]]))
          }
          
          summary.LRpos.temp = lapply(1:mod.level, function(i) .summary.samples(mean.LRpos.samples[[i]], level=quantiles))
          summary.LRneg.temp = lapply(1:mod.level, function(i) .summary.samples(mean.LRneg.samples[[i]], level=quantiles))
          summary.DOR.temp = lapply(1:mod.level, function(i) .summary.samples(mean.DOR.samples[[i]], level=quantiles))
          summary.RD.temp = lapply(1:mod.level, function(i) .summary.samples(mean.RD.samples[[i]], level=quantiles))
          summary.LDOR.temp = lapply(1:mod.level, function(i) .summary.samples(mean.LDOR.samples[[i]], level=quantiles))
          summary.LLRpos.temp = lapply(1:mod.level, function(i) .summary.samples(mean.LLRpos.samples[[i]], level=quantiles))
          summary.LLRneg.temp = lapply(1:mod.level, function(i) .summary.samples(mean.LLRneg.samples[[i]], level=quantiles))
          
          summary.summarized.statistics = lapply(1:mod.level, function(i) rbind(summary.LRpos.temp[[i]], 
                                                                                summary.LRneg.temp[[i]], 
                                                                                summary.DOR.temp[[i]], 
                                                                                summary.RD.temp[[i]], 
                                                                                summary.LDOR.temp[[i]], 
                                                                                summary.LLRpos.temp[[i]], 
                                                                                summary.LLRneg.temp[[i]]))
          for(j in 1:mod.level){
            rownames(summary.summarized.statistics[[j]]) = c("mean(LRpos)","mean(LRneg)","mean(DOR)", "mean(RD)", "mean(LDOR)", "mean(LLRpos)", "mean(LLRneg)")
          }
          names(summary.summarized.statistics) = unique(outdata$originaldata[,outdata$modality.setting])
          est$summary.summarized.statistics = summary.summarized.statistics
        }
      }
      
      ############# fitted samples to calculate LRpos, LRneg, DOR for each study
      fitted1.samples = model$inv.g(predictors.samples[fitted1.ind,])
      fitted2.samples = model$inv.g(predictors.samples[fitted2.ind,])
      
      if(model$model.type==1){
        LRpos.samples = fitted1.samples/(1-fitted2.samples)
        LRneg.samples = (1-fitted1.samples)/fitted2.samples
        DOR.samples = LRpos.samples/LRneg.samples
        RD.samples = fitted1.samples - (1-fitted2.samples)
        LLRpos.samples = log(LRpos.samples)
        LLRneg.samples = log(LRneg.samples)
        LDOR.samples = log(DOR.samples)
      }
      if(model$model.type==2){
        LRpos.samples = fitted1.samples/fitted2.samples
        LRneg.samples = (1-fitted1.samples)/(1-fitted2.samples)
        DOR.samples = LRpos.samples/LRneg.samples
        RD.samples = fitted1.samples - fitted2.samples
        LLRpos.samples = log(LRpos.samples)
        LLRneg.samples = log(LRneg.samples)
        LDOR.samples = log(DOR.samples)
      }
      if(model$model.type==3){
        LRpos.samples = (1-fitted1.samples)/(1-fitted2.samples)
        LRneg.samples = fitted1.samples/fitted2.samples
        DOR.samples = LRpos.samples/LRneg.samples
        RD.samples = (1-fitted1.samples) - (1-fitted2.samples)
        LLRpos.samples = log(LRpos.samples)
        LLRneg.samples = log(LRneg.samples)
        LDOR.samples = log(DOR.samples)
      }
      if(model$model.type==4){
        LRpos.samples = (1-fitted1.samples)/fitted2.samples
        LRneg.samples = fitted1.samples/(1-fitted2.samples)
        DOR.samples = LRpos.samples/LRneg.samples
        RD.samples = (1-fitted1.samples) - fitted2.samples
        LLRpos.samples = log(LRpos.samples)
        LLRneg.samples = log(LRneg.samples)
        LDOR.samples = log(DOR.samples)
      }
      
      summary.fitted.LRpos.temp = t(apply(LRpos.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.LRpos.temp) = studynames
      summary.fitted.LRneg.temp = t(apply(LRneg.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.LRneg.temp) = studynames
      summary.fitted.DOR.temp = t(apply(DOR.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.DOR.temp) = studynames
      summary.fitted.RD.temp = t(apply(RD.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.RD.temp) = studynames
      summary.fitted.LDOR.temp = t(apply(LDOR.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.LDOR.temp) = studynames
      summary.fitted.LLRpos.temp = t(apply(LLRpos.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.LLRpos.temp) = studynames
      summary.fitted.LLRneg.temp = t(apply(LLRneg.samples, 1, function(x) .summary.samples(x, level=quantiles)))
      rownames(summary.fitted.LLRneg.temp) = studynames
      
      
      est$summary.fitted.LRpos = summary.fitted.LRpos.temp
      est$summary.fitted.LRneg = summary.fitted.LRneg.temp
      est$summary.fitted.DOR = summary.fitted.DOR.temp
      est$summary.fitted.RD = summary.fitted.RD.temp
      est$summary.fitted.LDOR = summary.fitted.LDOR.temp
      est$summary.fitted.LLRpos = summary.fitted.LLRpos.temp
      est$summary.fitted.LLRneg = summary.fitted.LLRneg.temp
      
      if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
        if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
          if(model$model.type==1){
            est$Se.samples = fitted1.samples
            est$Sp.samples = fitted2.samples
            est[["mean(Se).samples"]] = mean.fitted1.samples
            est[["mean(Sp).samples"]] = mean.fitted2.samples
          }
          if(model$model.type==2){
            est$Se.samples = fitted1.samples
            est$Sp.samples = 1-fitted2.samples
            est[["mean(Se).samples"]] = mean.fitted1.samples
            est[["mean(Sp).samples"]] = 1-mean.fitted2.samples
          }
          if(model$model.type==3){
            est$Se.samples = 1-fitted1.samples
            est$Sp.samples = fitted2.samples
            est[["mean(Se).samples"]] = 1-mean.fitted1.samples
            est[["mean(Sp).samples"]] = mean.fitted2.samples
          }
          if(model$model.type==4){
            est$Se.samples = 1-fitted1.samples
            est$Sp.samples = 1-fitted2.samples
            est[["mean(Se).samples"]] = 1-mean.fitted1.samples
            est[["mean(Sp).samples"]] = 1-mean.fitted2.samples
          }
        }else{
          if(model$model.type==1){
            est$Se.samples = fitted1.samples
            est$Sp.samples = fitted2.samples
            est[["mean(Se).samples"]] = mean.fitted1.samples
            est[["mean(Sp).samples"]] = mean.fitted2.samples
          }
          if(model$model.type==2){
            est$Se.samples = fitted1.samples
            est$Sp.samples = 1-fitted2.samples
            est[["mean(Se).samples"]] = mean.fitted1.samples
            est[["mean(Sp).samples"]] = lapply(1:mod.level, function(i) 1-mean.fitted2.samples[[i]])
          }
          if(model$model.type==3){
            est$Se.samples = 1-fitted1.samples
            est$Sp.samples = fitted2.samples
            est[["mean(Se).samples"]] = lapply(1:mod.level, function(i) 1-mean.fitted1.samples[[i]])
            est[["mean(Sp).samples"]] = mean.fitted2.samples
          }
          if(model$model.type==4){
            est$Se.samples = 1-fitted1.samples
            est$Sp.samples = 1-fitted2.samples
            est[["mean(Se).samples"]] = lapply(1:mod.level, function(i) 1-mean.fitted1.samples[[i]])
            est[["mean(Sp).samples"]] = lapply(1:mod.level, function(i) 1-mean.fitted2.samples[[i]])
          }
        }
      }
    }
    
    variables.names = tolower(colnames(est$data))
    ############# scores
    est$waic = model[["waic"]]
    est$mlik = model[["mlik"]]
    est$cpo = model[["cpo"]]
    est$dic = model[["dic"]]
    
    ############# save inla result
    est$inla.result = model
    ############# save general setting
    est$misc$var.prior.list = outpriors$original.setting$var1
    est$misc$var2.prior.list = outpriors$original.setting$var2
    est$misc$cor.prior.list = outpriors$original.setting$cor
    est$misc$wishart.par = outpriors$original.setting$wishart.par
    est$misc$wishart.flag = model$wishart.flag
    
    if(is.null(outdata$covariates.setting) || outdata$covariates.setting==FALSE){
      est$misc$covariates.flag=FALSE
    } else{
      est$misc$covariates.flag=TRUE
      est$misc$covariates.name = outdata$covariates.setting
      est$misc$covariates.place.in.data = which(variables.names==outdata$covariates.setting)
    }
    if(is.null(outdata$modality.setting) || outdata$modality.setting==FALSE){
      est$misc$modality.flag=FALSE
    }else{
      est$misc$modality.flag=TRUE
      est$misc$modality.name = outdata$modality.setting
      est$misc$modality.place.in.data = which(variables.names==outdata$modality.setting)
      est$misc$modality.level = length(unique(outdata$originaldata[,outdata$modality.setting]))
    }
    est$misc$link = model$link
    est$misc$model.type = model$model.type
    est$misc$quantiles = model$quantiles
    est$misc$verbose = model$verbose
    est$misc$g = model$g
    est$misc$inv.g = model$inv.g
    ################### make new class
    class(est) = 'meta4diag'
    return(est)
  }else{
    stop("INLA need to be installed and loaded!\n
         Please use the following commants to install and load INLA,\n
         install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/testing\")
         library(INLA) \n")
  }
}