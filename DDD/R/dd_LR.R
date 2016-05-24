dd_LR = function(             
   brts,
   initparsoptDD,
   initparsoptCR,
   missnumspec,   
   outputfilename = NULL,
   seed = 42,
   endmc = 1000,
   alpha = 0.05,
   plotit = TRUE,
   res = 10 * (1 + length(brts) + missnumspec),
   ddmodel = 1,
   cond = 1,
   btorph = 1,
   soc = 2,
   tol = c(1E-3,1E-4,1E-6),
   maxiter = 2000,
   changeloglikifnoconv = FALSE,
   optimmethod = 'subplex',
   methode = 'analytical'   
   )
{
  if(!is.null(seed))
  {
     set.seed(roundn(seed))
  }
  age = max(brts)
  cat("Estimating parameters under the constant-rate model ...\n")  
  outCRO = dd_ML(brts = brts,initparsopt = initparsoptCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
  cat("Estimating parameters under the diversity-dependent model ...\n")  
  outDDO = dd_ML(brts = brts,initparsopt = initparsoptDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
  LRO = outDDO$loglik - outCRO$loglik
  out = cbind(NA,NA,outCRO,outDDO,NA,NA,NA,NA,NA,LRO)
  out = out[,-c(5,7,13)]
  newnames = c("model","mc","lambda_CR","mu_CR","LL_CR","conv_CR","lambda_DD1","mu_DD1","K_DD1","LL_DD1","conv_DD1","lambda_DD2","mu_DD2","K_DD2","LL_DD2","conv_DD2","LR")
  names(out) = newnames
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,file = outputfilename)
  }
  parsCR = as.numeric(outCRO[1:2])
  parsDD = as.numeric(outDDO[1:3])
  treeCR = list()
  treeDD = list()
  cat('Simulating trees under CR and DD models ...\n')
  for(mc in 1:endmc)
  {
     treeCR[[mc]] = dd_sim(pars = c(parsCR,Inf),age = age,ddmodel = 1)
     treeDD[[mc]] = dd_sim(pars = parsDD,age = age,ddmodel = 1)
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,treeCR,treeDD,file = outputfilename)
  }
  cat('Performing bootstrap for determining critical LR ...\n')  
  for(mc in 1:endmc)
  {
     cat('Analyzing simulation:',mc,'\n')
     brtsCR = branching.times(treeCR[[mc]][[1]])
     outCR = dd_ML(brtsCR,initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     outDD1 = dd_ML(brtsCR,initparsopt = parsDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     outDD2 = dd_ML(brtsCR,initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     if(outDD1$conv == -1 & outDD2$conv == -1)
     {
        maxLLDD = outCR$loglik
     } else if(outDD1$conv == -1 & outDD2$conv != -1)
     {
        maxLLDD = outDD2$loglik
     } else if(outDD1$conv != -1 & outDD2$conv == -1)
     {
        maxLLDD = outDD1$loglik
     } else {
        maxLLDD = max(outDD1$loglik,outDD2$loglik)
     }   
     LR = pmax(0,maxLLDD - outCR$loglik)
     outff = cbind(1,mc,outCR,outDD1,outDD2,LR)
     outff = outff[,-c(5,7,13,19)]     
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,treeCR,treeDD,file = outputfilename)
     }
  }
  opt = rep(0,endmc)
  cat('Performing bootstrap for determine power ...\n')
  for(mc in 1:endmc)
  {
     print(mc)
     brtsDD = branching.times(treeDD[[mc]][[1]])
     outCR = dd_ML(brtsDD,initparsopt = parsCR,idparsopt = 1:2, idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     outDD1 = dd_ML(brtsDD,initparsopt = parsDD,idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     outDD2 = dd_ML(brtsDD,initparsopt = c(parsCR + 0.05,length(brts) + 1000), idparsopt = 1:3,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv, optimmethod = optimmethod)
     if(outDD1$conv == -1 & outDD2$conv == -1)
     {
        maxLLDD = outCR$loglik
        opt[mc] = 1
     } else if(outDD1$conv != -1 & outDD2$conv == -1)
     {
        maxLLDD = outDD1$loglik
        opt[mc] = 2
     } else if(outDD1$conv == -1 & outDD2$conv != -1)
     {
        maxLLDD = outDD2$loglik
        opt[mc] = 3
     } else {
        maxLLDD = max(outDD1$loglik,outDD2$loglik)
        opt[mc] = 1 + min(which(c(outDD1$loglik,outDD2$loglik) == maxLLDD))
     }
     LR = pmax(0,maxLLDD - outCR$loglik)
     outff = cbind(1,mc,outCR,outDD1,outDD2,LR)
     outff = outff[,-c(5,7,13,19)]
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,opt,treeCR,treeDD,file = outputfilename)
     }
  }
  inverse_quantile = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     if(length(pup) > 0)
     {
        if(length(pup) < length(samplessort))
        {
           pup = min(pup)
           invquant = (pup + (x - samplessort[pup])/(samplessort[pup - 1] - samplessort[pup]))/length(samples)
        } else {
           invquant = 0
        }
     } else {
        invquant = 1
     }
     return(invquant)
  }
  funpvalue = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     pvalue = (length(pup) + 1)/ (length(samples) + 1)
     return(pvalue)
  }
  funpoweroftest = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     poweroftest = length(pup)/(length(samples) + 1)
     return(poweroftest)     
  }
  #pvalue = 1 - inverse_quantile(out$LR[2:(endmc + 1)],out$LR[1])
  pvalue = funpvalue(out$LR[2:(endmc + 1)],out$LR[1])
  LRalpha = as.numeric(quantile(out$LR[2:(endmc + 1)],1 - alpha,type = 4))
  #poweroftest = 1 - inverse_quantile(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  poweroftest = funpoweroftest(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  if(plotit == TRUE)
  {
      try(dev.off())
      try(dev.off())
      pdffilename = paste(getwd(),'/LR.pdf',sep = '')
      pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
      al = 0.03
      alw = 2
      alw2 = 1.7
      aa = 45
      par(mfrow = c(2,2),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      hist(out$LR[2:(1 + endmc)],main = 'Distribution of LLR under CR',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(out$LR[1:(endmc + 1)])))
      arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')
      box()
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", col = "white", col.lab = 'white', col.axis = 'white')
      hist(out$LR[(endmc + 2):(1 + 2 * endmc)], main = 'Distribution of LLR under DD',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30)
      box()
      arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')

      par(mfrow = c(2,3),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      lambda = out$lambda_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$lambda_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$lambda_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      mu = out$mu_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$mu_DD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$mu_DD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      K = 1E+120 * (opt == 1) + pmin(1E+120,out$K_DD1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K_DD2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
      hist(lambda,main = NULL, xlab = expression(lambda), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(lambda)))
      arrows(out$lambda_DD1[1],-1E+120, x1 = out$lambda_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(mu,main = NULL, xlab = expression(mu), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(mu)))
      arrows(out$mu_DD1[1],-1E+120, x1 = out$mu_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(K,main = NULL, xlab = 'K', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(min(K),max(K)))
      arrows(out$K_DD1[1],-1E+120, x1 = out$K_DD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      try(dev.off())
      try(dev.off())
      os = .Platform$OS.type
      if(os == "windows")
      {
          shell.exec(pdffilename)
      }
      if(os == "unix")
      {
          system(paste("open",pdffilename,sep = " "))
      }
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,opt,treeCR,treeDD,pvalue,LRalpha,poweroftest,file = outputfilename)
  }
  return(list(treeCR = treeCR,treeDD = treeDD,out = out,pvalue = pvalue,LRalpha = LRalpha,poweroftest = poweroftest))
}