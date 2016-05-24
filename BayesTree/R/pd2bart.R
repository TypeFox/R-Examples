pd2bart = function (
   x.train, y.train,
   xind=1:2, levs=NULL, levquants=c(.05,(1:9)/10,.95),
   pl=TRUE, plquants=c(.05,.95), 
   ...
)
{
   n = nrow(x.train)
   nlevels = rep(0,2)
   if(is.null(levs)) {
      levs = list()
      for(i in 1:2) {
         ux = unique(x.train[,xind[i]])
	 if(length(ux) <= length(levquants)) levs[[i]] = sort(ux)
	 else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
      }
   } 
   nlevels = unlist(lapply(levs,length))
   xvals <- as.matrix(expand.grid(levs[[1]],levs[[2]]))
   nxvals <- nrow(xvals)
   if (ncol(x.train)==2){
      cat('special case: only 2 xs\n')
      x.test = xvals
   } else {
      x.test=NULL
      for(v in 1:nxvals) {
         temp = x.train
         temp[,xind[1]] = xvals[v,1]
         temp[,xind[2]] = xvals[v,2]
         x.test = rbind(x.test,temp)
      }
   }
   pdbrt = bart(x.train,y.train,x.test,...)
   if (ncol(x.train)==2) {
      fdr = pdbrt$yhat.test
   } else {
      fdr = NULL 
      for(i in 1:nxvals) {
         cind =  ((i-1)*n+1):(i*n)
         fdr = cbind(fdr,(apply(pdbrt$yhat.test[,cind],1,mean)))
      }
   }
   if(is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
   else xlbs = colnames(x.train)[xind]
   if('sigma' %in% names(pdbrt)) {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,
      yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
   } else {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      y=pdbrt$y)
   }
   class(retval) = 'pd2bart'
   if(pl) plot(retval,plquants=plquants)
   return(retval)
}

