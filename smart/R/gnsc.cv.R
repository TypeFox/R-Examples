gnsc.cv <-
function(fit, x, y=NULL, z=NULL, nfold=NULL, folds=NULL, verbose = T){
     if(is.null(y)) stop("y must not be empty")
     n=length(y)
     nlambda = fit$nlambda
     lambda = fit$lambda
     yhat = matrix(0,n,nlambda)
     
     if(length(y)!=ncol(x)) stop("y must have the same length as ncol(x)")
     if(is.null(nfold)) nfold = ceiling(length(y)/min(table(y)))
     if(is.null(folds)) folds = sample(1:nfold,n,replace=T)
     if(is.null(z)) z = 1:dim(x)[1]

     if(max(folds)!=nfold) stop("folds must match nfold")
     x=x[order(z),];z=z[order(z)]        
     
     errors = rep(0,nlambda) 
     if(verbose) cat("Conducting GNSC crossvalidation: \n")     
     for(ci in 1:nfold){
         xx = x[,folds!=ci]; col.struc = y[folds!=ci]
         fit.tmp = gnsc.train(xx, col.struc=col.struc, row.struc=fit$row.struc, standardize = fit$standardize, 
         lambda = fit$lambda, verbose = F)  
         for(sim in 1:nlambda){      
            out.tmp = gnsc.predict(x[,folds==ci], x.class=y[folds==ci], fit.tmp$all.mean, fit.tmp$tilde.mu[,,sim], 
                 fit.tmp$row.struc, fit.tmp$col.struc, 
                 fit.tmp$icov, fit.tmp$path.se.id)
            errors[sim] = errors[sim]+out.tmp$error
         }
         if(verbose) cat("Conducting GNSC crossvalidation:", floor(100*ci/nfold),"%","\r")   
     }    

     out = list()
     errors = errors/nfold
     out$lambda.min = which.min(errors) 
     out$nonzero = fit$nonzero
     out$errors = errors
     out$lambda = lambda
     out$nlambda = nlambda 
     out$Thresh.mat = fit$Thresh.mat
     class(out) = "gnsccv"
     if(verbose) cat("done\n")
     
     rm(xx,fit,x,y,out.tmp,fit.tmp)
     return(out)

}
