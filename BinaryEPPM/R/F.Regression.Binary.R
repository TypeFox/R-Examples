F.Regression.Binary <-
function(parameter,model.type,model,link,ntrials,nsuccess,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef) {
   loglikelihood <- - LL.Regression.Binary(parameter,model.type,model,link,ntrials,nsuccess,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef) 
# numerical calculation of the gradient for nlm using latest numerical derivative routine
# and adding it as an attribute to loglikelihood. Because nlm is function
# minimization the negatives of the loglikelihood and its gradients are required
  attr(loglikelihood,"gradient") <- - grad(LL.Regression.Binary,
                x=parameter,method="Richardson",
                model.type=model.type,model=model,link=link,
                ntrials=ntrials,nsuccess=nsuccess,
                covariates.matrix.p=covariates.matrix.p,
                covariates.matrix.scalef=covariates.matrix.scalef,
                offset.p=offset.p,offset.scalef=offset.scalef)
   return(loglikelihood) }
