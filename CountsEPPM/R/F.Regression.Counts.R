F.Regression.Counts <-
function(parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b) {
   loglikelihood <- - LL.Regression.Counts(parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b) 

# numerical calculation of the gradient for nlm using latest numerical derivative routine
# and adding it as an attribute to loglikelihood. Because nlm is function
# minimization the negatives of the loglikelihood and its gradients are required

  attr(loglikelihood,"gradient") <- - grad(LL.Regression.Counts,
                x=parameter,method="Richardson",
                model.type=model.type,model=model,list.counts=list.counts,
                covariates.matrix.mean=covariates.matrix.mean,
                covariates.matrix.variance=covariates.matrix.variance,
                offset.mean=offset.mean,offset.variance=offset.variance,
                ltvalue=ltvalue,utvalue=utvalue,
                scale.factor.model=scale.factor.model,fixed.b=fixed.b)

   return(loglikelihood)                          }
