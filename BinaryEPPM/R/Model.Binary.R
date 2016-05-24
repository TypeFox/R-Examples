Model.Binary <-
function(parameter,model.type,model,link,ntrials,
                   covariates.matrix.p,covariates.matrix.scalef,
                   offset.p,offset.scalef) {
   if (model.type=="p only") { 
      if ((model=="binomial") | (model=="generalized binomial")) { 
            output <- Model.GB(parameter,model,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               offset.p=offset.p) 
                                         } # end of Model.GB
      if ((model=="beta binomial") | (model=="correlated binomial")) { 
            output <- Model.BCBinProb(parameter,model.type,model,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               offset.p=offset.p) 
                                         } # end of model.BCBinProb
                                } # end of p only models
   if (model.type=="p and scale-factor") { 
      if (model=="generalized binomial") { 
            output <- Model.JMVGB(parameter,model,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               covariates.matrix.scalef=covariates.matrix.scalef,
                               offset.p=offset.p,offset.scalef=offset.scalef) 
                                         } # end of model.JMVGB
      if ((model=="beta binomial") | (model=="correlated binomial")) { 
            output <- Model.BCBinProb(parameter,model.type,model,link,ntrials,
                               covariates.matrix.p=covariates.matrix.p,
                               covariates.matrix.scalef=covariates.matrix.scalef,
                               offset.p=offset.p,offset.scalef=offset.scalef) 
                                         } # end of model.BCBinProb
                                } # end of p and scale-factor models
   return(output) }
