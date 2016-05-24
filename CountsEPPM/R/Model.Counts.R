Model.Counts <-
function(parameter,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,
                   scale.factor.model,fixed.b,vnmax) {
   if (model.type=="mean only") { 
                   output <- Model.Faddy(parameter,model,covariates.matrix.mean,
                                       offset.mean,fixed.b,vnmax)
                                } # end of over-dispersed models
   if (model.type=="mean and variance") { 
          if ((model=="general") | (model=="general fixed b")) {
                   output <- Model.FaddyJMV.general(parameter,covariates.matrix.mean,
                                         covariates.matrix.variance,offset.mean,
                                         offset.variance,scale.factor.model,fixed.b,vnmax) }
          if (model=="limiting") { 
                   output <- Model.FaddyJMV.limiting(parameter,covariates.matrix.mean,
                                         covariates.matrix.variance,offset.mean,
                                         offset.variance,scale.factor.model,vnmax) }
                                } # end of under-dispersed model
   return(output)                          }
