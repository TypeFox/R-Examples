
BM_bernoulli_multiplex <-
  setRefClass("BM_bernoulli_multiplex",
              contains = "multivariate_model",
              methods = list(
                initialize = function(membership_type,adj,...)
                {
                  .self$initFields(membership_name = membership_type,
                                   model_name = "bernoulli_multiplex",
                                   adj = adj,
                                   ...)
                  .self$postinit()
                },
                plot_parameters = function(Q)
                {
                },
                prediction = function(Q)
                {
                  noms<-names(model_parameters[[Q]])
                  if(membership_name=='LBM')
                  {
                    pred <- list()
                    for(k in 1:length(adj))
                    {
                      pred[[k]]<-matrix(0,nrow=nrow(adj[[1]]),ncol=ncol(adj[[1]]))
                      
                      for (i in 1:length(noms))
                      {
                          if (substr(noms[i],k,k)=='1')
                          {  
                             pred[[k]]<-pred[[k]]+ memberships[[Q]]$Z1 %*% model_parameters[[Q]]$pi[[i]] %*% t(memberships[[Q]]$Z2)
                          }  
                      }
                    }
                    return(pred)
                  }
                  else
                  {
                    pred <- list()
                    for(k in 1:length(adj))
                    {
                      pred[[k]]<-matrix(0,nrow=nrow(adj[[1]]),ncol=nrow(adj[[1]]))
                      for (i in 1:length(noms))
                      {
                        if (substr(noms[i],k,k)=='1')
                        {  
                          pred[[k]]<-pred[[k]]+ memberships[[Q]]$Z %*% model_parameters[[Q]]$pi[[i]] %*% t(memberships[[Q]]$Z)
                        }  
                      }
                    }
                    return(pred)
                  }
                },
                residual = function(Q)
                {
                  pred <- .self$prediction(Q)
                  ret <- list()
                  for(k in 1:length(adj))
                  {
                    ret[[k]] <- adj[[k]]-pred[[k]]
                  }
                  return(ret)
                } 
              )
  )
