#' Summarises the bootstrap results.
#'
#' Creates summary statistics for each species. These consist of dataframes 
#' relating to summaries, abundance (N) and density (D) for both individuals 
#' and clusters. In addition, summary statistics for expected cluster size 
#' (Expected.S) are also calculated.
#'  
#' @param bootstrap.results list of arrays containg results from the repeated
#'   analyses.
#' @param model.index named character vector which acts as a look up table for
#'   duplicate detection function models.
#' @param clusters boolean whether observations are clusters of individuals
#' @param bootstrap.ddf.statistics array storing parameter estimates from ddf models
#' @param quantile.type numeric value describing which quantile algorithm to use 
#' @param analysis.options list describing the type of analysis carried out
#' @return ma object a list of summary statistics for each species
#' @note Internal functions not intended to be called by user.
#' @author Laura Marshall
#' @importFrom stats sd
#' @importFrom stats quantile
#'
process.bootstrap.results <- function(bootstrap.results, model.index, clusters, bootstrap.ddf.statistics, quantile.type, analysis.options = list(bootstrap, n, covariate.uncertainty, clusters, double.observer, unidentified.species, species.code.definitions, model.names)){
#process.bootstrap.results function to summarise bootstrap results
#
# Arguments:
#   bootstrap.results - list of arrays containg bootstrap results
#
# Value: 
#   ma object a list of summary statistics for each species
#
# Function Calls: none
#
  #analysis options
  n                        <- analysis.options$n
  bootstrap                <- analysis.options$bootstrap
  covariate.uncertainty    <- analysis.options$covariate.uncertainty 
  clusters                 <- analysis.options$clusters
  double.observer          <- analysis.options$double.observer
  unidentified.species     <- analysis.options$unidentified.species
  species.code.definitions <- analysis.options$species.code.definitions 
  model.names              <- analysis.options$model.names 
  
  #set up data storage for results
  results.summary <- list()
  quantile.type <- as.numeric(quantile.type)

  #add the analysis options to the result summary
  results.summary$analysis.options$bootstrap                <- bootstrap
  results.summary$analysis.options$n                        <- n
  results.summary$analysis.options$covariate.uncertainty    <- covariate.uncertainty
  results.summary$analysis.options$clusters                 <- clusters
  results.summary$analysis.options$double.observer          <- double.observer
  results.summary$analysis.options$unidentified.species     <- unidentified.species
  results.summary$analysis.options$species.code.definitions <- species.code.definitions
  results.summary$analysis.options$model.names              <- model.names
  
  #gather, array names, species.names, strata names, and number of strata       
  species.names <- dimnames(bootstrap.results[[1]])[[4]]
  strata.names <- dimnames(bootstrap.results[[1]])[[1]]
  no.strata  <- length(strata.names)
  all.species.codes <- names(model.index)
  unid.names <- all.species.codes[-which(all.species.codes%in%species.names)]
  
  
  #for each species summarise the bootstrap results
  for(sp in seq(along = species.names)){
  
    #object to store results in
    summary.element <- list()
    
    ### INDIVIDUALS ###
    #...$summary
    summary.element$individuals$summary <- data.frame(Region = strata.names, Area = rep(NA,no.strata), CoveredArea = rep(NA,no.strata), Effort = rep(NA,no.strata), n = rep(NA,no.strata), ER = rep(NA,no.strata), se.ER = rep(NA,no.strata), cv.ER = rep(NA,no.strata))
    for(st in seq(along = strata.names)){
      if(class(bootstrap.results[["individual.summary"]][strata.names[st],,,species.names[sp]]) == "matrix"){
        summary.element$individuals$summary[summary.element$individuals$summary$Region == strata.names[st],2:6] <- apply(bootstrap.results[["individual.summary"]][strata.names[st],,,species.names[sp]], 1, mean, na.rm = TRUE)
      }else if(class(bootstrap.results[["individual.summary"]][strata.names[st],,,species.names[sp]]) == "numeric"){
        summary.element$individuals$summary[summary.element$individuals$summary$Region == strata.names[st],2:6] <- bootstrap.results[["individual.summary"]][strata.names[st],,,species.names[sp]]
      }
      summary.element$individuals$summary[summary.element$individuals$summary$Region == strata.names[st],7] <- sd(bootstrap.results[["individual.summary"]][strata.names[st],5,,species.names[sp]], na.rm = TRUE)
    }
    summary.element$individuals$summary$cv.ER <- summary.element$individuals$summary$se.ER/summary.element$individuals$summary$ER
    
    #...$N 
    summary.element$individuals$N <- data.frame(Label = strata.names, Estimate = rep(NA,no.strata), se = rep(NA,no.strata), cv = rep(NA,no.strata), lcl = rep(NA,no.strata), ucl = rep(NA,no.strata), df = rep(NA,no.strata), pctUnid = rep(NA,no.strata), pctUnid.se = rep(NA,no.strata), pctUnid.cv = rep(NA,no.strata), pctUnid.lcl = rep(NA,no.strata), pctUnid.ucl = rep(NA,no.strata))
    for(st in seq(along = strata.names)){
      if(class(bootstrap.results[["individual.N"]][strata.names[st],,,species.names[sp]]) == "matrix"){
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(2,7,8)] <- apply(bootstrap.results[["individual.N"]][strata.names[st],,,species.names[sp]], 1, mean, na.rm = TRUE)
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(3,9)] <- apply(bootstrap.results[["individual.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, sd, na.rm = TRUE)
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(5,11)] <- apply(bootstrap.results[["individual.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, quantile, na.rm = TRUE, type = quantile.type, probs = 0.025)
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(6,12)] <- apply(bootstrap.results[["individual.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, quantile, na.rm = TRUE, type = quantile.type, probs = 0.975)
      }else if(class(bootstrap.results[["individual.N"]][strata.names[st],,,species.names[sp]]) == "numeric"){
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(2,7,8)] <- bootstrap.results[["individual.N"]][strata.names[st],,,species.names[sp]]
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(3,9)] <- rep(NA,2)
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(5,11)] <- bootstrap.results[["individual.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]]
        summary.element$individuals$N[summary.element$individuals$N$Label == strata.names[st], c(6,12)] <- bootstrap.results[["individual.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]]
      }      
    }
    summary.element$individuals$N$cv <- summary.element$individuals$N$se/summary.element$individuals$N$Estimate 
    summary.element$individuals$N$pctUnid.cv <- summary.element$individuals$N$pctUnid.se/summary.element$individuals$N$pctUnid
    
    #...$D
    summary.element$individuals$D <- data.frame(Label = strata.names, Estimate = rep(NA,no.strata), se = rep(NA,no.strata), cv = rep(NA,no.strata), lcl = rep(NA,no.strata), ucl = rep(NA,no.strata), df = rep(NA,no.strata))
    summary.element$individuals$D$Estimate <- summary.element$individuals$N$Estimate/summary.element$individuals$summary$Area
    summary.element$individuals$D$cv <- summary.element$individuals$N$cv
    summary.element$individuals$D$lcl <- summary.element$individuals$N$lcl/summary.element$individuals$summary$Area
    summary.element$individuals$D$ucl <- summary.element$individuals$N$ucl/summary.element$individuals$summary$Area
    summary.element$individuals$D$df <- summary.element$individuals$N$df
    summary.element$individuals$D$se <- summary.element$individuals$D$Estimate*summary.element$individuals$N$cv
  
    if(clusters){    
      ### CLUSTERS ###
      #...$summary
     summary.element$clusters$summary <- data.frame(Region = strata.names, Area = rep(NA,no.strata), CoveredArea = rep(NA,no.strata), Effort = rep(NA,no.strata), n = rep(NA,no.strata), k = rep(NA,no.strata), ER = rep(NA,no.strata), se.ER = rep(NA,no.strata), cv.ER = rep(NA,no.strata))
      for(st in seq(along = strata.names)){
        if(class(bootstrap.results[["clusters.summary"]][strata.names[st],,,species.names[sp]]) == "matrix"){
          summary.element$clusters$summary[summary.element$clusters$summary$Region == strata.names[st],c(2:7)] <- apply(bootstrap.results[["clusters.summary"]][strata.names[st],,,species.names[sp]], 1, mean, na.rm = TRUE)
        }else if(class(bootstrap.results[["clusters.summary"]][strata.names[st],,,species.names[sp]]) == "numeric"){
          summary.element$clusters$summary[summary.element$clusters$summary$Region == strata.names[st],c(2:7)] <- bootstrap.results[["clusters.summary"]][strata.names[st],,,species.names[sp]]
        }                         
        summary.element$clusters$summary[summary.element$clusters$summary$Region == strata.names[st],8] <- sd(bootstrap.results[["clusters.summary"]][strata.names[st],6,,species.names[sp]], na.rm = TRUE)
      }
      summary.element$clusters$summary$cv.ER <- summary.element$clusters$summary$se.ER/summary.element$clusters$summary$ER
      
      #...$N
      summary.element$clusters$N <- data.frame(Label = strata.names, Estimate = rep(NA,no.strata), se = rep(NA,no.strata), cv = rep(NA,no.strata), lcl = rep(NA,no.strata), ucl = rep(NA,no.strata), df = rep(NA,no.strata), pctUnid = rep(NA,no.strata), pctUnid.se = rep(NA,no.strata), pctUnid.cv = rep(NA,no.strata), pctUnid.lcl = rep(NA,no.strata), pctUnid.ucl = rep(NA,no.strata))
      for(st in seq(along = strata.names)){
        if(class(bootstrap.results[["clusters.N"]][strata.names[st],,,species.names[sp]]) == "matrix"){
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(2,7,8)] <- apply(bootstrap.results[["clusters.N"]][strata.names[st],,,species.names[sp]], 1, mean, na.rm = TRUE)
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(3,9)] <- apply(bootstrap.results[["clusters.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, sd, na.rm = TRUE)
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(5,11)] <- apply(bootstrap.results[["clusters.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, quantile, na.rm = TRUE, type = quantile.type, probs = 0.025)
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(6,12)] <- apply(bootstrap.results[["clusters.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]], 1, quantile, na.rm = TRUE, type = quantile.type, probs = 0.975)
        }else if(class(bootstrap.results[["clusters.N"]][strata.names[st],,,species.names[sp]]) == "numeric"){
           summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(2,7,8)] <- bootstrap.results[["clusters.N"]][strata.names[st],,,species.names[sp]]
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(3,9)] <- rep(NA,2)          
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(5,11)] <- bootstrap.results[["clusters.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]]
          summary.element$clusters$N[summary.element$clusters$N$Label == strata.names[st], c(6,12)] <- bootstrap.results[["clusters.N"]][strata.names[st],c("Estimate","PercentUnidentified"),,species.names[sp]]
        }
      }
      summary.element$clusters$N$cv <- summary.element$clusters$N$se/summary.element$clusters$N$Estimate 
      summary.element$clusters$N$pctUnid.cv <- summary.element$clusters$N$pctUnid.se/summary.element$clusters$N$pctUnid
      
      #...$D
      summary.element$clusters$D <- data.frame(Label = strata.names, Estimate = rep(NA,no.strata), se = rep(NA,no.strata), cv = rep(NA,no.strata), lcl = rep(NA,no.strata), ucl = rep(NA,no.strata), df = rep(NA,no.strata))
      summary.element$clusters$D$Estimate <- summary.element$clusters$N$Estimate/summary.element$clusters$summary$Area
      summary.element$clusters$D$cv <- summary.element$clusters$N$cv
      summary.element$clusters$D$lcl <- summary.element$clusters$N$lcl/summary.element$clusters$summary$Area
      summary.element$clusters$D$ucl <- summary.element$clusters$N$ucl/summary.element$clusters$summary$Area
      summary.element$clusters$D$df <- summary.element$clusters$N$df
      summary.element$clusters$D$se <- summary.element$clusters$D$Estimate*summary.element$clusters$N$cv
      
      ### Expected.S ###
      summary.element$Expected.S <- data.frame(Region = strata.names, Expected.S = rep(NA,no.strata), se.Expected.S = rep(NA,no.strata), pro.Expected.S = rep(NA,no.strata), pro.se.Expected.S = rep(NA,no.strata))
      for(st in seq(along = strata.names)){
        summary.element$Expected.S[summary.element$Expected.S$Region == strata.names[st], 2] <- mean(bootstrap.results[["Expected.S"]][strata.names[st],"Expected.S",,species.names[sp]], na.rm = TRUE)
        summary.element$Expected.S[summary.element$Expected.S$Region == strata.names[st], 3] <- sd(bootstrap.results[["Expected.S"]][strata.names[st],"Expected.S",,species.names[sp]], na.rm = TRUE)
        summary.element$Expected.S[summary.element$Expected.S$Region == strata.names[st], 4] <- mean(bootstrap.results[["Expected.S"]][strata.names[st],"new.Expected.S",,species.names[sp]], na.rm = TRUE)
        summary.element$Expected.S[summary.element$Expected.S$Region == strata.names[st], 5] <- sd(bootstrap.results[["Expected.S"]][strata.names[st],"new.Expected.S",,species.names[sp]], na.rm = TRUE)
      }
    }
    #summarise ddf statistics
    ddf.code <- model.index[[species.names[sp]]]
    summary.element$ddf$convergence <- bootstrap.ddf.statistics[[ddf.code]]$convergence
    model.names <- dimnames(bootstrap.ddf.statistics[[ddf.code]]$convergence)[[2]] 
    for(m in seq(along = model.names)){
      if(!is.null(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)){
        if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[2] > 1){
          summary.element$ddf[[model.names[m]]]$ds.params <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
        }else if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[2] == 1){
          param.name <- dimnames(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[[2]]
          param.data <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
          if(length(param.data) > 0){
            summary.element$ddf[[model.names[m]]]$ds.params <- matrix(param.data, ncol = 1, dimnames = list(1:length(param.data), param.name))
          }else{
            summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
            next
          }
        }
      }
      if(!is.null(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)){
        if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[2] > 1){
          summary.element$ddf[[model.names[m]]]$mr.params <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
        }else if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[2] == 1){
          param.name <- dimnames(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[[2]]
          param.data <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
          if(length(param.data) > 0){
            summary.element$ddf[[model.names[m]]]$mr.params <- matrix(param.data, ncol = 1, dimnames = list(1:length(param.data), param.name))
          }else{
            summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
            next
          }
        }
      }
      criteria <- names(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]])[2]
      summary.element$ddf[[model.names[m]]][[criteria]] <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]][[criteria]][bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1] 
      summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
    }
    #add summary element to results list
    results.summary$species[[species.names[sp]]] <- summary.element     
  }#next species
  #add in summaries for unids
  if(length(unid.names) > 0){
    for(u in seq(along = unid.names)){
      ddf.code <- model.index[[unid.names[u]]]
      summary.element <- list()
      summary.element$ddf$convergence <- bootstrap.ddf.statistics[[ddf.code]]$convergence
      model.names <- dimnames(bootstrap.ddf.statistics[[ddf.code]]$convergence)[[2]] 
      for(m in seq(along = model.names)){
        if(!is.null(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)){
          if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[2] > 1){
            summary.element$ddf[[model.names[m]]]$ds.params <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
          }else if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[2] == 1){
            param.name <- dimnames(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param)[[2]]
            param.data <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$ds.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
            if(length(param.data) > 0){
              summary.element$ddf[[model.names[m]]]$ds.params <- matrix(param.data, ncol = 1, dimnames = list(1:length(param.data), param.name))
            }else{
              summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
              next
            }
          }
        }
        if(!is.null(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)){
          if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[2] > 1){
            summary.element$ddf[[model.names[m]]]$mr.params <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
          }else if(dim(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[2] == 1){
            param.name <- dimnames(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param)[[2]]
            param.data <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$mr.param[bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1,]
            if(length(param.data) > 0){
              summary.element$ddf[[model.names[m]]]$mr.params <- matrix(param.data, ncol = 1, dimnames = list(1:length(param.data), param.name))
            }else{
              summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
              next
            }
          }
        }
        criteria <- names(bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]])[2]
        summary.element$ddf[[model.names[m]]][[criteria]] <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]][[criteria]][bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$selected == 1] 
        summary.element$ddf[[model.names[m]]]$model.description <- bootstrap.ddf.statistics[[ddf.code]][[model.names[m]]]$model.description
      }
      results.summary$unidentified[[unid.names[u]]] <- summary.element
    }
  } 
  return(results.summary) 
}




