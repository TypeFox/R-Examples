mySummary <- function(results, study.name, drug.name){

  # pull out summary tables and models
  summary_results = NULL
  model_results = list()
  for (result_index in 1:length(results)){
    summary_results = rbind(summary_results, results[[result_index]][[1]])
    model_results = c(model_results, results[[result_index]][[2]]) 
  }  
    
  # rbind all dataframes
  #results <- ldply(results, data.frame)

  # eliminate duplicates for etas-slope which repeat accross Eta-EC50/EMAX
  summary_results <- unique(summary_results)

  # in the summary, give all models tested order from the best (lower AIC) to worse (higher AIC)
  out_summary <- summary_results[with(summary_results, order(AIC)),]

  # order models in the same way
  out_models = list()
  for (model in out_summary$MODEL){
    out_models[[model]] = model_results[[model]]
  }
  
  # remove first part of ERROR.MESSAGE so that function names don't show
  out_summary$ERROR.MESSAGE <- gsub(".* : \\n  ", "", out_summary$ERROR.MESSAGE)
  out_summary$ERROR.MESSAGE <- gsub("\\n$", "", out_summary$ERROR.MESSAGE)
  
  # add DRUG and STUDY columns to the summary table
  out_summary$DRUG <- drug.name
  out_summary$STUDY <- study.name
  
  # rearrange columns of the summary table
  out_summary <- out_summary[, c(ncol(out_summary), ncol(out_summary)-1, 1:(ncol(out_summary)-2))]
  return(list(summary = out_summary, models = out_models))
}
