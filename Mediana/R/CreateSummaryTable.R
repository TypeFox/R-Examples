######################################################################################################################

# Function: CreateSummaryTable.
# Argument: Results returned by the CSE function.
# Description: This function is used to create a summary table with all results

CreateSummaryTable = function(evaluation.result){
  nscenario = length(evaluation.result)
  table.evaluation.result = list()
  for (i in 1:nscenario){
    ncriterion = length(evaluation.result[[i]]$criterion)
    table.list = list()
    for (j in 1:ncriterion){
      scenario = i
      id = evaluation.result[[i]]$criterion[[j]]$id
      res = format(round(evaluation.result[[i]]$criterion[[j]]$result, digits = 4), digits = 4, nsmall = 4)
      rownames(res) = colnames(res) = NULL
      test = rownames(evaluation.result[[i]]$criterion[[j]]$result)
      table.list[[j]] = data.frame(scenario = scenario, id=rep(id,nrow(res)), test=test, result=res)
    }
    table.evaluation.result[[i]] = do.call(rbind, table.list)
  }
  return(do.call(rbind, table.evaluation.result))
}
# End of CreateSummaryTable