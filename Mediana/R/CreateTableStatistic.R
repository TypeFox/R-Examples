############################################################################################################################

# Function: CreateTableStatistic.
# Argument: analysis.strucure and label (optional).
# Description: Generate a summary table of statistic for the report.

CreateTableStatistic = function(analysis.structure, label = NULL) {

  # Number of statistic
  n.statistic = length(analysis.structure$statistic)

  statistic.table = matrix(nrow = n.statistic, ncol = 4)
  nsample = rep(0,n.statistic)
  for (i in 1:n.statistic) {
    statistic.table[i, 1] = analysis.structure$statistic[[i]]$id
    statistic.desc = do.call(analysis.structure$statistic[[i]]$method,list(c(),list("Description",analysis.structure$statistic[[i]]$par)))
    statistic.table[i, 2] = statistic.desc[[1]]
    if (length(statistic.desc)>1) {
      statistic.table[i, 3] = paste0(statistic.desc[[2]],analysis.structure$statistic[[i]]$par)
    } else {
      statistic.table[i, 3] = analysis.structure$statistic[[i]]$par
    }
    nsample[i]=length(analysis.structure$statistic[[i]]$samples)
    npersample=rep(0,nsample[i])
    sample.id=rep("",nsample[i])
    text=""
    for (j in 1:nsample[i]) {
      npersample[j]=length(analysis.structure$statistic[[i]]$samples[[j]])
      for (k in 1:npersample[j]) {
        sample.id[j]=paste0(sample.id[j],", ", analysis.structure$statistic[[i]]$samples[[j]][[k]])
      }
      sample.id[j]=paste0("{",sub(", ","",sample.id[j]),"}")
      text=paste0(text,", ",sample.id[j])
    }
    statistic.table[i, 4] = sub(", ","",text)
  }

  statistic.table = as.data.frame(statistic.table)
  colnames(statistic.table) = c("Statistic ID", "Statistic type", "Statistic parameters", "Samples")
  return(statistic.table)
}
# End of CreateTableStatistic