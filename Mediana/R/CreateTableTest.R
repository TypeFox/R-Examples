############################################################################################################################

# Function: CreateTableTest.
# Argument: analysis.strucure and label (optional).
# Description: Generate a summary table of test for the report.

CreateTableTest = function(analysis.structure, label = NULL) {

  # Number of test
  n.test = length(analysis.structure$test)

  test.table = matrix(nrow = n.test, ncol = 4)
  nsample = rep(0,n.test)
  for (i in 1:n.test) {
    test.table[i, 1] = analysis.structure$test[[i]]$id
    test.desc = do.call(analysis.structure$test[[i]]$method,list(c(),list("Description",analysis.structure$test[[i]]$par)))
    test.table[i, 2] = test.desc[[1]]
    if (length(test.desc)>1) {
      test.table[i, 3] = paste0(test.desc[[2]],analysis.structure$test[[i]]$par)
    } else {
      test.table[i, 3] = analysis.structure$test[[i]]$par
    }
    nsample[i]=length(analysis.structure$test[[i]]$samples)
    npersample=rep(0,nsample[i])
    sample.id=rep("",nsample[i])
    text=""
    for (j in 1:nsample[i]) {
      npersample[j]=length(analysis.structure$test[[i]]$samples[[j]])
      for (k in 1:npersample[j]) {
        sample.id[j]=paste0(sample.id[j],", ", analysis.structure$test[[i]]$samples[[j]][[k]])
      }
      sample.id[j]=paste0("{",sub(", ","",sample.id[j]),"}")
      text=paste0(text,", ",sample.id[j])
    }
    test.table[i, 4] = sub(", ","",text)
  }

  test.table = as.data.frame(test.table)
  colnames(test.table) = c("Test ID", "Test type", "Test parameters", "Samples")
  return(test.table)
}
# End of CreateTableTest