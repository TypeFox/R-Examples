#developed by Max Greene

library(RSNNS)
basePath <- ("./")
data(iris)

set.seed(2)

#shuffle the vector
iris <- iris[sample(nrow(iris)),]

#normalize data
inputs <- normalizeData(iris[,1:4], "norm")

#outputs <- decodeClassLabels(iris[,5])
outputs <- decodeClassLabels(iris[,5], valTrue=0.9, valFalse=0.1)

numHiddenUnits <- c(10,10)

snnsObject <- SnnsRObjectFactory()

snnsObject$setFFLearnFunc("Std_Backpropagation")
snnsObject$setPrunFunc("OptimalBrainSurgeon")

#snnsObject$setLearnFunc("PruningFeedForward")
#snnsObject$setPrunFunc("OptimalBrainSurgeon","Std_Backpropagation", 10.0, 5.0, TRUE, 1000, 100, 1.0, 1e-6, TRUE, TRUE)

snnsObject$setUpdateFunc('Topological_Order')

snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')

snnsObject$createNet(c(ncol(inputs),numHiddenUnits,ncol(outputs)), TRUE)

patset <- snnsObject$createPatSet(inputs, outputs)

snnsObject$setCurrPatSet(patset$set_no)

snnsObject$shufflePatterns(TRUE)

snnsObject$DefTrainSubPat()

snnsObject$initializeNet(-1)

snnsObject$saveNet(paste(basePath,"mlp_irisSnnsR_untrained.net",sep=""),"mlp_irisSnnsR_untrained.net")

parameters <- c(0.2, 0, 0, 0, 0)

no_of_cycles <- 100000

error <- vector()


for(i in 1:no_of_cycles) {
  res <- snnsObject$learnAllPatternsFF(parameters)
  error[i] <- res[[2]]
  if (i %% 100 == 0) print(c((format(i,digits=2)),error[i]))
}


plot(error, type="l")

predictions <- snnsObject$predictCurrPatSet("output", c(0))

confusionMatrix(outputs,predictions)

snnsObject$getFFLearnFunc()
snnsObject$getPrunFunc()
print(error[no_of_cycles])

snnsObject$saveNet(paste(basePath,"mlp_irisSnnsR_trained.net",sep=""),"mlp_irisSnnsR_trained")

snnsObject$saveNewPatterns(paste(basePath,"mlp_irisSnnsR.pat",sep=""), patset$set_no)




# pruning 


# assign parameters
first_error <- error[no_of_cycles]
max_pr_error_increase <- 10.0
pr_accepted_error <- 1.0
#pr_recreate <- TRUE
no_of_pr_retrain_cycles <- 1000
min_error_to_stop <- 0.01
init_matrix_value <- 1e-6
input_pruning <- TRUE
hidden_pruning <- TRUE


# maximum error
max_error = first_error * (1 + max_pr_error_increase / 100);
if(max_error < pr_accepted_error)
  max_error = pr_accepted_error

PR_ALL_PATTERNS <- -1
net_error <- vector()
retrain_loop_count <- 0

# pruning until net error <= maximum error
repeat {
  
  # delete some weights
  pr_res <- snnsObject$callPrunFunc(PR_ALL_PATTERNS)
  
  # calculate net error
  pr_res <- snnsObject$calcMeanDeviation(PR_ALL_PATTERNS) 
  net_error[1] <- pr_res[[2]]
  
  
  # retrain network
  if (net_error[1] > min_error_to_stop){ 
    
    
    for(j in 1:no_of_pr_retrain_cycles) {
      re_res <- snnsObject$learnAllPatternsFF(parameters)
      net_error[j] <- re_res[[2]]
      if (j %% 100 == 0) print(c(retrain_loop_count, (format(j,digits=2)), net_error[j]))
    }
    
    plot(net_error, type="l")
    
  }
  
  
  if (net_error[j] <= max_error) {
    retrain_loop_count <- retrain_loop_count + 1
  } else { 
    break 
  }
  
  
}

snnsObject$saveNet(paste(basePath,"mlp_irisSnnsR_pruned.net",sep=""),"mlp_irisSnnsR_pruned.net") 