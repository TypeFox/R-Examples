#############################################################################
#
# Copyright Alexandru Olteanu, Patrick Meyer and Sébastien Bigaret, 2015
#
# Contributors:
#   Alexandru Olteanu <al.olteanu@gmail.com>
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sébastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
#  	
# This software, MCDA, is a package for the R statistical software which 
# allows to use MCDA algorithms and methods. 
# 
# This software is governed by the CeCILL license (v2) under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/ or redistribute the software under the terms of the
# CeCILL license as circulated by CEA, CNRS and INRIA at the following
# URL "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#		
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#		
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

MRSortInferenceApprox <- function(performanceTable, assignments, categoriesRanks, criteriaMinMax, alg_total_time = 90, alg_repeats = 3, alg_repeat_time = 30, 
                                         alg_repeat_iterations = 30,mh_max_temp_step = 0.2, mh_min_temp_step = 0.02, mh_temp_step_increase = 1.25,
                                         mh_temp_step_decrease = 0.8, veto = FALSE, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## check the input data
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.vector(assignments)))
    stop("assignments should be a vector")
  
  if (!(is.vector(categoriesRanks)))
    stop("categoriesRanks should be a vector")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!is.logical(veto))
    stop("veto should be a boolean")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be a vector")
  
  if (!is.numeric(alg_total_time))
    stop("alg_total_time should be numeric")
  else if (alg_total_time%%1!=0)
    stop("alg_total_time should be an integer")
  else if (alg_total_time<=0)
    stop("alg_total_time should be strictly pozitive")
  
  if (!is.numeric(alg_repeats))
    stop("alg_repeats should be numeric")
  else if (alg_repeats%%1!=0)
    stop("alg_repeats should be an integer")
  else if (alg_repeats<=0)
    stop("alg_repeats should be strictly pozitive")
  
  if (!is.numeric(alg_repeat_time))
    stop("alg_repeat_time should be numeric")
  else if (alg_repeat_time%%1!=0)
    stop("alg_repeat_time should be an integer")
  else if (alg_repeat_time<=0)
    stop("alg_repeat_time should be strictly pozitive")

  if (!is.numeric(alg_repeat_iterations))
    stop("alg_repeat_iterations should be numeric")
  else if (alg_repeat_iterations%%1!=0)
    stop("alg_repeat_iterations should be an integer")
  else if (alg_repeat_iterations<=0)
    stop("alg_repeat_iterations should be strictly pozitive")
    
  if (!is.numeric(mh_max_temp_step))
    stop("mh_max_temp_step should be numeric")
  else if (mh_max_temp_step<=0)
    stop("mh_max_temp_step should be strictly pozitive")
  
  if (!is.numeric(mh_min_temp_step))
    stop("mh_min_temp_step should be numeric")
  else if (mh_min_temp_step<=0)
    stop("mh_min_temp_step should be strictly pozitive")
  
  if (!is.numeric(mh_temp_step_increase))
    stop("mh_temp_step_increase should be numeric")
  else if (mh_temp_step_increase<=1)
    stop("mh_temp_step_increase should be strictly above 1")
  
  if (!is.numeric(mh_temp_step_decrease))
    stop("mh_temp_step_decrease should be numeric")
  else if (mh_temp_step_decrease<=0 || mh_temp_step_decrease>=1)
    stop("mh_temp_step_decrease should be betweem 0 and 1")
  
  ## filter the data according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
    assignments <- assignments[alternativesIDs]
  } 
  
  if (!is.null(criteriaIDs)){
    performanceTable <- performanceTable[,criteriaIDs]
    criteriaMinMax <- criteriaMinMax[criteriaIDs]
  }
  
  # data is filtered, check for some data consistency
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (is.null(dim(performanceTable))) 
    stop("less than 2 criteria or 2 alternatives")
  
  # -------------------------------------------------------
  
  # init total starting time
  
  start.time.total <- Sys.time()
  
  # get number of alternatives, criteria and categories
  
  numAlt <- dim(performanceTable)[1]
  
  numCrit <- dim(performanceTable)[2]
  
  numCat <- length(categoriesRanks)
  
  # initialize model parameters
  
  model.params <- InitModel(performanceTable, assignments, categoriesRanks, criteriaMinMax)
  
  # initialize best parameters
  
  best.params <- list(gamma = model.params$gamma, lambda = model.params$lambda, weights = model.params$weights, profilesPerformances = model.params$profilesPerformances, vetoPerformances = model.params$vetoPerformances)
  
  best.fitness <- Fitness(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, best.params)
  
  maxed.fitness <- FALSE
  
  # check if we've spent all the allocated time
  
  time.taken <- Sys.time() - start.time.total
  
  if (time.taken >= alg_total_time)
    return(best.params)
  
  # repeat the algorithm several times
  
  for (i in 1:alg_repeats)
  {
    # init local starting time
    
    start.time.local <- Sys.time()
    
    # start with best parameters
    
    model.params <- list(gamma = best.params$gamma, lambda = best.params$lambda, weights = best.params$weights, profilesPerformances = best.params$profilesPerformances, vetoPerformances = best.params$vetoPerformances)
    
    # initialize temperature
    
    temp.step <- mh_max_temp_step
    
    # go through the algorithm iterations
    
    for (j in 1:alg_repeat_iterations)
    {
      # get lambda and weights
      
      temp.params <- InferLW(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params)
      
      if(!is.null(temp.params))
      {
        model.params$lambda <- temp.params$lambda
        
        model.params$weights <- temp.params$weights
      }
      
      # get profiles
      
      temp.params <- InferP(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params, temp.step)
      
      if(!is.null(temp.params))
      {
        model.params$profilesPerformances <- temp.params$profilesPerformances
        
        model.params$vetoPerformances <- temp.params$vetoPerformances
      }
      
      # evaluate parameters
      
      fitness <- Fitness(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params)
      
      print(c(fitness,best.fitness))
      
      # check for overall improvement
      
      if (fitness >= best.fitness)
      {
        # increase temperature step so that the MH will do less iterations
        
        temp.step <- temp.step * mh_temp_step_increase
        
        # stay within limits
        
        if (temp.step > mh_max_temp_step)
          temp.step <- mh_max_temp_step
        
        # update model parameters only if fitness has improved or randomly if it remained stagnant
        
        if (fitness > best.fitness || sample(c(TRUE,FALSE),1))
        {          
          best.params <- list(gamma = model.params$gamma, lambda = model.params$lambda, weights = model.params$weights, profilesPerformances = model.params$profilesPerformances, vetoPerformances = model.params$vetoPerformances)
          
          # record best fitness
          
          best.fitness <- fitness
        }
        
        # check if we've maxed out the fitness -> algorithm stops
        
        if (best.fitness == 1.0)
          maxed.fitness <- TRUE
      }
      else
      {
        # decrease temperature step so that the MH will do more iterations
        
        temp.step <- temp.step * mh_temp_step_decrease
        
        # stay within limits
        
        if (temp.step < mh_min_temp_step)
          temp.step <- mh_min_temp_step
      }
      
      # check if we've spent all the allocated time for this repeat
      
      time.taken <- Sys.time() - start.time.local
      
      if (time.taken >= alg_repeat_time)
        break
      
      # check if we've maxed out the fitness -> algorithm stops
      
      if (maxed.fitness)
        break
    }
    
    # check if we've spent all the allocated time
    
    time.taken <- Sys.time() - start.time.total
    
    if (time.taken >= alg_total_time)
      break
    
    # check if we've maxed out the fitness -> algorithm stops
    
    if (maxed.fitness)
      break
  }
  
  # add bottom profile
  
  bottomprofile = rep(-Inf,numCrit)
  
  for (i in 1:numCrit)
    if(criteriaMinMax[i] == "min")
      bottomprofile[i] <- Inf
  
  best.params$profilesPerformances <- rbind(model.params$profilesPerformances,bottomprofile)
  
  best.params$vetoPerformances <- rbind(model.params$vetoPerformances,bottomprofile)
  
  rownames(best.params$profilesPerformances) <- names(categoriesRanks)
  
  rownames(best.params$vetoPerformances) <- names(categoriesRanks)
  
  # return result
  
  return(best.params)
}

InitModel <-function(performanceTable, assignments, categoriesRanks, criteriaMinMax){
  
  # get number of alternatives, criteria and categories
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  # init parameters
  
  model.params <- list(gamma = 0.001, lambda = 0.5, weights = rep(1/numCrit, times = numCrit), profilesPerformances = matrix(0,numCat-1,numCrit), vetoPerformances = matrix(0,numCat-1,numCrit))
  
  colnames(model.params$profilesPerformances) <- colnames(performanceTable)
  
  colnames(model.params$vetoPerformances) <- colnames(performanceTable)
  
  # init vetoes
  
  for (j in 1:numCrit)
  {
    # get criterion preference direction
    if(criteriaMinMax[j] == "max")
      model.params$vetoPerformances[,j] <- rep(apply(performanceTable, 2, min)[j] - model.params$gamma, times = numCat-1)
    else
      model.params$vetoPerformances[j] <- rep(apply(performanceTable, 2, max)[j] + model.params$gamma, times = numCat-1)
  }
  
  # go thorough each criterion
  
  for (j in 1:numCrit)
  {
    # get criterion preference direction
    
    critdir <- 1
    
    if(criteriaMinMax[j] == "min")
      critdir <- -1
    
    # get all values and a list of the categories in which the alternatives containing that value are assigned
    
    values <- c()
    
    valuecategories <- list()
    
    for (i in 1:numAlt)
    {
      if(!(performanceTable[i,j] %in% values))
      {
        values <- c(values,performanceTable[i,j])
        
        valueindex <- match(performanceTable[i,j],values)
        
        valuecategories[[valueindex]] <- c(categoriesRanks[assignments[i]])
      }
      else
      {
        valueindex <- match(performanceTable[i,j],values)
        
        valuecategories[[valueindex]] <- c(valuecategories[[valueindex]],categoriesRanks[assignments[i]])
      }
    }
    
    # order the values from worst to best
    
    valuecategories <- valuecategories[order(critdir * values)]
    
    values <- values[order(critdir *values)]
    
    # get profiles values
    
    startvalindex <- 1
    
    # go from the worst profile to the best
    
    for (i in (numCat-1):1)
    {
      # are there still values to explore?
      
      if(startvalindex <= length(values))
      {
        # get current value
        
        value <- values[startvalindex]
        
        # find its initial fitness for the current profile
        
        f <- 0
        
        # go through values below the one at startvalindex
        
        if(startvalindex > 1)
          for(k in 1:(startvalindex-1))
            for(l in valuecategories[[k]])
            {
              # all values belonging to alternatives that are classified in a category below the profile (index is higher than that of the profile) affect pozitively the fitness; all others negatively
              
              if(l > i)
                f <- f + 1
              else
                f <- f - 1
            }
        
        if(startvalindex <= length(values))
          for(k in startvalindex:length(values))
            for(l in valuecategories[[k]])
            {
              # all values belonging to alternatives that are classified in a category above the profile (index is lower or equal to that of the profile) affect pozitively the fitness; all others negatively
              
              if(l <= i)
                f <- f + 1
              else
                f <- f - 1
            }
        
        # go through following values and update f
        
        newf <- f
        
        currentvalindex <- startvalindex
        
        if(startvalindex < length(values))
          for(k in (startvalindex+1):length(values))
          {
            # update f
            
            for(l in valuecategories[[k]])
            {
              if(l > i)
                newf <- newf + 1
              else
                newf <- newf - 1
            }
            
            # if f is improved we store it and the value that the profile should take
            
            if(newf > f)
            {
              f <- newf
              
              value <-values[k]
              
              currentvalindex <- k
            }
          }
        
        # same as above but for when we have reached the last value for this criterion and we consider a value above
        
        for(l in valuecategories[[length(values)]])
        {
          if(l > i)
            newf <- newf + 1
          else
            newf <- newf - 1
        }
        
        # if f is improved we store it and the value that the profile should take
        
        if(newf > f)
        {
          f <- newf
          
          value <- values[length(values)] + critdir * model.params$gamma
          
          currentvalindex <- length(values)
        }
        
        # set profile value
        
        model.params$profilesPerformances[i,j] <- value
        
        # set the starting value index for the next profile
        
        startvalindex <- currentvalindex# + 1
      }
      else
      {
        # we have are at the top of the scale
        
        model.params$profilesPerformances[i,j] <- values[length(values)] + critdir * model.params$gamma
      }
    }
  }
  return(model.params)
}

InferLW <-function(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params){

  # get number of alternatives, criteria and categories
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  # take out alternatives with vetoes
  
  model.discordance <- GetDiscordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
  
  alternatives <- c(1:numAlt)[model.discordance %==% 0]
  
  numAlt <- length(alternatives)
  
  if(numAlt == 0)
    return(NULL)
  
  # get temp directory
  
  tempPath <- tempdir()
  
  # get model file
  
  modelFile <- system.file("extdata","MRSortInferenceLW.gmpl", package="MCDA")
  
  # create temporary data file
  
  dataFile <- tempfile()
  
  # copy file to temp directory
  
  file.copy(modelFile, dataFile)
  
  # open writing channel
  
  sink(dataFile, append=TRUE)
  
  # write data
  
  cat("data;\n")
  cat("param X := ")
  cat(numAlt)
  cat(";\n\n")
  
  cat("param F := ")
  cat(numCrit)
  cat(";\n\n")
  
  cat("param lClow : ")
  cat(1:numCrit)
  cat(" := \n")
  for (i in 1:numAlt){
    cat(i)
    cat("\t")
    for (j in 1:numCrit)
    {
      critdir <- 1
      if(criteriaMinMax[j] == "min")
        critdir <- -1
      categ <- categoriesRanks[assignments[alternatives[i]]]
      if(categ %==% numCat)
        cat("1")
      else
      {
        if((critdir * performanceTable[alternatives[i],j]) %>=% (critdir * model.params$profilesPerformances[categ,j]))
          cat("1")
        else
          cat("0")
      }
      if(j != numCrit)
        cat("\t")
      else
      {
        if(i != numAlt)
          cat("\n")
        else
          cat(";\n\n")
      }
    }
  }
  
  cat("param lCupp : ")
  cat(1:numCrit)
  cat(" := \n")
  for (i in 1:numAlt){
    cat(i)
    cat("\t")
    for (j in 1:numCrit)
    {
      critdir <- 1
      if(criteriaMinMax[j] == "min")
        critdir <- -1
      categ <- categoriesRanks[assignments[alternatives[i]]]
      if(categ %==% 1)
        cat("0")
      else
      {
        if((critdir * performanceTable[alternatives[i],j]) %>=% (critdir * model.params$profilesPerformances[categ-1,j]))
          cat("1")
        else
          cat("0")
      }
      if(j != numCrit)
        cat("\t")
      else
      {
        if(i != numAlt)
          cat("\n")
        else
          cat(";\n\n")
      }
    }
  }
  
  cat("param gamma:=")
  cat(model.params$gamma)
  cat(";\n")
  
  cat("end;\n")
  
  sink()
  
  lp<-initProbGLPK()
  
  tran<-mplAllocWkspGLPK()
  
  setMIPParmGLPK(PRESOLVE, GLP_ON)
  
  termOutGLPK(GLP_OFF)
  
  out<-mplReadModelGLPK(tran, dataFile, skip=0)
  
  if (is.null(out))
    out <- mplGenerateGLPK(tran)
  else 
    stop(return_codeGLPK(out))
  
  if (is.null(out))
    mplBuildProbGLPK(tran,lp)
  else 
    stop(return_codeGLPK(out))
  
  # solve the problem
  
  solveMIPGLPK(lp)
  
  if(mipStatusGLPK(lp)==5){
    mplPostsolveGLPK(tran, lp, sol = GLP_MIP)
    
    solution <- mipColsValGLPK(lp)
    
    # get results
    
    varnames <- c()
    
    for (i in 1:length(solution))
      varnames <- c(varnames,getColNameGLPK(lp,i))
    
    lambda <- solution[varnames=="lambda"]
    
    weightsnames <- c()
    
    for (i in 1:numCrit)
    {
      weightsnames <- c(weightsnames,paste("w[",i,"]",sep=""))
    }
    
    weights <- c()
    
    for (i in 1:numCrit)
      weights <- c(weights,solution[varnames==weightsnames[i]])
    
    return(list(lambda = lambda, weights = weights))
    
  }
  else
    return(NULL)
}

InferP <-function(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params, temp.step){
  
  # get number of alternatives, criteria and categories
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  # init assignments and concordance
  
  model.concordance <- GetConcordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
  
  model.discordance <- GetDiscordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
  
  model.assignments <- GetAssignments(model.concordance, model.discordance, model.params)
  
  # init temperature
  
  t <- 1.0
  
  while(t > 0)
  {
    # go through each profile at random
    
    for(k in sample(1:(numCat-1)))
    {
      # go through each criterion at random
      
      for(j in sample(1:numCrit))
      {
        # get range within which the profile can move
        
        valmin <- apply(performanceTable, 2, min)[j] - model.params$gamma
        
        valmax <- apply(performanceTable, 2, max)[j] + model.params$gamma
        
        if(criteriaMinMax[j] == "max")
        {
          if(k > 1)
            valmax <- model.params$profilesPerformances[k-1,j]
          
          if(k < numCat - 1)
            valmin <- model.params$profilesPerformances[k+1,j]
          
          if(valmin %<=% model.params$vetoPerformances[k,j])
            valmin <- model.params$vetoPerformances[k,j] + model.params$gamma
        }
        else
        {
          if(k > 1)
            valmin <- model.params$profilesPerformances[k-1,j]
          
          if(k < numCat - 1)
            valmax <- model.params$profilesPerformances[k+1,j]
          
          if(valmax %>=% model.params$vetoPerformances[k,j])
            valmax <- model.params$vetoPerformances[k,j] - model.params$gamma
        }
        
        # get new value
        
        val <- model.params$profilesPerformances[k,j]
        
        h <- c(-Inf,-Inf)
        
        # try several random values and select the one with maximum heuristic value
        
        for(i in 1:10)
        {
          newval <- runif(1, valmin, valmax)
          
          newh <- Heuristic(k, j, newval, performanceTable, assignments, categoriesRanks, criteriaMinMax, model.assignments, model.concordance, model.discordance, model.params)
          
          if(newh[1] > h[1] || (newh[1] %==% h[1] && newh[2] > h[2]) || (newh %==% h && sample(c(TRUE,FALSE),1)))
          {
            h <- newh
            
            val <- newval
          }
        }
        
        # simulated annealing condition for accepting the change
        
        if(h[1] > 0 || (h[1] %==% 0 && h[2] > 0) || runif(1,0,1) < exp(-1/t))
        {
          # update profile value
          
          model.params$profilesPerformances[k,j] <- val
          
          # update assignments and concordance (could be done more smartly as only a column in the concordance matrix changes)
          
          model.concordance <- GetConcordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
          
          model.assignments <- GetAssignments(model.concordance, model.discordance, model.params)
        }
      }
    }

    if(veto)
    {
      for(k in sample(1:(numCat-1)))
      {
        # go through each criterion at random
        
        for(j in sample(1:numCrit))
        {
          # get range within which the profile can move
          
          valmin <- apply(performanceTable, 2, min)[j] - model.params$gamma
          
          valmax <- apply(performanceTable, 2, max)[j] + model.params$gamma
          
          if(criteriaMinMax[j] == "max")
          {
            if(k > 1)
              valmax <- model.params$vetoPerformances[k-1,j]
            
            if(k < numCat - 1)
              valmin <- model.params$vetoPerformances[k+1,j]
            
            if(valmax %>=% model.params$profilesPerformances[k,j])
              valmax <- model.params$profilesPerformances[k,j] - model.params$gamma
          }
          else
          {
            if(k > 1)
              valmin <- model.params$vetoPerformances[k-1,j]
            
            if(k < numCat - 1)
              valmax <- model.params$vetoPerformances[k+1,j]
            
            if(valmin %<=% model.params$profilesPerformances[k,j])
              valmin <- model.params$profilesPerformances[k,j] + model.params$gamma
          }

          # get new value
          
          val <- model.params$vetoPerformances[k,j]
          
          h <- c(-Inf,-Inf)
          
          # try several random values and select the one with maximum heuristic value
          
          for(i in 1:10)
          {
            newval <- runif(1, valmin, valmax)
            
            newh <- HeuristicV(k, j, newval, performanceTable, assignments, categoriesRanks, criteriaMinMax, model.assignments, model.concordance, model.discordance, model.params)
            
            if(newh[1] > h[1] || (newh[1] %==% h[1] && newh[2] > h[2]) || (newh %==% h && sample(c(TRUE,FALSE),1)))
            {
              h <- newh
              
              val <- newval
            }
          }
          
          # simulated annealing condition for accepting the change
          
          if(h[1] > 0 || (h[1] %==% 0 && h[2] > 0) || runif(1,0,1) < exp(-1/t))
          {
            # update profile value
            
            model.params$vetoPerformances[k,j] <- val
            
            # update assignments and concordance (could be done more smartly as only a column in the concordance matrix changes)
            
            model.concordance <- GetConcordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
            
            model.disncordance <- GetDiscordance(performanceTable, categoriesRanks, criteriaMinMax, model.params)
            
            model.assignments <- GetAssignments(model.concordance, model.discordance, model.params)
          }
        }
      }
    }
    
    # reduce temperature
    
    t <- t - temp.step
  }
  
  return(list(profilesPerformances = model.params$profilesPerformances, vetoPerformances = model.params$vetoPerformances))
}

GetAssignments <-function(model.concordance, model.discordance, model.params){
  
  # get a list of alternatives assignments to categories using the given model parameters
  
  numAlt <- dim(model.concordance)[1]
  
  numCat <- dim(model.concordance)[2] + 1
  
  model.assignments = c(rep(1,times=numAlt))
  
  for(i in 1:numAlt)
  {
    for(k in (numCat-1):1)
    {
      # compare support to majority threshold
      if(!(model.concordance[i,k] %>=% model.params$lambda) || model.discordance[i,k] > 0)
      {
        # insufficient support -> k is the upper profile of the category we're in
        
        model.assignments[i] <- k + 1
        
        break
      }
    }
  }
  
  return(model.assignments)
}

GetConcordance <-function(performanceTable, categoriesRanks, criteriaMinMax, model.params){
  
  # get the overall concordance indices for each alternative and each profile
  
  numAlt <- dim(performanceTable)[1]
  
  numCrit <- dim(performanceTable)[2]
  
  numCat <- length(categoriesRanks)
  
  model.concordance = matrix(0,numAlt,numCat-1)
  
  for(i in 1:numAlt)
  {
    for(k in 1:(numCat-1))
    {
      for(j in 1:numCrit)
      {
        critdir <- 1
        
        if(criteriaMinMax[j] == "min")
          critdir <- -1
        
        if((critdir * performanceTable[i,j][[1]]) %>=% (critdir * model.params$profilesPerformances[k,j][[1]]))
          model.concordance[i,k] <- model.concordance[i,k] + model.params$weights[j]
      }
    }
  }
  return(model.concordance)
}

GetDiscordance <-function(performanceTable, categoriesRanks, criteriaMinMax, model.params){
  
  # get the overall concordance indices for each alternative and each profile
  
  numAlt <- dim(performanceTable)[1]
  
  numCrit <- dim(performanceTable)[2]
  
  numCat <- length(categoriesRanks)
  
  model.discordance = matrix(0,numAlt,numCat-1)
  
  for(i in 1:numAlt)
  {
    for(k in 1:(numCat-1))
    {
      for(j in 1:numCrit)
      {
        critdir <- 1
        
        if(criteriaMinMax[j] == "min")
          critdir <- -1
        if((critdir * performanceTable[i,j][[1]]) %<=% (critdir * model.params$vetoPerformances[k,j][[1]]))
          model.discordance[i,k] <- model.discordance[i,k] + 1
      }
    }
  }
  return(model.discordance)
}

Heuristic <-function(k, j, value ,performanceTable, assignments, categoriesRanks, criteriaMinMax, model.assignments, model.concordance, model.discordance, model.params){
  
  # get number of alternatives, criteria and categories
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  # get criterion preference direction
  
  critdir <- 1
  
  if(criteriaMinMax[j] == "min")
    critdir <- -1
  
  # init heuristic
  
  h = c(0,0)
  
  # go thourough each alternative
  
  for(i in 1:numAlt)
  {
    # get categories  to which it is assigned by the DM and by the model
    
    given.k <- categoriesRanks[assignments[i]][[1]]
    
    found.k <- model.assignments[i]
    
    # get object, old profile and new profile values multiplying by critdir in order to use one set of conditions for both cases
    
    old.val <- critdir * model.params$profilesPerformances[k,j]
    
    new.val <- critdir * value
    
    obj.val <- critdir * performanceTable[i,j]
    
    # object misclassified in k or above instead of k + 1 -> model.concordance >= l and no veto
    
    if(given.k %==% (k + 1) && found.k %<=% k)
    {
      # moving profile above object corrects classification
      
      if(new.val > obj.val && obj.val %>=% old.val && !((model.concordance[i,k] - model.params$weights[j]) %>=% model.params$l))
        h[1] <- h[1] + 1
      
      # moving profile above object does not improve classification but reduces concordance
      
      if(new.val > obj.val && obj.val %>=% old.val && (model.concordance[i,k] - model.params$weights[j]) %>=% model.params$l)
        h[2] <- h[2] + 1
      
      # moving profile below object does not improve classification and increases concordance -> maybe add a third component to the fitness
    }
    
    # object misclassified in k + 1 or below instead of k -> model.concordance < l or veto
    
    if(given.k %==% k && found.k %>=% (k + 1))
    {
      # if object misclassified due to veto then nothing can be done here
      if(model.discordance[i,k] %==% 0)
      {      
        # moving profile below object corrects classification
        
        if(old.val > obj.val && obj.val %>=% new.val && (model.concordance[i,k] + model.params$weights[j]) %>=% model.params$l)
          h[1] <- h[1] + 1
        
        # moving profile below object does not improve classification but increases concordance
        if(old.val > obj.val && obj.val %>=% new.val && !((model.concordance[i,k] + model.params$weights[j]) %>=% model.params$l))
          h[2] <- h[2] + 1
        
        # moving profile above object does not improve classification and decreases concordance -> maybe add a third component to the fitness
      }
    }
    
    # object correctly classified in k + 1 -> model.concordance < l or veto
    
    if(given.k %==% (k + 1) && found.k %==% (k + 1))
    {
      if(model.discordance[i,k] %==% 0)
      {      
        # moving profile below object results in misclassification
        if(old.val > obj.val && obj.val %>=% new.val && (model.concordance[i,k] + model.params$weights[j]) %>=% model.params$l)
          h[1] <- h[1] - 1
        
        # moving profile below object keeps correct classification but increases concordance
        if(old.val > obj.val && obj.val %>=% new.val && !((model.concordance[i,k] + model.params$weights[j]) %>=% model.params$l))
          h[2] <- h[2] - 1
        
        # moving profile above object keeps correct classification and decreases concordance
      }
    }
    
    # object correctly classified in k -> model.concordance >= l and no veto
    
    if(given.k %==% k && found.k %==% k)
    {
      # moving profile above object results in misclassification
      
      if(new.val > obj.val && obj.val %>=% old.val && (model.concordance[i,k] - model.params$weights[j]) %>=% model.params$l)
        h[1] <- h[1] - 1
      
      # moving profile above object keeps correct classification but decreases concordance
      
      if(new.val > obj.val && obj.val %>=% old.val && !((model.concordance[i,k] - model.params$weights[j]) %>=% model.params$l))
        h[2] <- h[2] - 1
      
      # moving profile below object keeps correct classification and increases concordance
    }
  }
  return(h)
}

HeuristicV <-function(k, j, value ,performanceTable, assignments, categoriesRanks, criteriaMinMax, model.assignments, model.concordance, model.discordance, model.params){
  
  # get number of alternatives, criteria and categories
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  # get criterion preference direction
  
  critdir <- 1
  
  if(criteriaMinMax[j] == "min")
    critdir <- -1
  
  # init heuristic
  
  h = c(0,0)
  
  # go thourough each alternative
  
  for(i in 1:numAlt)
  {
    # get categories  to which it is assigned by the DM and by the model
    
    given.k <- categoriesRanks[assignments[i]][[1]]
    
    found.k <- model.assignments[i]
    
    # get object, old profile and new profile values multiplying by critdir in order to use one set of conditions for both cases
    
    old.val <- critdir * model.params$vetoPerformances[k,j]
    
    new.val <- critdir * value
    
    obj.val <- critdir * performanceTable[i,j]
    
    # object misclassified in k or above instead of k + 1 -> model.concordance >= l and no veto
    
    if(given.k %==% (k + 1) && found.k %<=% k)
    {
      # moving profile above object corrects classification
      
      if(new.val %>=% obj.val && obj.val > old.val)
        h[1] <- h[1] + 1
      
      # moving profile below object does not improve classification
    }
    
    # object misclassified in k + 1 or below instead of k -> model.concordance < l or veto
    
    if(given.k %==% k && found.k %>=% (k + 1))
    {   
      # moving profile below object corrects classification
      
      if(old.val %>=% obj.val && obj.val > new.val && model.concordance[i,k] %>=% model.params$l && model.discordance[i,k] %==% 1)
        h[1] <- h[1] + 1
      
      # moving profile below object does not improve classification but decreases vetoes
      if(old.val %>=% obj.val && obj.val > new.val && (model.discordance[i,k] > 1 || (!(model.concordance[i,k] %>=% model.params$l) && model.discordance[i,k] %==% 1)))
        h[2] <- h[2] + 1
      
      # moving profile above object does not improve classification and increases discordance
      if(new.val %>=% obj.val && obj.val > old.val)
        h[2] <- h[2] - 1
    }
    
    # object correctly classified in k + 1 -> model.concordance < l or veto
    
    if(given.k %==% (k + 1) && found.k %==% (k + 1))
    {
      # moving profile below object results in misclassification
      if(old.val %>=% obj.val && obj.val > new.val && model.concordance[i,k] %>=% model.params$l && model.discordance[i,k] %==% 1)
        h[1] <- h[1] - 1
      
      # moving profile below object keeps correct classification and reduces discordance
      if(old.val %>=% obj.val && obj.val > new.val && (!(model.concordance[i,k] %>=% model.params$l) || model.discordance[i,k] > 1))
        h[2] <- h[2] + 1
      
      # moving profile above object keeps correct classification but increases concordance
      if(new.val %>=% obj.val && obj.val > old.val)
        h[2] <- h[2] - 1
    }
    
    # object correctly classified in k -> model.concordance >= l and no veto
    
    if(given.k %==% k && found.k %==% k)
    {
      # moving profile above object results in misclassification
      
      if(new.val %>=% obj.val && obj.val > old.val)
        h[1] <- h[1] - 1
    }
  }
  return(h)
}

GetCategory <- function(performance, criteriaMinMax, veto, model.params){
  
  numCrit <- length(performance)
  
  numCat <- dim(model.params$profilesPerformances)[1] + 1
  
  # go through profiles from the worst (highest index) to the best (lowest index)
  
  for(k in (numCat-1):1)
  {
    # compute relation
    
    C <- 0
    
    V <- FALSE
    
    for(j in 1:numCrit)
    {
      critdir <- 1
      
      if(criteriaMinMax[j] == "min")
        critdir <- -1
      
      if((critdir * performance[j]) %>=% (critdir * model.params$profilesPerformances[k,j]))
        C <- C + model.params$weights[j]
      
      if(veto)
        if((critdir * performance[j]) %<=% (critdir * model.params$vetoPerformances[k,j]))
        {
          V <- TRUE
          break
        }
    }
    # compare support to majority threshold
    if(!(C %>=% model.params$lambda) || V)
    {
      # insufficient support -> k is the upper profile of the category we're in
      
      return(k + 1)
    }
  }
  return(1)
}

Fitness <-function(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto, model.params){
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  ca <- 0
  
  # go through each alternative
  
  for(i in 1:numAlt)
  {
    catfound <- GetCategory(performanceTable[i,], criteriaMinMax, veto, model.params)
    
    if(categoriesRanks[assignments[i]] == catfound)
      ca <- ca + 1
  }
  
  ca <- ca / numAlt
  
  return(ca)
}
