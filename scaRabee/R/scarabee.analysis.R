
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

scarabee.analysis <- function(files=NULL,
                              method='population',
                              runtype=NULL,
                              debugmode=FALSE,
                              estim.options=NULL,
                              npts=NULL,
                              alpha=NULL,
                              solver.options=list(method='lsoda')
                              ) {
  
  oldwd <- getwd()
  
  anatry <- try({
  
    # Check files
    if (!is.null(files)){
      file.names <- c('data','param','model')
      if (is.list(files)){
        if (any(is.na(match(names(files),file.names))))
          stop(paste('files argument must be a list with the following ',
                     'items: \'data\', \'param\',\n  and \'model\'', sep=''))
        if (any(sapply(file,is.null))){
          stop('one or more element of the files argument are NULL.')
        }
      } else {
        stop('files argument must be a list.')
      }
    } else {
      stop('files argument cannot be NULL.')
    }
    if (!file.exists(files$data)){
      stop('the specified dataset file does not exist.')
    }
    if (!file.exists(files$model)){
      stop('the specified model file does not exist.')
    }
    
    # Check method
    if (is.null(method)){
      method <- 'population'
    } else {
      if (length(method)!=1 | !method%in%c('population','subject')){
        stop('method argument should be \'population\' or \'subject\'.')
      }
    }
    
    # Check runtype
    if (length(runtype)!=1 | 
        !runtype%in%c('estimation','simulation','gridsearch')){
      stop(paste('runtype argument must be \'estimation\', \'simulation\', ',
                 '\'gridsearch\'.',sep=''))
    }
      # Default is estimation (no gridsearch)
    issim <- 0
    isgrid <- 0
    if (runtype=='simulation') issim <- 1
    if (runtype=='gridsearch') {
      isgrid <- 1
      issim <- 1
      method <- 'population'
    }
    
    # Check debugmode
    if (!is.logical(debugmode)){
      stop('debugmode argument must either be a logical variable.')
    }
    
    # Check estim values
    if (!is.numeric(estim.options$maxiter) | !is.numeric(estim.options$maxfunc))
      stop(paste('estim.options$maxiter and estim.options$maxfunc should ',
                 'be numerical values.',sep=''))
    if (length(estim.options$maxiter)!=1 | length(estim.options$maxfunc)!=1)
      stop(paste('estim.options$maxiter and estim.options$maxfunc should ',
                 'only contain one value.',sep=''))
    
    # Check npts
    if (!is.null(npts)){
      if (!is.numeric(npts) | length(npts)!=1)
        stop('npts argument should be a single numerical value.')
      npts <-as.integer(npts)
    }
    
    # Check alpha
    if (!is.null(alpha)){
      if (any(!is.numeric(alpha)))
        stop('alpha argument is not numeric.')
      if (any(alpha<=1))
        stop('alpha cannot contain value(s) below 1.')
    }
    
    # Check solver.options
    if (is.null(solver.options)){
      solver.options <- list(method='lsoda')
    } else {
      if (!is.list(solver.options))
        stop('solver.options argument must be a list.')
      
      if (!'method'%in%names(solver.options))
        stop('solver.options argument must contain a \'method\' level.')
      
      if (is.null(solver.options$method))
      stop('method level of the solver.options argument cannot be NULL.')
      
    }
  ##############################################################################
    
    # 1- process mode file
    code <- scarabee.read.model(files=files, runtype=runtype)
      analysis <- code$name
      code$name <- NULL
    
    # create file names
    files$newdata <- NULL
    files$iter   <- paste(analysis,'.iterations.csv',sep='')
    files$report <- paste(analysis,'.report.txt',sep='')
    files$pred   <- paste(analysis,'.predictions.csv',sep='')
    files$est    <- paste(analysis,'.estimates.csv',sep='')
    files$sim    <- paste(analysis,'.simulations.csv',sep='')
    
    # 2- import parameter file (names + initial guess)
    param <- scarabee.read.parms(files=files)
    
    # 3- process data to be fitted
    tmp <- scarabee.read.data(files=files,
                              method=method)
    
    data <- tmp$data
    new.data.file <- tmp$new
    if (new.data.file) {
      files$newdata <- paste(files$data, '.new', sep='')
    }
    
    scarabee.check.reserved(names=param[,1],names(data[[1]][[1]]$cov)[-1])
    
    # 4- create working/backup directory
    scarabee.directory(curwd=getwd(),
                       files=files,
                       runtype=runtype,
                       analysis=analysis)
    
    #########################################################################
    # Definition of some variables to pass to model and estimation functions
    problem           <- list()
    problem$code      <- code
    problem$data      <- data
    problem$method    <- method
    problem$init      <- param
    problem$debugmode <- debugmode
    problem$modfun    <- code$template
    problem$solver.options <- solver.options
    
    # Divert the output if necessary
    if (!interactive())
      sink(file=paste(analysis,'.log',sep=''))
    
    if (isgrid == 1) {
      # Initializing report files
      initialize.report(problem=problem,param=param,files=files,isgrid=isgrid)
      
      # Evaluate grid
      cat(paste('* Direct grid search started at: ', Sys.time(),'\n\n',sep=''))
      fgrid <- scarabee.gridsearch(problem=problem,npts=npts,
                                   alpha=alpha,files=files)
      
      # Edition of the report file
      finalize.grid.report(problem=problem,fgrid=fgrid,files=files)
      
      # Update initial estimates in problem
      problem$init$value <- as.numeric(fgrid[1,1:(dim(fgrid)[2]-2)])
      
      cat('\n* Simulating model using the best solution from grid search\n\n')
    }
    
    if (issim == 0) {
      # Estimation run
      
      # Check if all parameters are fixed
      if (all(problem$init$isfix==1))
        stop(paste('model optimization cannot be performed ',
                   'when all parameters are fixed. Estimation aborted.',sep=''))
      
      # Check models for all subjects
      for (id in problem$data$ids){
        # Create subproblem for id
        idproblem <- problem[c('code','method','init','debugmode','modfun',
                               'solver.options')]
        idproblem$data <- problem$data[[id]]
        
        # check model
        scarabee.check.model(problem=idproblem,files=files)
      }
      
      # Initializing and preparing iteration log and report files
      initialize.report(problem=problem,param=param,files=files)
      
      # Fittings
      for (id in problem$data$ids){
        # Estimation start
        cat('\n','* Subject: ',id,'\n')
        cat('Optimization started\n')
        
        # Update report file
        write(paste('\n','* Subject: ',id,sep=''),
              file=files$report,
              append=TRUE,
              sep='\n')
        
        # Create subproblem for id
        idproblem <- problem[c('code','method','init','debugmode','modfun',
                               'solver.options')]
        idproblem$data <- problem$data[[id]]
        
        # parameter optimization
        Fit <- fitmle(problem=idproblem,estim.options=estim.options,files=files)
        
        # covariance step
        Fit <- fitmle.cov(problem=idproblem,Fit=Fit)
        
        # Edition of the report file
        Fit <- finalize.report(problem=idproblem,Fit=Fit,files=files)
        
        # Creation of a file containing final model predictions, residuals and
        # weighted residual for all states
        residual.report(problem=idproblem,Fit=Fit,files=files)
        
        cat('\nOptimization completed.\n')
      }
      # Creation of some diagnostic plots
      estimation.plot(problem=problem,Fit=Fit,files=files)
      
      cat('\nPlease, see report.txt file for more details.\n')
    
    } else {
      # Simulation run
      cat(paste('Simulation started at: ', Sys.time(),
                '. This operation can take a few momemt.\n',sep=''))
      
      # Evaluates the model given the parameter estimates and saves to file
      simdf <- simulation.report(problem=problem,files=files)
      
      # Plot simulations
      simulation.plot(problem=problem,simdf=simdf,files=files)
      
      cat(paste('\nSimulation complete at: ',Sys.time(),'\n',sep=''))
      cat(paste('Model predictions were saved into ',files$sim,'.\n\n',sep=''))
    }
    
    # Wrap-up things
    scarabee.clean(files=files,analysis=analysis)
    
    if (!interactive())
      sink()
  })
  
  if (class(anatry)=="try-error") setwd(oldwd)
  
}

