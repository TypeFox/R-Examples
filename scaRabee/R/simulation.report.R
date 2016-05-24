
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

simulation.report <- function(problem=NULL,files=NULL){
  
  issim <- 1
  
  # Create the matrix to be save to file
  simdf <- data.frame(ID=numeric(0),
                      TRT=numeric(0),
                      DVID=numeric(0),
                      TIME=numeric(0),
                      SIM=numeric(0),
                      OBS=numeric(0))
  
  ids <- problem$data$ids
  for (id in ids){
    # Evaluate model at each dose level
    trts <- problem$data[[id]]$trts
    fsim <- list()
    
    for (trt in as.character(trts)){
      # Create subproblem
      subproblem <- problem[c('code','method','init','debugmode','modfun',
                              'solver.options')]
      subproblem$data$xdata  <- sort(unique(problem$data[[id]][[trt]]$ana$TIME))
      subproblem$cov <- problem$data[[id]][[trt]]$cov
      subproblem$bolus <- problem$data[[id]][[trt]]$bolus
      subproblem$infusion <- problem$data[[id]][[trt]]$infusion
      
      # Retrieve primary parameters
      x <- subproblem$init
      parms <- c(get.parms.data(x=x,which='value',type='P'),
                 get.parms.data(x=x,which='value',type='L'),
                 get.parms.data(x=x,which='value',type='IC'),
                 get.parms.data(x=x,which='value',type='V'))
      names(parms) <- c(get.parms.data(x=x,which='names',type='P'),
                        get.parms.data(x=x,which='names',type='L'),
                        get.parms.data(x=x,which='names',type='IC'),
                        get.parms.data(x=x,which='names',type='V'))
      attr(parms,'type') <- c(get.parms.data(x=x,which='type',type='P'),
                              get.parms.data(x=x,which='type',type='L'),
                              get.parms.data(x=x,which='type',type='IC'),
                              get.parms.data(x=x,which='type',type='V'))
                      
      # Retrieve derived parameters
      if (subproblem$modfun%in%c('ode.model','dde.model')){
        derparms <- derived.parms(parms=parms,
                                  covdata=subproblem$cov,
                                  codederiv=subproblem$code$deriv,
                                  check=TRUE)
      } else {
        derparms <- NULL
      } 
      
      # Evaluate the model given the parameter estimates
      if (subproblem$modfun%in%c('ode.model','dde.model')){
        fsim[[trt]] <- do.call(eval(parse(text=subproblem$modfun)),
                               list(parms=parms,
                                    derparms=derparms,
                                    code=subproblem$code,
                                    bolus=subproblem$bolus,
                                    infusion=subproblem$infusion,
                                    xdata=subproblem$data$xdata,
                                    covdata=subproblem$cov,
                                    issim=issim,
                                    check=TRUE,
                                    options=subproblem$solver.options))
      } else {
        fsim[[trt]] <- do.call(eval(parse(text=subproblem$modfun)),
                               list(parms=parms,
                                    derparms=derparms,
                                    code=subproblem$code,
                                    bolus=subproblem$bolus,
                                    infusion=subproblem$infusion,
                                    xdata=subproblem$data$xdata,
                                    covdata=subproblem$cov,
                                    issim=issim,
                                    check=TRUE))
      }
    }
    
    # Check if there are predictions for all observation states 
    nstate <- size(fsim[[1]],1)-1
    obs.states <- unique(problem$data[[id]][[trt]]$ana$DVID)
    
    if (nstate < max(obs.states)){
      stop(paste('the number of model outputs is lower than',
                 ' the number of states in the\n  observation dataset.', 
                 ' Simulation aborted.', sep=''))
    }
    
    # Update data frame
    for (trt in as.character(trts)){
      times <- fsim[[trt]][1,]
      ntime <- length(times)
      for (cmt in 1:nstate){
        tmpdf <- data.frame(ID=rep(id,ntime),
                            TRT=rep(trt,ntime),
                            DVID=rep(cmt,ntime),
                            TIME=times,
                            SIM=fsim[[trt]][cmt+1,],
                            OBS=rep(NA,ntime))
        simdf <- rbind(simdf, tmpdf)
      } 
    }
    
    # Add observations to simdf
    for (trt in as.character(trts)){
      tmpdf <- problem$data[[id]][[trt]]$ana
      names(tmpdf) <- c('TIME','DVID','OBS')
      tmpdf$ID <- id ; tmpdf$TRT <- trt ; tmpdf$SIM <- NA
      tmpdf <- tmpdf[,c('ID','TRT','DVID','TIME','SIM','OBS')]
      
      simdf <- rbind(simdf, tmpdf)
    }
    
    # Report completion
    cat(sprintf('Subject %d completed.\n',id))
  }
  
  # Create or overwrite a simulation csv file
  write.table(simdf,
              file=files$sim,
              sep=',',
              na='NA',
              row.names=FALSE,
              quote=FALSE)
  
  return(simdf)
  
}

