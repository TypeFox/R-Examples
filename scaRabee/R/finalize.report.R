
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

finalize.report <- function(problem=NULL,Fit=NULL,files=NULL){
  
# Adds some termination information
  # Calculates the dimension of the longest parameter name
  maxlen <- max(9,max(nchar(problem$init$names)))

  # Builds and writes the first 2 lines of the estimation table title
  tmp <- sprintf('Akaike Information Criterium (AIC): %s\n',
                 sprintf('%f',Fit$AIC))
  tmp <- c(tmp,
           sprintf('Total number of parameters: %d',
                   length(problem$init$names)))
  tmp <- c(tmp,
           sprintf('Number of estimated parameters: %d\n',
                   length(Fit$estimations)))
  tmp <- c(tmp,'\nPrimary parameter(s)\n')
  tmp <- c(tmp,
           sprintf(paste(c(rep(' ',maxlen+6),'%s       %s'),
                         collapse=''),'Initial','Final'))
  tmp <- c(tmp,
           sprintf(paste(c(' %s',rep(' ',maxlen-9+6),'%s        ',
                           '%s         ','%s      %s'),collapse=''),
                   'Parameter','value','value','CV%','Confidence interval (95%)'))

  write(tmp,file=files$report,append=TRUE,sep='\n')
  
# Need to reorder the fixed parameters as fitmle_cov did for estimated
# parameters (Fit$orderedestimations)
  # Determines which parameters were estimated or fixed
  ordered <- order.parms.list(x=problem$init)
  estparam <- problem$init[which(problem$init$isfix==0),]
  fixparam <- problem$init[which(problem$init$isfix==1),]
  estorder <- ordered[which(ordered$isfix==0),]
  fixorder <- ordered[which(ordered$isfix==1),]

  # Calculates the number of fixed model and variance parameters
    # p = nb of model parametres
  p <- length(get.parms.data(x=fixparam,which='type',type='P')) +
       length(get.parms.data(x=fixparam,which='type',type='L')) +
       length(get.parms.data(x=fixparam,which='type',type='IC'))
    # q = nb of variance parametres
  q <- length(get.parms.data(x=fixparam,which='type',type='V'))
  
  # Determines in fixparam the corresponding parameter indices from fixorder
  indices <- NULL
  if ((p+q)!=0){
    for (i in 1:(p+q)){
      for (j in 1:(p+q)){
         if (fixparam$names[j]==fixorder$names[i]){
           indices <- c(indices,j)
         }
      }
    }
  }
  
  # Gets the reordered parameter value and names
  if (is.null(indices)){
    Fit$orderedfixed <- data.frame(names=character(0),
                                   value=numeric(0),
                                   type =character(0),
                                   isfix=numeric(0),
                                   stringsAsFactors=FALSE)
  } else {
    Fit$orderedfixed <- data.frame(names=fixparam$names[indices],
                                   value=fixparam$value[indices],
                                   type =fixparam$type[indices],
                                   isfix=rep(1,length(indices)),
                                   stringsAsFactors=FALSE)
  }
 
# Need to reorder the initial data structure for reporting reasons
  orderparnames <- c(Fit$orderedestimations$names,Fit$orderedfixed$names)
  indices <- c()
  for (i in 1:length(problem$init$names)){
    for (j in 1:length(problem$init$names)){
       if (problem$init$names[j]==orderparnames[i]){
         indices <- c(indices,j)
       }
    }    
  }
  Fit$orderedinitial <- data.frame(names=problem$init$names[indices],
                                   type =problem$init$type[indices],
                                   value=problem$init$value[indices],
                                   isfix=problem$init$isfix[indices],
                                   stringsAsFactors=FALSE)
 
# Builds and writes each line of the estimation table
  # first the block of estimated parameters
    # nb of model parametres
  nemodpar <- length(get.parms.data(x=Fit$orderedestimations,which='type',type='P'))+
              length(get.parms.data(x=Fit$orderedestimations,which='type',type='L'))+
              length(get.parms.data(x=Fit$orderedestimations,which='type',type='IC'))
    # nb of variance parametres
  nevarpar <- length(get.parms.data(x=Fit$orderedestimations,which='type',type='V'))
  
  nepar   <- nemodpar+nevarpar
  
  if (nepar>=1){
    for (i in 1:nepar){
      lenname <- nchar(Fit$orderedestimations$names[i])
      tabline <- sprintf(paste(c(' %s',rep(' ',maxlen-lenname+2),
                                 '%s   %s  %s     [%s, %s]'),collapse=''),
                         Fit$orderedestimations$names[i],
                         sprintf('%10.4g', Fit$orderedinitial$value[i]),
                         sprintf('%10.4g', Fit$orderedestimations$value[i]),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.3g', Fit$cv[i])),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.4g', Fit$ci[i,1])),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.4g', Fit$ci[i,2])))
      write(tabline,file=files$report,append=TRUE,sep='\n')
    }
  }
  
  # then the block of fixed parameters
    # nb of model parametres
  nfmodpar <- length(get.parms.data(x=Fit$orderedfixed,which='type',type='P'))+
              length(get.parms.data(x=Fit$orderedfixed,which='type',type='L'))+
              length(get.parms.data(x=Fit$orderedfixed,which='type',type='IC'))
    # nb of variance parametres
  nfvarpar <- length(get.parms.data(x=Fit$orderedfixed,which='type',type='V'))
  
  nfpar <- nfmodpar+nfvarpar
  
  if (nfpar>=1){
    for (i in 1:nfpar){
      lenname <- nchar(Fit$orderedfixed$names[i])
      tabline <- sprintf(paste(c(' %s',rep(' ',maxlen-lenname+2),'%s %s'),
                               collapse=''),
                         Fit$orderedfixed$names[i],
                         sprintf('%10.4g',Fit$orderedfixed$value[i]),
                         '       Fixed')
      write(tabline,file=files$report,append=TRUE,sep='\n')
    }
  }
  
  # Adds the covariance matrix
  write(sprintf('\n\n%s\n', 'Covariance matrix'),
        file=files$report,append=TRUE,sep='\n')
  
  if (nepar==1){
    write(' Not available when only one parameter is estimated.\n',
          file=files$report,append=TRUE,sep='\n')
  } else {
    if (any(Fit$cov=='singular')){
      write(' M matrix singular - Covariance matrix unobtainable.\n',
            file=files$report,append=TRUE,sep='\n')
    } else {
      tmpform <- rep(' ',maxlen)
    
      for (i in 1:(nepar-1)){
        tmpmax <- max(10,nchar(Fit$orderedestimations$names[i]))
        tmpform <- paste(c(tmpform,
                           rep(' ',
                               tmpmax-nchar(Fit$orderedestimations$names[i])+2),
                           '%s'),
                         collapse='')
      }
      tmpmax <- max(10,nchar(Fit$orderedestimations$names[nepar]))
      tmpform <- paste(c(tmpform,
                         rep(' ',
                             tmpmax-nchar(Fit$orderedestimations$names[nepar])+2),
                         '%s'),
                       collapse='')
      write(do.call(sprintf,
                    c(list(tmpform),
                      as.list(Fit$orderedestimations$names))),
            file=files$report,append=TRUE,sep='\n')
      
      covmat <- Fit$cov
      
      for (i in 1:nepar){
        tmp <- list()
        lenname <- nchar(Fit$orderedestimations$names[i])
        tmp[[1]] <- paste(c(Fit$orderedestimations$names[i],
                            rep(' ',maxlen-lenname+2)),collapse='')
        
        for (j in 2:(nepar+1)){
          lenname <- nchar(Fit$orderedestimations$names[j-1])
          tmpmax <- max(10,lenname)
          tmpmin <- min(maxlen,lenname)
          # Report only by block, assuming no covariance between model and rv 
          # parameters
          if ((i<=nemodpar & (j-1)<=nemodpar) |
              ((i>nemodpar & (j-1)>nemodpar))){
            if (i<(j-1)){
              tmp[j] <- paste(c(rep(' ',tmpmax-1),'-'),collapse='')
            } else {
              if (abs(covmat[i,j-1])<.Machine$double.eps){
                # Round to 0 if below precision limit
                tmp[j] <- paste(c(rep(' ',tmpmax-1),'0'),collapse='')
              } else {
                tmp[j] <- sprintf('%10.3g', covmat[i,j-1])
                if ((tmpmin-nchar(tmp[j]))>=0){
                  tmp[j] <- paste(c(rep(' ',tmpmin-nchar(tmp[j])),tmp[j]),
                                  collapse='')
                }
              }
            }
            tmp[j] <- paste(c(tmp[j],'  '),collapse='')
          } else {
            tmp[j] <- paste(c(rep(' ',tmpmax-1),'-  '),collapse='')
          }
        }
        write(do.call(sprintf,
                      c(list(paste(rep('%s',nepar+1),collapse='')),
                        as.list(tmp))),
              file=files$report,append=TRUE,sep='\n')
      }
    }
  }
      
  # Adds the correlation matrix
  write(sprintf('\n\n%s\n', 'Correlation matrix'),
        file=files$report,append=TRUE,sep='\n')
  
  if (nepar==1){
    write(' Not available when only one parameter is estimated.\n',
          file=files$report,append=TRUE,sep='\n')
  } else {
    if (any(Fit$cov=='singular')){
      write(' M matrix singular - correlation matrix unobtainable.\n',
            file=files$report,append=TRUE,sep='\n')
    } else { 
      tmpform <- rep(' ',maxlen)
      for (i in 1:(nepar-1)){
        tmpmax <- max(10,nchar(Fit$orderedestimations$names[i]))
        tmpform <- paste(c(tmpform,
                           rep(' ',
                               tmpmax-nchar(Fit$orderedestimations$names[i])+2),
                           '%s'),
                         collapse='')
      }
      tmpmax <- max(10,nchar(Fit$orderedestimations$names[nepar]))
      tmpform <- paste(c(tmpform,
                         rep(' ',
                             tmpmax-nchar(Fit$orderedestimations$names[nepar])+2),
                         '%s'),
                       collapse='')
      write(do.call(sprintf,
                    c(list(tmpform),
                    as.list(Fit$orderedestimations$names))),
            file=files$report,append=TRUE,sep='\n')
          
      cormat <- transpose(Fit$cor)
      
      for (i in 1:nepar){
        tmp <- c()
        lenname <- nchar(Fit$orderedestimations$names[i])
        tmp[1] <- paste(c(Fit$orderedestimations$names[i],
                        rep(' ',maxlen-lenname+2)),collapse='')
        
        for (j in 2:(nepar+1)){
          lenname <- nchar(Fit$orderedestimations$names[j-1])
          tmpmax <- max(10,lenname)
          tmpmin <- min(maxlen,lenname)
          # Report only by block, assuming no covariance between model and rv parameters
          if ((i<=nemodpar & (j-1)<=nemodpar) |
              ((i>nemodpar & (j-1)>nemodpar))){
            if (i<(j-1)){
              tmp[j] <- paste(c(rep(' ',tmpmax-1),'-'),collapse='')
            } else {
              if (abs(cormat[i,j-1])<.Machine$double.eps){
                # Round to 0 if below precision limit
                tmp[j] <- paste(c(rep(' ',tmpmax-1),'0'),collapse='')
              } else {
                tmp[j] <- sprintf('%10.3g', cormat[i,j-1])
                if ((tmpmin-nchar(tmp[j]))>=0){
                  tmp[j] <- paste(c(rep(' ',tmpmin-nchar(tmp[j])),tmp[j]),
                    collapse='')
                }
              }
            }
            tmp[j] <- paste(c(tmp[j],'  '),collapse='')
          } else {
            tmp[j] <- paste(c(rep(' ',tmpmax-1),'-  '),collapse='')
          }
        }
        write(do.call(sprintf,
                      c(list(paste(rep('%s',nepar+1),collapse='')),
                        as.list(tmp))),
              file=files$report,append=TRUE,sep='\n')
      }
    }
  }
  
# Add secondary parameter stats
  if (size(Fit$sec$names,1)!=0){
    # Calculates the dimension of the longest parameter name
    maxlen <- max(9,max(nchar(Fit$sec$names)))
    
    write(sprintf('\n\n%s\n', 'Secondary parameter(s)'),
          file=files$report,append=TRUE,sep='\n')
    
    # Builds and writes the 2 lines of the estimation table title 
    tmp <- sprintf(paste(c(rep(' ',maxlen+6),'%s       %s'),
                         collapse=''),'Initial','Final')
    tmp <- c(tmp,
             sprintf(paste(c(' %s',rep(' ',maxlen-9+6),'%s        ',
                             '%s         ','%s      %s'),collapse=''),
                     'Parameter','value','value','CV%','Confidence interval (95%)'))
    
    write(tmp,file=files$report,append=TRUE,sep='\n')
    
    # Builds and writes each line of the computation table
    for (i in 1:length(Fit$sec$names)){
      lenname <- nchar(Fit$sec$names[i])
      tabline <- sprintf(paste(c(' %s',rep(' ',maxlen-lenname+2),
                                 '%s   %s  %s     [%s, %s]'),collapse=''),
                         Fit$sec$names[i],
                         sprintf('%10.4g', Fit$sec$init[i]),
                         sprintf('%10.4g', Fit$sec$estimates[i]),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.3g', Fit$sec$cv[i])),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.4g', Fit$sec$ci[i,1])),
                         ifelse(any(Fit$cov=='singular'),
                                '        NA',
                                sprintf('%10.4g', Fit$sec$ci[i,2])))
      write(tabline,file=files$report,append=TRUE,sep='\n')    
    }
    
  } else {
    
    write(sprintf('\n\n%s\n', 'No secondary parameter'),
          file=files$report,append=TRUE,sep='\n')
  
  }

# Output parameter estimates to estimate file
  tmp <- problem$init[,c('names','value','isfix')]
  tmp$order <- 1:size(tmp,1)
  estimates <- merge(tmp,Fit$orderedestimations,by='names',all=TRUE,sort=FALSE)
  estimates$value <- ifelse(is.na(estimates$value.y),
                            estimates$value.x,
                            estimates$value.y)
  estimates <- estimates[order(estimates$order),]
  
  # Creates and populates the estimate file if necessary
  if (problem$data$id==1)  
    write(paste(c('ID',estimates$names,Fit$sec$names),collapse=','),
          file=files$est,append=FALSE,sep='\n')
  
  tmpform <- paste(c('%d,',
                     rep('%0.5g,',size(estimates,1)-1+length(Fit$sec$estimates)),
                     '%0.5g'), collapse='')
  estimates <- c(problem$data$id,estimates$value,Fit$sec$estimates)
  
  write(do.call(sprintf,c(list(tmpform),as.list(estimates))),
        file=files$est,append=TRUE,sep='\n')
  
  return(Fit)
  
}

