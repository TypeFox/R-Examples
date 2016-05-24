calc_css <- function(parameters=NULL,
                    chem.name=NULL,
                    chem.cas=NULL, 
                    species='Human',
                    f = .01,
                    daily.dose=1,
                    doses.per.day=3,
                    days = 10,
                    output.units = "uM",
                    concentration='plasma',
                    suppress.messages=F,
                    model='pbtk',
                    default.to.human=F,
                    f.change = 0.00001,
                    ...)
{
  
  if(is.null(parameters)){
    if(tolower(model)=='pbtk'){
      parameters <- parameterize_pbtk(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
    }else if(tolower(model)=='3compartment'){
      parameters <- parameterize_3comp(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
    }else if(tolower(model)=='1compartment'){
      parameters <- parameterize_1comp(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
    }
  } 

  css <- calc_analytic_css(parameters=parameters,daily.dose=daily.dose,concentration='plasma',model=model,suppress.messages=T) 
  conc <- (1 - f) * css 

  if(tolower(model) == 'pbtk'){
    out <- solve_pbtk(parameters=parameters, daily.dose=daily.dose,doses.per.day=doses.per.day,days = days,suppress.messages=T,...)
    Final_Conc <- out[dim(out)[1],c("Agutlumen","Cart","Cven","Clung","Cgut","Cliver","Ckidney","Crest")]
  }else if(tolower(model) =='3compartment'){
    out <- solve_3comp(parameters=parameters, daily.dose=daily.dose,doses.per.day=doses.per.day, days = days,suppress.messages=T,...)
    Final_Conc <- out[dim(out)[1],c("Agutlumen","Cgut","Cliver","Crest")]
  }else if(tolower(model)=='1compartment'){
    out <- solve_1comp(parameters=parameters,daily.dose=daily.dose,doses.per.day=doses.per.day, days = days,suppress.messages=T,...)
    Final_Conc <- out[dim(out)[1],c("Agutlumen","Ccompartment")]
  }else stop('The model options are only: 1compartment, 3compartment, and pbtk.')
  
  day <- days
  past.out <- NULL
  while(out[dim(out)[1],'AUC'] - out[match(days - 1,out[,'time']),'AUC'] < conc & (out[dim(out)[1],'AUC'] - out[match(days - 1,out[,'time']),'AUC'])/
       (out[match(days - 1,out[,'time']),'AUC'] - out[match(days - 2,out[,'time']),'AUC']) > 1 + f.change)
  {
    if(day < 3600)
    {
      days <- days * 3
    }else{
      days <- days * 2
    }
    day <- day + days
    
  if(tolower(model) == 'pbtk'){
    past.out <- as.data.frame(out)
    out <- solve_pbtk(parameters=parameters,initial.values = Final_Conc, daily.dose=daily.dose,doses.per.day=doses.per.day, days = days,suppress.messages=T,...)
    Final_Conc <- out[dim(out)[1],c("Agutlumen","Cart","Cven","Clung","Cgut","Cliver","Ckidney","Crest")]
  }else if(tolower(model) =='3compartment'){
    past.out <- as.data.frame(out)
    out <- solve_3comp(parameters=parameters,initial.values = Final_Conc, daily.dose=daily.dose,doses.per.day=doses.per.day, days = days,suppress.messages=T,...)
    Final_Conc <- out[dim(out)[1],c("Agutlumen","Cgut","Cliver","Crest")]
  }else if(tolower(model)=='1compartment'){
    past.out <- as.data.frame(out)
    out <- solve_1comp(parameters=parameters,daily.dose=daily.dose,doses.per.day=doses.per.day, days = days,suppress.messages=T,initial.values=Final_Conc,...)
    Final_Conc <- out[dim(out)[1],c('Agutlumen','Ccompartment')]
  }
  
    if(day > 36500) break 
  }
  
  if(out[dim(out)[1],'AUC'] - out[match(days - 1,out[,'time']),'AUC'] > conc)
  {
    dif <- out[match(1,out[,'time']):dim(out)[1],'AUC'] - out[1:match(days - 1,out[,'time']),'AUC']   #area under curve for 24 hours previous to each point.  length of days-1
    if(dif[1] > conc)the.day <- day - days + 1              #first day not included in dif
    else the.day <- day - days + 1 + ceiling(min(which(dif > conc)) / length(dif) * (days - 1))   #add one since first day is not included.  then add the day where we see the first point such that the last 24 hours is above the concentration
  }else{
    if((out[dim(out)[1],'AUC'] - out[match(days - 1,out[,'time']),'AUC'])/
      (out[match(days - 1,out[,'time']),'AUC'] - out[match(days - 2,out[,'time']),'AUC']) < 1 + f.change){
      out[,'AUC'] <- out[,'AUC'] + past.out[dim(past.out)[1],'AUC']           #make AUC consistent with the previous output
      end.past.out <- subset(past.out,time >= past.out[[dim(past.out)[1],'time']] - 1)   #only take the last day
      end.past.out[,'time'] <- end.past.out[,'time'] - end.past.out[dim(end.past.out)[1],'time']       #make the time negative
      complete.out <- rbind(end.past.out[1:(dim(end.past.out)[1]-1),],out)                             #combine them
      #same as dif but divided by the area under the curve of the 24 hours prior to each point.  past.out is added to account for the second day
      div.dif <- (complete.out[match(1,complete.out[,'time']):dim(complete.out)[1],'AUC'] - complete.out[match(0,complete.out[,'time']):match(days - 1,complete.out[,'time']),'AUC']) / 
               (complete.out[match(0,complete.out[,'time']):match(days - 1,complete.out[,'time']),'AUC'] - complete.out[1:match(days - 2,complete.out[,'time']),'AUC'])
      if(div.dif[1] < 1 + f.change)the.day <- day - days + 1 
      else the.day <- day - days + 1 + ceiling(min(which(div.dif < 1 + f.change)) / length(div.dif) * (days - 1))
    }else{ 
     if(!suppress.messages)cat("Steady state not reached after 100 years.")
     the.day <- 36500   
    }     
  }
  
  if (tolower(output.units) == tolower("mg/L")) 
  {
      out[,'AUC'] <- out[,'AUC']/1e+06 * parameters[["MW"]] * 1000
      css <- css /1e+06 * parameters[["MW"]] * 1000
      if(tolower(model)=='1compartment'){
        out[,'Ccompartment'] <- out[,'Ccompartment']/1e+06 * parameters[["MW"]] * 1000
      }else{  
        out[,'Cplasma'] <- out[,'Cplasma']/1e+06 * parameters[["MW"]] * 1000
      }
  } else if (tolower(output.units) != tolower("uM")) stop("Currently can only return units of mg/L and uM")
  
  if(tolower(concentration)=='plasma'){
    if(tolower(model)=='1compartment'){
      max=as.numeric(max(out[,'Ccompartment']))
    }else{
      max=as.numeric(max(out[,'Cplasma']))
    }
    avg=as.numeric(out[dim(out)[1],'AUC'] - out[match(days-1,out[,'time']),'AUC'])
  }else if(tolower(concentration)=='blood'){
    if(tolower(model)=='pbtk'){
      max=as.numeric(max(out[,'Cven']))
    }else if(tolower(model) == '3compartment'){
      max=as.numeric(max(out[,'Cplasma'] * parameters[['Rblood2plasma']]))
    }else{
      max=as.numeric(max(out[,'Ccompartment'] * parameters[['Rblood2plasma']]))
    }   
   avg=as.numeric((out[dim(out)[1],'AUC'] - out[match(days-1,out[,'time']),'AUC'])*parameters[['Rblood2plasma']])
  }else stop("Only blood and plasma concentrations are calculated.")
  if(!suppress.messages){
    if(is.null(chem.cas) & is.null(chem.name)){
      cat(paste(toupper(substr(concentration,1,1)),substr(concentration,2,nchar(concentration)),sep=''),"concentrations returned in",output.units,"units.\n")
    }else cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),concentration,"concentrations returned in",output.units,"units.\n")
  }
  return(list(avg=avg,frac=as.numeric((out[match(days - (day - the.day),out[,'time']),'AUC'] - out[match(days - 1 - (day - the.day),out[,'time']),'AUC']) / css), 
    max=max,
    the.day =as.numeric(the.day)))
}