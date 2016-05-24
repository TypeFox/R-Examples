solve_1comp <- function(chem.cas=NULL,
                        chem.name=NULL,
                        times=NULL,
                        parameters=NULL,
                        daily.dose=1,
                        dose = NULL,
                        doses.per.day=NULL,
                        days=10,
                        tsteps=4,
                        suppress.messages=F,
                        species='Human',
                        output.units='uM',
                        plots=F,
                        initial.values=NULL,
                        iv.dose = F,
                        method="lsoda",rtol=1e-8,atol=1e-12,
                        default.to.human=F,
                        dosing.matrix = NULL,
                        recalc.elimination=F,                        
                        ...)
{     
   Agutlumen <- Acompartment <- Ccompartment <- NULL 
  if(is.null(chem.cas) & is.null(chem.name) & is.null(parameters)) stop('Parameters, chem.name, or chem.cas must be specified.')
  if(is.null(parameters)){  parameters <- parameterize_1comp(chem.name=chem.name,chem.cas=chem.cas,species=species,default.to.human=default.to.human) 
  }else{
    name.list <- c("Vdist","kelim","kgutabs","Rblood2plasma","MW","million.cells.per.gliver","hematocrit","Fgutabs")
    if(!all(name.list %in% names(parameters)))stop(paste("Missing parameters:",paste(name.list[which(!name.list %in% names(parameters))],collapse=', '),".  Use parameters from parameterize_1comp."))
  }
  Rb2p <- parameters[['Rblood2plasma']]
    


  if(is.null(times)) times <- round(seq(0, days, 1/(24*tsteps)),8)
  start <- times[1]
  end <- times[length(times)] 
  if(iv.dose){
    doses.per.day <- NULL
    dosing.matrix <- NULL
    if(is.null(dose)) dose <- daily.dose
  }else{   
    if(is.null(dosing.matrix)){ 
      if(is.null(dose)){
        if(!is.null(doses.per.day)){
          dose <- daily.dose / doses.per.day * parameters$Fgutabs
        }else dose <- daily.dose * parameters$Fgutabs
      }else dose <- dose * parameters$Fgutabs
    }else{
      if(!is.null(dim(dosing.matrix))){
        rc <- which(dim(dosing.matrix) == 2)
        if(rc == 1){
          if(is.null(rownames(dosing.matrix)) | any(!(rownames(dosing.matrix) %in% c('time','dose')))) stop('dosing.matrix must have column or row names of \"time\" and \"dose\" or be a vector of times.')
          dosing.times <- dosing.matrix['time',]
          dose.vector <- dosing.matrix['dose',] * parameters$Fgutabs
        }else{
          if(is.null(colnames(dosing.matrix)) | any(!(colnames(dosing.matrix) %in% c('time','dose')))) stop('dosing.matrix must have column or row names of \"time\" and \"dose\" or be a vector of times.')
          dosing.times <- dosing.matrix[,'time']
          dose.vector <- dosing.matrix[,'dose'] * parameters$Fgutabs   
        }
      }else{
        if(is.null(dose)) stop("Argument dose must be entered to overwrite daily.dose when a time vector is entered into dosing.matrix.")
        dosing.times <- dosing.matrix
        dose.vector <- rep(dose * parameters$Fgutabs,length(dosing.matrix))
      } 
      if(start == dosing.times[1]){
        dose <- dose.vector[[1]]
        dosing.times <- dosing.times[2:length(dosing.times)]
        dose.vector <- dose.vector[2:length(dose.vector)]
      }else dose <- 0  
    } 
  }
   
  if(tolower(output.units)=='um' |  tolower(output.units) == 'mg/l') use.amounts <- F
  if(tolower(output.units)=='umol' |  tolower(output.units) == 'mg') use.amounts <- T 
   
  if(tolower(output.units)=='um' | tolower(output.units) == 'umol')
  {
    dose <- as.numeric(dose / 1000 / parameters[["MW"]] * 1000000)
    if(!is.null(dosing.matrix)) dose.vector <- as.numeric(dose.vector / 1000 / parameters[["MW"]] * 1000000)
  } else if(!(tolower(output.units) == 'mg/l' | tolower(output.units) == 'mg'))
  {
    stop('Output.units can only be uM, umol, mg, or mg/L.')

  }
  
   if (use.amounts)
  {
    CompartmentsToInitialize <-c("Agutlumen","Acompartment")
  } else {
    CompartmentsToInitialize <-c("Agutlumen","Ccompartment")
  }

  for (this.compartment in CompartmentsToInitialize)
  {
  # If the compartment has a value specified in the list initial.values, then set it to that value:
    if (this.compartment %in% names(initial.values))
    {
      eval(parse(text=paste(this.compartment,"<-",initial.values[[this.compartment]])))
      
    }
  # Otherwise set the value to zero:
    else eval(parse(text=paste(this.compartment,"<- 0")))
  }

  if(use.amounts)
  {
    if(iv.dose)
    {
      state <- c(Agutlumen=Agutlumen,Acompartment = Acompartment + dose ,Ametabolized = 0,AUC=0)
    }else{
      state <- c(Agutlumen= Agutlumen + dose,Acompartment = Acompartment ,Ametabolized = 0,AUC=0)
    }
  }else{
    if(iv.dose)
    {
      state <- c(Agutlumen=Agutlumen,Acompartment = dose + Ccompartment * parameters[['Vdist']],Ametabolized = 0,AUC=0)
    } else{
      state <- c(Agutlumen= Agutlumen + dose,Acompartment = Ccompartment * parameters[['Vdist']],Ametabolized = 0,AUC=0)
    }
  }
  
  if(recalc.elimination){
  parameters$kelim <- calc_elimination_rate(parameters=parameters,chem.cas=chem.cas,chem.name=chem.name,species=species,suppress.messages=T,default.to.human=default.to.human)

  }
  parameters[['ke']] <- parameters[['kelim']]   
  parameters[['vdist']] <- parameters[['Vdist']]
  parameters <- initparms1comp(parameters[!(names(parameters) %in% c("kelim","Rblood2plasma","MW",'Vdist','million.cells.per.gliver','hematocrit','Fgutabs'))])
  
  state <-initState1comp(parameters,state)
  

  if(is.null(dosing.matrix)){
    if(is.null(doses.per.day)){
      out <- ode(y = state, times = times,func="derivs1comp", parms=parameters, method=method,rtol=rtol,atol=atol, dllname="httk",initfunc="initmod1comp", nout=length(Outputs1comp),outnames=Outputs1comp,...)
    }else{
      dosing <- seq(start + 1/doses.per.day,end-1/doses.per.day,1/doses.per.day)
      length <- length(dosing)
      eventdata <- data.frame(var=rep('Agutlumen',length),time = round(dosing,8),value = rep(dose,length), method = rep("add",length))                          
      out <- ode(y = state, times = times, func="derivs1comp", parms = parameters, method=method,rtol=rtol,atol=atol, dllname="httk",initfunc="initmod1comp", nout=length(Outputs1comp),outnames=Outputs1comp,events=list(data=eventdata),...)
    }
  }else{
    eventdata <- data.frame(var=rep('Agutlumen',length(dosing.times)),time = dosing.times,value = dose.vector, method = rep("add",length(dosing.times)))                          
    out <- ode(y = state, times = times, func="derivs1comp", parms = parameters, method=method,rtol=rtol,atol=atol, dllname="httk",initfunc="initmod1comp", nout=length(Outputs1comp),outnames=Outputs1comp,events=list(data=eventdata),...)
  }
   


 if(plots==T)
  {
    plot(out,select=c(CompartmentsToInitialize,"Ametabolized","AUC"))
  }
  
  out <- out[,c("time",CompartmentsToInitialize,"Ametabolized","AUC")]
  class(out) <- c('matrix','deSolve')
  
  if(!suppress.messages){
    if(is.null(chem.cas) & is.null(chem.name)){
      if(use.amounts){
        cat("Values returned in",output.units,"/ kg BW units.\n")
      }else{
        if(tolower(output.units) == 'um'){
          out.amount <- 'umol'
        }else out.amount <- 'mg'
        cat("Amounts returned in",out.amount,"/ kg BW and concentration returned in",output.units,"units.\n")
      }
    }else{
      if(use.amounts){
        
        cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),"values returned in",output.units,"/ kg BW units.\n")
      }else{
        if(tolower(output.units) == 'um'){
          out.amount <- 'umol'
        }else out.amount <- 'mg'
        cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),"amounts returned in",out.amount,"/ kg BW and concentration returned in",output.units,"units.\n")
      }
    }
    if(tolower(output.units) == 'mg'){
      cat("AUC is area under compartment concentration in mg/L * days units with Rblood2plasma =",Rb2p,".\n")
    }else if(tolower(output.units) == 'umol'){
      cat("AUC is area under compartment concentration in uM * days units with Rblood2plasma =",Rb2p,".\n")
    }else cat("AUC is area under plasma concentration curve in",output.units,"* days units with Rblood2plasma =",Rb2p,".\n")
  }

  return(out)
}