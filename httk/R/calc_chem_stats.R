# This function simulates a daily repeated dosing treatment and summarizes the PK statistics from that treatment:
calc_chem_stats <- function(chem.name=NULL,chem.cas=NULL,parameters=NULL,days,stats=c("AUC","mean","peak"),doses.per.day=NULL,daily.dose=1,dose=NULL,species="Human",output.units="uM",concentration='plasma',model='pbtk',suppress.messages=F,...)
{        
  good.units <- c("uM","mg/L")
  if (!(tolower(output.units) %in% tolower(good.units))) stop(paste("Do not know how to calculate units",output.units,". Please select from: ",paste(good.units,collapse=" ")))
  valid.stats <- c("AUC","mean","peak")
  
  if(is.null(parameters)){
    if(tolower(model) == 'pbtk'){
      parameters <- parameterize_pbtk(chem.name=chem.name,chem.cas=chem.cas,species=species)
    }else if(tolower(model) == '3compartment'){
      parameters <- parameterize_3comp(chem.name=chem.name,chem.cas=chem.cas,species=species)
    }else if(tolower(model) == '1compartment'){
      parameters <- parameterize_1comp(chem.name=chem.name,chem.cas=chem.cas,species=species)
    }else stop("Model can only be '3compartment', '1compartment', or 'pbtk'.")
  }
  
  if(tolower(model) == 'pbtk'){
    PKtimecourse <- solve_pbtk(parameters=parameters,days = days,species=species,dose=dose,doses.per.day=doses.per.day,suppress.messages=T,output.units=output.units,daily.dose=daily.dose,...)
  }else if(tolower(model) == '3compartment'){
    PKtimecourse <- solve_3comp(parameters=parameters,days = days,species=species,dose=dose,doses.per.day=doses.per.day,suppress.messages=T,output.units=output.units,daily.dose=daily.dose,...)
  }else if(tolower(model) == '1compartment'){
    PKtimecourse <- solve_1comp(parameters=parameters,days = days,species=species,dose=dose,doses.per.day=doses.per.day,suppress.messages=T,output.units=output.units,daily.dose=daily.dose,...)
    colnames(PKtimecourse)[which(colnames(PKtimecourse) == 'Ccompartment')] <- 'Cplasma'
  } 
  
  #if(tolower(output.units) == tolower("mg/L")) PKtimecourse[,c('AUC','Cserum')] <- PKtimecourse[,c('AUC','Cserum')] * parameters[['MW']] / 1000 
  output <- list()
  
  # If mean is requested, calculate it last in case AUC is also requested:
  if ("mean" %in% stats) stats <- c(stats[stats!="mean"],"mean")
  for (this.stat in stats)
  {
    if (!(tolower(this.stat) %in% tolower(valid.stats))) stop(paste("calc_stats cannot calculate",this.stat,". Valid stats are:",paste(valid.stats,collapse=" "),"."))
    if (tolower(this.stat) == "auc") output[["AUC"]] <- as.numeric(PKtimecourse[dim(PKtimecourse)[1],'AUC'])
    if (tolower(this.stat) == "peak") output[["peak"]] <- calc_timecourse_peak(PKtimecourse[,c("time",'Cplasma')])
    if (tolower(this.stat) == "mean")
    {
      if(!is.null(output[["AUC"]])) output[["mean"]] <- output[["AUC"]]/days
      else output[["mean"]] <- as.numeric(PKtimecourse[dim(PKtimecourse)[1],'AUC']/days) 
    }
  }

  # If only one stat was asked for, don't return a list, return just the first entry in the list:
  if (length(output) == 1) output <- output[[1]]
  if(tolower(concentration)=='blood'){
    if(length(output) == 1){
      output <- output * parameters[['Rblood2plasma']]
    }else{
      for(this.stat in stats){
        output[[this.stat]] <- output[[this.stat]] *  parameters[['Rblood2plasma']]
      }
    }
  }else if(tolower(concentration) != 'plasma') stop("Only blood and plasma concentrations are calculated.")

  if(!suppress.messages){
    if(is.null(chem.cas) & is.null(chem.name)){
      cat(paste(toupper(substr(concentration,1,1)),substr(concentration,2,nchar(concentration)),sep=''),"values returned in",output.units,"units.\n")
    }else cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),concentration,"concentrations returned in",output.units,"units.\n")
    
    if('AUC' %in% stats) cat("AUC is area under plasma concentration curve in",output.units,"* days units with Rblood2plasma =",parameters[['Rblood2plasma']],".\n")
  }    
  
  return(output)
}