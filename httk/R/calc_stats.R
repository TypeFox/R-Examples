calc_stats <-function(days,chem.name=NULL,chem.cas=NULL,parameters=NULL,stats=c("AUC","peak","mean"),species='Human',exclude.fub.zero=F,daily.dose=1,dose=NULL,doses.per.day=NULL,output.units='uM',concentration='plasma',model='pbtk',suppress.messages=F,...)
{
  AUC <- NULL
  peak <- NULL
  mean <- NULL
  out <- NULL
  
  if(is.null(chem.name) & is.null(chem.cas) & is.null(parameters)){
    for(this.CAS in get_cheminfo(species=species,exclude.fub.zero = exclude.fub.zero,model=model)){
      stat <- calc_chem_stats(chem.cas=this.CAS,days=days,stats=stats,species=species,dose=dose,daily.dose=daily.dose,doses.per.day=doses.per.day,concentration=concentration,output.units=output.units,model=model,suppress.messages=T,...)

      if(length(stat)==1){
        out[this.CAS] <-  stat 
      }else{
        AUC[this.CAS] <- stat[["AUC"]]
        peak[this.CAS] <- stat[["peak"]] 
        mean[this.CAS] <- stat[["mean"]] 
      }
    }
    if(length(stat)!=1){
      if(!is.null(AUC) & !is.null(peak) & is.null(mean)) out <- list(AUC=AUC,peak=peak)
      else if(!is.null(AUC) & is.null(peak) & !is.null(mean)) out <- list(AUC=AUC,mean=mean)
      else if(is.null(AUC) & !is.null(peak) & !is.null(mean)) out <- list(mean=mean,peak=peak)
      else out <- list(AUC=AUC,peak=peak,mean=mean)
    }
    if(!suppress.messages){
      cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),concentration,"concentrations returned in",output.units,"units.\n")
      if('auc' %in% tolower(stats)) cat("AUC is area under plasma concentration curve in",output.units,"* days units.\n")
    }
  }else{
    out <- calc_chem_stats(chem.name=chem.name,chem.cas=chem.cas,parameters=parameters,days=days,stats=stats,species=species,daily.dose=daily.dose,dose=dose,doses.per.day=doses.per.day,concentration=concentration,output.units=output.units,model=model,suppress.messages=suppress.messages,...)
  } 
  return(out)
}