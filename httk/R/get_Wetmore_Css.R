# This function retrives a steady-state chemical concentration from the Wetmore et al. (2012) and (2013) publications
get_wetmore_css <- function(chem.cas=NULL,chem.name=NULL,daily.dose=1,which.quantile=0.95,species="Human",clearance.assay.conc=NULL,output.units="mg/L",suppress.messages=F)
{
  Wetmore.data <- Wetmore.data
  if (species == "Human") available.quantiles <- c(0.05,0.5, 0.95)
  else available.quantiles <- 0.5
  if (!all(which.quantile %in% available.quantiles)) stop("Wetmore papers only includes 5%, 50%, and 95% quantiles for human and 50% for rat.")
  if (!(tolower(output.units) %in% c("mg/l","um"))) stop("Wetmore papers only includes mg/L and uM values for Css")
  out <- get_chem_id(chem.cas=chem.cas,chem.name=chem.name)
  chem.cas <- out$chem.cas
  chem.name <- out$chem.name
    
  this.data <- subset(Wetmore.data,Wetmore.data[,"CAS"]==chem.cas&toupper(Wetmore.data[,"Species"])==toupper(species))
   
    if (!is.null(clearance.assay.conc)) 
    {
      this.data <- subset(this.data,this.data[,"Concentration..uM."]==clearance.assay.conc)[1,]
      if (dim(this.data)[1]!=1) stop(paste("No",clearance.assay.conc,"uM clearance assay data for",chem.name,"in",species))
    }else{
      if(1 %in% this.data[,"Concentration..uM."]){
        this.data <- subset(this.data,this.data[,"Concentration..uM."]== 1)[1,] 
        clearance.assay.conc <- 1
      }else{
        this.data <- this.data[1,]
        clearance.assay.conc <- this.data[,"Concentration..uM."][[1]]
      }
    }
    out <- NULL
    if (tolower(output.units)=="mg/l")
    {
      if (0.05 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_lower_5th_perc.mg.L."])
      if (0.5 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_median_perc.mg.L."])
      if (0.95 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_upper_95th_perc.mg.L."])
    } else if(tolower(output.units) == 'um') {
      if (0.05 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_lower_5th_perc.uM."])
      if (0.5 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_median_perc.uM."])
      if (0.95 %in% which.quantile) out <- c(out,daily.dose*this.data[,"Css_upper_95th_perc.uM."])
    } else{
     stop('Output.units can only be uM or mg/L.')
    }
  
  if(!suppress.messages){
    cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),"plasma concentrations returned in",output.units,"units.\n")
    cat("Retrieving Css from Wetmore et al. based on ",clearance.assay.conc," uM intrinsic clearance data for the ",which.quantile," quantile in ",species,".\n")
  }
  return(out)
}