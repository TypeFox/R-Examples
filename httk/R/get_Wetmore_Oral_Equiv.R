# This functions converts a uM chemical concetrations ("conc") to an oral equivalent dose using the values from the Wetmore et al. (2012) and (2013) publications
get_wetmore_oral_equiv <- function(conc,chem.name=NULL,chem.cas=NULL,suppress.messages=F,which.quantile=0.95,species="Human",input.units='uM',output.units='mg',clearance.assay.conc=NULL,...)
{
  Wetmore.data <- Wetmore.data
  if(is.null(chem.cas)) chem.cas <- get_chem_id(chem.name=chem.name)[['chem.cas']]
  if(tolower(input.units) =='mg/l' | tolower(output.units) == 'mol'){
    MW <- get_physchem_param("MW",chem.CAS=chem.cas)
  }   
  if(tolower(input.units) == 'mg/l'){
    conc <- conc / 1000 / MW * 1000000
  }else if(tolower(input.units)!='um') stop('Input units can only be in mg/L or uM.')
  if(is.null(clearance.assay.conc)){
    this.data <- subset(Wetmore.data,Wetmore.data[,"CAS"]==chem.cas&toupper(Wetmore.data[,"Species"])==toupper(species))
    this.conc <- this.data[which.min(abs(this.data[,'Concentration..uM.'] - conc)),'Concentration..uM.']
  }else{
    this.conc <- clearance.assay.conc
  }
  Css <- try(get_wetmore_css(daily.dose=1,chem.cas=chem.cas,which.quantile=which.quantile,species=species,clearance.assay.conc=this.conc,suppress.messages=T,output.units='uM',...))
  dose <- conc / Css
  if(tolower(output.units) == 'mol'){
    dose <- dose /1000 / MW * 1000000 
  }else if(tolower(output.units) != 'mg') stop("Output units can only be in mg or uM.")
   if (!suppress.messages){
    cat(paste("Retrieving Css from Wetmore et al. based on ",this.conc," uM intrinsic clearance data for the ",which.quantile," quantile in ",species,".\n",sep=""))
    cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),input.units,"concentration converted to",output.units,"/kg bw/day dose.\n")
   }
   if (class(Css) == "try-error"){
     return(NA)
   }else{
# Converting Css for 1 mg/kg/day to dose needed to produce concentration conc uM:   
     return(dose)
   }
}