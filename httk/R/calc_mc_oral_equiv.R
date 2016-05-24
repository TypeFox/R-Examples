calc_mc_oral_equiv <- function(conc,chem.name=NULL,chem.cas=NULL,which.quantile=0.95,species="Human",input.units='uM',output.units='mg',suppress.messages=F,return.samples=F,...)
{
  if(!(tolower(input.units) %in% c('um','mg/l'))) stop("Input units can only be uM or mg/L.")
  Css <- try(calc_mc_css(daily.dose=1,chem.name=chem.name,chem.cas=chem.cas,which.quantile=which.quantile,species=species,output.units=input.units,suppress.messages=T,return.samples=return.samples,...))
  dose <- conc/Css
  if(tolower(output.units) == 'mol'){
    if(is.null(chem.cas)) chem.cas <- get_chem_id(chem.name=chem.name)[['chem.cas']]
    MW <- get_physchem_param("MW",chem.CAS=chem.cas)
    dose <- dose /1000 / MW * 1000000 
  }else if(tolower(output.units) != 'mg') stop("Output units can only be in mg or mol.")
  if(!suppress.messages & !return.samples){
    cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),input.units,"concentration converted to",output.units,"/kg bw/day dose for",which.quantile,"quantile.\n")
  }
	if (class(Css) == "try-error"){
    return(NA)
  }else{
    return(dose)
  }
  
}