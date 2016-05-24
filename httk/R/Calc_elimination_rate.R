# This function calculats an elimination rate for a one compartment model where 
# eelimination is entirely due to metablism by the liver and glomerular filtration
# in the kidneys.
calc_elimination_rate <- function(chem.cas=NULL,chem.name=NULL,parameters=NULL,species="Human",suppress.messages=F,
                                 default.to.human=F)
{
    name.list <- c("Clint","Funbound.plasma","Qtotal.liverc","Qgfrc","BW","million.cells.per.gliver","Vliverc","liver.density",'Fhep.assay.correction')
    if(is.null(parameters)){
      parameters <- parameterize_steadystate(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
      Vd <- calc_vdist(chem.cas=chem.cas,chem.name=chem.name,species=species,suppress.messages=T,default.to.human=default.to.human) 
    }else{ 
      if(!all(name.list %in% names(parameters))){
        if(is.null(chem.cas) & is.null(chem.name))stop('chem.cas or chem.name must be specified when not including all 3compartmentss parameters.')
        params <- parameterize_steadystate(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
        parameters <- c(parameters,params[name.list[!(name.list %in% names(parameters))]])
      }
      if('Vdist' %in% names(parameters)){
        Vd <- parameters[['Vdist']]
      }else{
        if(is.null(chem.name) & is.null(chem.cas))stop('chem.cas or chem.name must be specified when Vdist is not included in parameters.')
        Vd <- calc_vdist(parameters=parameters,chem.cas=chem.cas,chem.name=chem.name,species=species,suppress.messages=T,default.to.human=default.to.human) 
      }    
    } 
    clearance <- calc_total_clearance(parameters=parameters,suppress.messages=T,default.to.human=default.to.human) #L/h/kgBW

    if(!suppress.messages)cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),"elimination rate returned in units of 1/h.\n")
    return(as.numeric(clearance/Vd))
}
