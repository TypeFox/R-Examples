 calc_analytic_css <- function(chem.name=NULL,chem.cas = NULL,parameters=NULL,daily.dose=1,output.units='uM',model = 'pbtk',species='Human',concentration='plasma',suppress.messages=F,recalc.blood2plasma=F,default.to.human=F)
 {
    if(is.null(chem.cas) & is.null(chem.name) & is.null(parameters)) stop('Must specify chem.cas, chem.name, or parameters.')
    
    dose <- daily.dose
    good.units <- c("uM","mg/L")
    if (!(tolower(output.units) %in% tolower(good.units))) stop(paste("Do not know how to calculate units",output.units,". Please select from: ",paste(good.units,collapse=", ")))
    dose <- dose / 24 
    if(is.null(parameters)){
      if(is.null(chem.cas)){
        out <- get_chem_id(chem.name=chem.name)
        chem.cas <- out$chem.cas
      }
      MW <- get_physchem_param('MW',chem.CAS=chem.cas)
    }else{
      MW <- parameters[['MW']]
    } 
    if(tolower(output.units)=='um'){ 
         dose <- dose / 1000 / MW * 1000000
    }else if(tolower(output.units) != 'mg/l') stop('Output.units can only be uM or mg/L.')
     
    if(tolower(model)=='pbtk')
    {
       if(is.null(parameters)){
         parameters <- parameterize_pbtk(chem.cas=chem.cas,species=species,default.to.human=default.to.human) 
       }else{
         name.list <- c("BW","Clmetabolismc","Funbound.plasma","Fgutabs","Fhep.assay.correction","hematocrit","kdermabs","Kgut2pu","kgutabs","kinhabs","Kkidney2pu","Kliver2pu","Klung2pu","Krbc2pu","Krest2pu","million.cells.per.gliver","MW","Qcardiacc" ,"Qgfrc","Qgutf","Qkidneyf","Qliverf","Rblood2plasma","Vartc","Vgutc","Vkidneyc","Vliverc","Vlungc","Vrestc","Vvenc")
         if(!all(name.list %in% names(parameters)))stop(paste("Missing parameters:",paste(name.list[which(!name.list %in% names(parameters))],collapse=', '),".  Use parameters from parameterize_pbtk."))
       }  
       if(recalc.blood2plasma) parameters[['Rblood2plasma']] <- 1 - parameters[['hematocrit']] + parameters[['hematocrit']] * parameters[['Krbc2pu']] * parameters[['Funbound.plasma']]  


       
       Qcardiac <-  parameters[["Qcardiacc"]] / parameters[['BW']]^0.25  
       Qgfr <-  parameters[["Qgfrc"]] / parameters[['BW']]^0.25    
       Clmetabolism <-  parameters[["Clmetabolismc"]]  
       Kliver2pu <- parameters[['Kliver2pu']]

       Qgut <- parameters[["Qgutf"]] * Qcardiac
       Qliver <- parameters[["Qliverf"]] * Qcardiac
       Qkidney <- parameters[['Qkidneyf']] * Qcardiac
       Qrest <- Qcardiac-Qgut-Qliver-Qkidney
       Rblood2plasma <- parameters[['Rblood2plasma']]
       fub <- parameters[["Funbound.plasma"]]

       dose <- dose * parameters$Fgutabs
       
       Css <- (dose * (Qliver + Qgut) / (fub * Clmetabolism / Rblood2plasma + (Qliver + Qgut))) / (Qcardiac - (Qliver + Qgut)**2 /(fub * Clmetabolism / Rblood2plasma + (Qliver + Qgut)) - Qkidney**2 / (Qgfr * fub / Rblood2plasma + Qkidney) - Qrest)
     if (tolower(concentration)=='plasma')
      {
        Css <- Css / parameters[['Rblood2plasma']]
      } else if (tolower(concentration)!='blood') stop("Only blood and plasma concentrations are calculated.")
    }
    else if (tolower(model)=='3compartmentss')
    {
      if (is.null(parameters)){
        parameters <- parameterize_steadystate(chem.cas=chem.cas,species=species,default.to.human=default.to.human)
      }else{
        name.list <- c("Clint","Funbound.plasma","Fhep.assay.correction","Qtotal.liverc","Qgfrc","BW","MW","Fgutabs","million.cells.per.gliver","Vliverc","liver.density")
        if(!all(name.list %in% names(parameters)))stop(paste("Missing parameters:",paste(name.list[which(!name.list %in% names(parameters))],collapse=', '),".  Use parameters from parameterize_steadystate."))
      }
      if(parameters$Funbound.plasma == 0) stop('Fraction unbound plasma cannot be zero.  Use calc_mc_css or get_wetmore_css to predict steady state for this chemical with three compartment steady state model.')
      
      dose <- dose * parameters$Fgutabs
      Css <- dose/(parameters$Qgfrc/parameters[['BW']]^.25 * parameters$Funbound.plasma + calc_hepatic_clearance(parameters=parameters,suppress.messages=T))
      if (tolower(concentration)=='blood')
      {
        Rb2p <- calc_rblood2plasma(chem.name=chem.name,chem.cas=chem.cas)
        Css <- Css * Rb2p
      } else if (tolower(concentration)!='plasma') stop("Only blood and plasma concentrations are calculated.")      
    }else if(tolower(model) == '1compartment'){
      if(is.null(parameters)){
        parameters <- parameterize_1comp(chem.cas=chem.cas,species=species,default.to.human=default.to.human)
      }else{
        name.list <- c("Vdist","million.cells.per.gliver","kelim","kgutabs","Rblood2plasma","MW","hematocrit","Fgutabs")
        if(!all(name.list %in% names(parameters)))stop(paste("Missing parameters:",paste(name.list[which(!name.list %in% names(parameters))],collapse=', '),".  Use parameters from parameterize_1comp."))
      }
      dose <- dose * parameters$Fgutabs
      Css <- dose / parameters$kelim / parameters$Vdist
      if (tolower(concentration)=='blood')
      {
        Css <- Css * parameters[['Rblood2plasma']]
      } else if (tolower(concentration)!='plasma') stop("Only blood and plasma concentrations are calculated.")
    }else if(tolower(model) == '3compartment'){
      if (is.null(parameters)){
        parameters <- parameterize_3comp(chem.cas=chem.cas,species=species,default.to.human=default.to.human)
      }else{
        name.list <- c("BW","Clmetabolismc","Funbound.plasma","Fgutabs","Fhep.assay.correction","hematocrit","Kgut2pu","Krbc2pu","kgutabs","Kliver2pu","Krest2pu","million.cells.per.gliver","MW","Qcardiacc","Qgfrc","Qgutf","Qliverf","Rblood2plasma","Vgutc","Vliverc","Vrestc")
        if(!all(name.list %in% names(parameters)))stop(paste("Missing parameters:",paste(name.list[which(!name.list %in% names(parameters))],collapse=', '),".  Use parameters from parameterize_3comp."))
        name.list2 <- c("BW","Clmetabolismc","Fgutabs","Funbound.plasma","Fhep.assay.correction","hematocrit","kdermabs","Kgut2pu","kgutabs","kinhabs","Kkidney2pu","Kliver2pu","Klung2pu","Krbc2pu","Krest2pu","million.cells.per.gliver","MW","Qcardiacc" ,"Qgfrc","Qgutf","Qkidneyf","Qliverf","Rblood2plasma","Vartc","Vgutc","Vkidneyc","Vliverc","Vlungc","Vrestc","Vvenc")
        if(any(name.list2[which(!name.list2 %in% name.list)] %in% names(parameters)))stop("Parameters are from parameterize_pbtk.  Use parameters from parameterize_3comp.")
      }
      if(recalc.blood2plasma) parameters[['Rblood2plasma']] <- 1 - parameters[['hematocrit']] + parameters[['hematocrit']] * parameters[['Krbc2pu']] * parameters[['Funbound.plasma']]
      

      dose <- dose * parameters$Fgutabs
      Css <- dose * parameters[['BW']]^0.25  / (parameters$Clmetabolismc * parameters[['BW']]^0.25  + parameters$Qgfrc * (parameters$Qliverf + parameters$Qgutf) * parameters$Qcardiacc / ((parameters$Qliverf + parameters$Qgutf) * parameters$Qcardiacc + parameters$Funbound.plasma * parameters$Qgfrc / parameters$Rblood2plasma)) / parameters$Funbound.plasma
      if (tolower(concentration)=='blood')
      {
         Css <- Css * parameters[['Rblood2plasma']]
      } else if (tolower(concentration)!='plasma') stop("Only blood and plasma concentrations are calculated.")
    }else stop("Model must be either \"pbtk\", \"1compartment\", \"3compartmentss\", or \"3compartment\".")
    
  if(!suppress.messages){
    if(is.null(chem.cas) & is.null(chem.name)){
      cat(paste(toupper(substr(concentration,1,1)),substr(concentration,2,nchar(concentration)),sep=''),"concentrations returned in",output.units,"units.\n")
    }else cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),concentration,"concentrations returned in",output.units,"units.\n")
  }    
       

 #  Css <- dose * Rblood2plasma * (Qliver + Qgut) * (Qkidney * Rblood2plasma  + Qgfr * fub)  / ((##(Rblood2plasma * Qkidney + Qgfr * fub) * (Qcardiac - Qrest) - Qkidney**2 * Rblood2plasma) * ((Qliver + Qgut) * Rblood2plasma +  Clmetabolism * fub * Kliver2plasma) - (Qliver + Qgut)**2 * (Qkidney *Rblood2plasma + Qgfr * fub)*Rblood2plasma)
  
    #Css <- dose/(QGFRc*fub+calc_Hepatic_Clearance(Params))
    return(as.numeric(Css))
}