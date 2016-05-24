# from Ito and Houston (2004)
calc_hepatic_clearance <- function(chem.name=NULL,chem.cas=NULL,parameters=NULL,species='Human',hepatic.model='well-stirred',suppress.messages=F)
{
  model <- hepatic.model
  name.list <- c("Clint","Funbound.plasma","Qtotal.liverc","million.cells.per.gliver","Vliverc","BW","liver.density",'Fhep.assay.correction')
  if(is.null(parameters)){
  parameters <- parameterize_steadystate(chem.cas=chem.cas,chem.name=chem.name,species=species)
  }else if(!all(name.list %in% names(parameters))){
    if(is.null(chem.cas) & is.null(chem.name))stop('chem.cas or chem.name must be specified when not including all necessary 3compartmentss parameters.')
    params <- parameterize_steadystate(chem.cas=chem.cas,chem.name=chem.name,species=species)
    parameters <- c(parameters,params[name.list[!(name.list %in% names(parameters))]])
  }
  Clint <- get_param("Clint",parameters,"calc_Hepatic_Clearance") # uL/min/10^6 cells

    fu_hep <- get_param("Fhep.assay.correction",parameters,"calc_Hepatic_Clearance") 
   #try(get_param("Fraction_unbound_hepatocyteassay",parameters,"calc_Hepatic_Clearance")) # fraction set if paramaterize function called with fu_hep_correct=TRUE
  #if (class(fu_hep) == "try-error") fu_hep <- 1
# Correct for fraction of chemical unbound in in vitro hepatocyte assay:
  Clint <- Clint / fu_hep

  fub <- get_param("Funbound.plasma",parameters,"calc_Hepatic_Clearance") # unitless fraction
  Qtotal.liverc <- get_param("Qtotal.liverc",parameters,"calc_Hepatic_Clearance",default=1.24) # L/h/kgBW
  Vliverc <- get_param("Vliverc",parameters,"calc_Hepatic_Clearance") #  L/kg BW
  liver.density <- get_param("liver.density",parameters,"calc_Hepatic_Clearance") # g/mL
  Dn <- get_param("Dn",parameters,"calc_Hepatic_Clearance",default=0.17) #
  #model <- get_param("model",parameters,"calc_Hepatic_Clearance",default="well-stirred")
  million.cells.per.gliver <- get_param("million.cells.per.gliver",parameters,"calc_Hepatic_Clearance") # 10^6 cells/g-liver
  
  if (!(tolower(model) %in% c("well-stirred","parallel tube","dispersion","unscaled")))
    stop("Model other than \"well-stirred,\" \"parallel tube,\", \"dispersion\", or \"unscaled\" specified.")

  # Convert from uL/min/10^6 cells to uL/min/g-liver to uL/min/kg BW
  Clint <- Clint*million.cells.per.gliver
  # Convert from uL/min/g-liver to uL/min/kg BW
  Clint <- Clint*(Vliverc*1000*liver.density)
  # Convert from uL/min/kg BW to L/h/kg BW
  Clint <- Clint/10^6*60 

  Qtotal.liverc <- Qtotal.liverc / parameters[['BW']]^0.25
  if (tolower(model) == "unscaled")
    CLh <- Clint
  else if (tolower(model) == "well-stirred")
    CLh <- Qtotal.liverc*fub*Clint/(Qtotal.liverc+fub*Clint)  
  else if (tolower(model) == "parallel tube")
    CLh <- Qtotal.liverc*(1-exp(-fub*Clint/Qtotal.liverc))
  else if (tolower(model) == "dispersion")
  {
    Rn <- fub*Clint/Qtotal.liverc
    a <- sqrt(1 + 4*Rn*Dn)
    CLh <- Qtotal.liverc*(1 - 4*a/((1+a)^2*exp((a-1)/2/Dn)-(1-a)^2*exp(-(a+1)/2/Dn)))
  }
  
  if(!suppress.messages) cat("Hepatic clearance calculated with the",hepatic.model,"model in units of L/h/kg.\n")
  
  return(as.numeric(CLh))
}

