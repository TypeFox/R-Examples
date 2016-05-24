# This function retrives physico-chemical properties ("param") for the chemical specified by chem.name or chem.CAS from the vLiver tables.
get_physchem_param <- function(param,chem.name=NULL,chem.CAS=NULL)
{
  Wetmore.data <- Wetmore.data
  chem.physical_and_invitro.data <- chem.physical_and_invitro.data
  if (is.null(chem.CAS) & is.null(chem.name))
  {
    stop("Must specifiy compound name or CAS.\n")
  } else if ((!is.null(chem.CAS) & !any(Wetmore.data[,"CAS"]==chem.CAS)) & (!is.null(chem.name) & !any(Wetmore.data[,"Compound"]==chem.name)))
  {
    stop("Compound not found.\n")
  } else {
    if (!is.null(chem.CAS)) chem.name <- chem.physical_and_invitro.data[chem.physical_and_invitro.data[,"CAS"]==chem.CAS,"Compound"][1]
    else chem.CAS <- chem.physical_and_invitro.data[chem.physical_and_invitro.data[,"Compound"]==chem.name,"CAS"][1]

    if (!(param %in% c("MW","logP","pKa_Donor","pKa_Accept",'logMA'))) stop(paste("Parameter",param,"not among \"MW\", \"logP\", \"logMA\", \"pKa_Donor\", and \"pKa_Accept\".\n"))
    chem.physical_and_invitro.data.index <- which(chem.physical_and_invitro.data$CAS==chem.CAS)
    
    if (!is.na(as.numeric(chem.physical_and_invitro.data[chem.physical_and_invitro.data.index,param])) | (param %in% c("pKa_Donor","pKa_Accept","logMA")))
    {
      value <- chem.physical_and_invitro.data[chem.physical_and_invitro.data.index,param] 
      if(!is.na(value) & param %in% c('pKa_Donor','pKa_Accept') & regexpr(",",value)!=-1) value <- sort(as.numeric(strsplit(value,",")[[1]]))
      return(as.numeric(value))
    }else stop(paste("Incomplete phys-chem data for ",chem.name," -- missing ",param,"."))
  }
}

