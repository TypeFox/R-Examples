get_chem_id <- function(chem.cas=NULL,
                        chem.name=NULL)
{
  chem.physical_and_invitro.data <- chem.physical_and_invitro.data
  if (is.null(chem.cas) & is.null(chem.name))
  {
    stop("Must specifiy compound name or CAS.\n")
  } else if ((!is.null(chem.cas) & !any(chem.physical_and_invitro.data$CAS==chem.cas)) & (!is.null(chem.name) & !any(chem.physical_and_invitro.data$Compound==chem.name)))
  {
    stop("Compound not found.\n")
  }

  if (!is.null(chem.cas))
  {
#If chemical is identified by CAS, we must make sure its a valid CAS:
    if (!(chem.cas %in% chem.physical_and_invitro.data$CAS)) stop("CAS number not found, use get_cheminfo() for valid CAS numbers or set chem.name= argument.\n")
#Set the chemical name:
    found.chem.name <- chem.physical_and_invitro.data[chem.physical_and_invitro.data[,"CAS"]==chem.cas,"Compound"]
  }

  if (!is.null(chem.name))
  {
#If called by name, need to do a search to find the CAS number:
    names.index <- gsub("\\s","",tolower(chem.physical_and_invitro.data$Compound))
    names.index <- gsub("\\-","",names.index)
    name.key <- gsub("\\s","",tolower(chem.name))
    name.key <- gsub("\\-","",name.key)
    if (!any(names.index==name.key)) stop ("Chemical name not found.")
#Set the chemical CAS:
    found.chem.cas <- chem.physical_and_invitro.data[names.index==name.key,"CAS"]
    found.chem.cas <- found.chem.cas[!is.na(found.chem.cas)]
  }

  if (!is.null(chem.cas) & !is.null(chem.name))
  {
    if (chem.cas != found.chem.cas) stop(paste("Both CAS",chem.cas,"and name",chem.name,"were provided as arguments, but found other CAS -- ",found.chem.cas," -- for chemical,",chem.name))
    else if (tolower(chem.name) != tolower(found.chem.name)) warning(paste("Both CAS",chem.cas,"and name",chem.name,"were provided as arguments, but CAS also corresponds to name",found.chem.name))
  }

  if (is.null(chem.cas)) chem.cas <- found.chem.cas
  if (is.null(chem.name)) chem.name <- found.chem.name
  return(list(chem.cas=chem.cas,chem.name=chem.name))
}