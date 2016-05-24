# This function retrives in vitro PK data (e.g. intrinsic metabolic clearance or fraction unbound in plasma) from the vLiver tables.
get_invitroPK_param <- function(param,species,chem.name=NULL,chem.CAS=NULL)
{
  chem.physical_and_invitro.data <- chem.physical_and_invitro.data
  if (is.null(chem.CAS) & is.null(chem.name))
  {
    stop("Must specifiy compound name or CAS.\n")
  } else if ((!is.null(chem.CAS) & !any(chem.physical_and_invitro.data[,"CAS"]==chem.CAS)) & (!is.null(chem.name) & !any(chem.physical_and_invitro.data[,"Compound"]==chem.name)))
  {
    stop("Compound not found.\n")
  } else {
    if (!is.null(chem.CAS)) chem.name <- chem.physical_and_invitro.data[chem.physical_and_invitro.data[,"CAS"]==chem.CAS,"Compound"][1]
    else chem.CAS <- chem.physical_and_invitro.data[chem.physical_and_invitro.data[,"Compound"]==chem.name,"CAS"][1]
    if (!(param %in% c("Clint","Funbound.plasma","Clint.pValue","Fgutabs"))) stop(paste("Parameter",param,"not among \"Clint\", \"Clint.pValue\", and \"Funbound.plasma\".\n"))

    chem.physical_and_invitro.data.index <- which(chem.physical_and_invitro.data$CAS==chem.CAS)
    this.col.name <- tolower(paste(species,param,sep="."))
    if (!(this.col.name %in% tolower(colnames(chem.physical_and_invitro.data))))
    {
      warning(paste("No in vitro ",param," data for ",chem.name," in ",species,".",sep=""))
      for (alternate.species in c("Human","Rat","Mouse","Dog","Monkey","Rabbit"))
      {
        this.col.name <- tolower(paste(alternate.species,param,sep="."))
        if (this.col.name %in% tolower(colnames(chem.physical_and_invitro.data)))
        {
          warning(paste("Substituting ",alternate.species," in vitro ",param," data for ",chem.name," ",species,".",sep=""))
          break()
        }
      }
    }
    if (this.col.name %in% tolower(colnames(chem.physical_and_invitro.data)))
    {
      this.col.index <- which(tolower(colnames(chem.physical_and_invitro.data))==this.col.name)
      if (!is.na(as.numeric(chem.physical_and_invitro.data[chem.physical_and_invitro.data.index,this.col.index])) | param=="Clint.pValue")
        return(as.numeric(chem.physical_and_invitro.data[chem.physical_and_invitro.data.index,this.col.index]))
    }
    stop(paste("Incomplete in vitro PK data for ",chem.name," in ",species," -- missing ",param,".",sep=""))
  }
}