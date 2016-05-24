in.list <- function(chem.cas=NULL,
                    which.list="ToxCast")
{
  chem.lists <- chem.lists
  if (!(tolower(which.list) %in% tolower(names(chem.lists))))
    stop(paste("List",which.list,"not in available lists:",paste(names(chem.lists),collapse=", ")))

  if (!suppressWarnings(CAS.checksum(chem.cas))) stop (paste("Chemical CAS",chem.cas,"failes checksum."))
  return(chem.cas %in% chem.lists[tolower(names(chem.lists))==tolower(which.list)][[1]][,"CAS"])
}

is.nhanes.serum.parent <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.serum.parent"))
is.nhanes.serum.analyte <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.serum.analyte"))
is.nhanes.blood.parent <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.blood.parent"))
is.nhanes.blood.analyte <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.blood.analyte"))
is.nhanes.urine.parent <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.urine.parent"))
is.nhanes.urine.analyte <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES.urine.analyte"))
is.tox21 <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="Tox21"))
is.toxcast <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="ToxCast"))
is.expocast <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="ExpoCast"))
is.nhanes <- function(chem.cas) return(in.list(chem.cas=chem.cas,which.list="NHANES"))
is.httk <- function(chem.cas,species="Human", model = "3compartmentss") return(chem.cas %in% get_cheminfo(species=species,model = model))