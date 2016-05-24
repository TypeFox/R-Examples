##                       getMETA                      ##
##      This code is part of the rusda package        ##
##      F.-S. Krah 2015 (last update: 2015-07-11)     ##

getMETA <- function(x)
{
  dbs <- grep("This report contains data from the following databases:", x)
  dbs <- x[dbs]
  n <- grep("Nomenclature", dbs)
  hf <- grep("Fungus-Host", dbs)
  sp <- grep("Specimens", dbs)
  lit <- grep("Literature", dbs)
  
  n2 <- grep("Nomenclature data for ", x)
  sp2 <- grep("The Specimens database has ", x)
  hf2 <- grep("The Fungus-Host Distributions database has ", x)
  lit2 <- grep("The Literature database has ", x)
  n <- length(c(n, n2))
  lit <- length(c(lit, lit2 ))
  sp <- length(c(sp,sp2))
  hf <- length(c(hf, hf2))
  c(Nomenclature = (n > 0), Specimens = (sp > 0), 
    Host_Fungus=(hf > 0), Literature=(lit > 0))*1
}