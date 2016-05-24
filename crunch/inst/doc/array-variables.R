## ---- results='hide', echo=FALSE, message=FALSE--------------------------
library(crunch)
load("vignettes.RData")

## ---- eval=FALSE---------------------------------------------------------
#  grep("^imiss_", names(ds), value=TRUE)

## ---- echo=FALSE---------------------------------------------------------
grep("^imiss_", names(start_make_array), value=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  ds$imiss_b

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss_b, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$imiss <- makeArray(pattern="^imiss_", dataset=ds, name="Issue importance")
#  ds$imiss

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  subvariables(ds$imiss)

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss_subvars, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$imiss$imiss_b

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss_b, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  names(subvariables(ds$imiss))

## ---- echo=FALSE---------------------------------------------------------
names_imiss_subvars

## ---- eval=FALSE---------------------------------------------------------
#  names(subvariables(ds$imiss)) <- c("The economy", "Immigration",
#      "The environment", "Terrorism", "Gay rights", "Education",
#      "Health care", "Social security", "The budget deficit",
#      "The war in Afghanistan", "Taxes", "Medicare", "Abortion")
#  subvariables(ds$imiss)

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss_subvars2, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  sorting <- order(names(subvariables(ds$imiss)))
#  subvariables(ds$imiss) <- subvariables(ds$imiss)[sorting]
#  subvariables(ds$imiss)

## ---- echo=FALSE---------------------------------------------------------
cat(show_imiss_subvars3, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$boap_4

## ---- echo=FALSE---------------------------------------------------------
cat(show_boap_4, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$boap <- makeMR(pattern="^boap_[0-9]+", dataset=ds,
#      name="Approval of Obama on issues",
#      selections=c("Strongly approve", "Somewhat approve"))
#  ds$boap

## ---- echo=FALSE---------------------------------------------------------
cat(show_boap, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$boap <- undichotomize(ds$boap)
#  ds$boap

## ---- echo=FALSE---------------------------------------------------------
cat(show_boap2, sep="\n")

## ---- eval=FALSE---------------------------------------------------------
#  ds$boap <- dichotomize(ds$boap, "Strongly approve")
#  ds$boap

## ---- echo=FALSE---------------------------------------------------------
cat(show_boap3, sep="\n")

## ------------------------------------------------------------------------
grep("boap", names(ds), value=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  unbind(ds$boap)
#  ds <- refresh(ds)

## ---- results='hide', echo=FALSE, message=FALSE--------------------------
ds <- start_make_array

## ------------------------------------------------------------------------
grep("boap", names(ds), value=TRUE)

