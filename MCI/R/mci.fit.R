mci.fit <-
function (mcidataset, submarkets, suppliers, shares, ..., origin=TRUE) {
mciworkfile <- mci.transmat (mcidataset, submarkets, suppliers, shares, ...)
mciworkfile_columns <- ncol(mciworkfile)
mciworkfile_rows <- nrow(mciworkfile)
mciworkfile_names <- names(mciworkfile)
mci_depvar <- mciworkfile_names[3]
mci_expvars <- mciworkfile_names[4:mciworkfile_columns]
mci_expvars_formula <- paste(mci_expvars, collapse=" + ")

if (origin == TRUE) {
mci_formula <- paste(mci_depvar, "~ 0 +", mci_expvars_formula)
}
else {
mci_formula <- paste(mci_depvar, "~", mci_expvars_formula)
}
mci_model <- lm (mci_formula, data = mciworkfile)
return(mci_model)
}
