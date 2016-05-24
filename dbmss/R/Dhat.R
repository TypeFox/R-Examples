Dhat <-
function(X, r = NULL, Cases, Controls = NULL, Intertype = FALSE, CheckArguments = TRUE) {
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  # K of cases.
  KCases <- Khat(X, r, Cases, Cases, CheckArguments = FALSE)
  # Default controls are all points except cases. Reserved name is "CoNtRoLs_"
  Y <- X
  if (is.null(Controls)) {
    Controls <- "CoNtRoLs_"
    if (Controls %in% levels(Y$marks$PointType)) stop("A point type is named 'CoNtRoLs_'. It must be changed to use the 'Controls = NULL' option of Dhat.")
    levels(Y$marks$PointType) <- c(levels(Y$marks$PointType), Controls)
    Y$marks$PointType[Y$marks$PointType != Cases] <- Controls
  }
  # K of controls. r must be those of cases.
  if (Intertype) {
    KControls <- Khat(Y, KCases$r, Cases, Controls, CheckArguments = FALSE)   
  } else {
    KControls <- Khat(Y, KCases$r, Controls, Controls, CheckArguments = FALSE)     
  }
  # Calculate the difference (a difference between fv's yields a dataframe)
  Dvalues <- KCases-KControls
  
  return (fv(cbind(as.data.frame(KCases)[1], Dvalues[2:3]), argu="r", ylab=quote(D(r)), valu=attr(KCases, "valu"), fmla=attr(KCases, "fmla"), alim=attr(KCases, "alim"), labl=c("r", "%s[ind](r)", "hat(%s)[iso](r)"), desc=attr(KCases, "desc"), unitname=attr(KCases, "unitname"), fname="D"))
}
