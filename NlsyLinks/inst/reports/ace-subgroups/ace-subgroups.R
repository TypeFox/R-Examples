rm(list=ls(all=TRUE)) #Clear variables from previous runs.

# @knitr load_sources ------------------------------------------------------------

# @knitr load_packages -----------------------------------------------------------
library(magrittr)
library(knitr)
requireNamespace("NlsyLinks", quietly=T)
requireNamespace("xtable", quietly=T)
requireNamespace("dplyr", quietly=T)
requireNamespace("plyr", quietly=T)
requireNamespace("scales", quietly=T)

# @knitr define_globals ----------------------------------------------------------
oName <- "HeightZGenderAge" # o' stands for outcomes
# oName <- "WeightZGenderAge"

relationshipPaths <- 2
# relationshipPaths <- c(1, 2, 3, 4, 5)

rVersions <- c("R", "RFull", "RExplicit", "RImplicit",  "RImplicit2004")
dropIfHousematesAreNotSameGeneration <- FALSE
# rGroupsToDrop <- c()
rGroupsToDrop <- c(.375)
suppressGroupTables <- TRUE
# determinantThreshold <- 1e-5

# sql <- paste(
#   "SELECT 
#       Process.tblRelatedValuesArchive.AlgorithmVersion, 
#       Process.tblRelatedStructure.RelationshipPath, 
#       Process.tblRelatedValuesArchive.SubjectTag_S1, 
#       Process.tblRelatedValuesArchive.SubjectTag_S2,
#       Process.tblRelatedValuesArchive.RImplicitPass1, 
#       Process.tblRelatedValuesArchive.RImplicit, 
#       Process.tblRelatedValuesArchive.RImplicitSubject, 
#       Process.tblRelatedValuesArchive.RImplicitMother, 
#       Process.tblRelatedValuesArchive.RImplicit2004, 
#       Process.tblRelatedValuesArchive.RExplicitPass1, 
#       Process.tblRelatedValuesArchive.RExplicit, 
#       Process.tblRelatedValuesArchive.RPass1, 
#       Process.tblRelatedValuesArchive.R,
#       Process.tblRelatedValuesArchive.RFull, 
#       SameGeneration
#   FROM Process.tblRelatedValuesArchive 
#       INNER JOIN Process.tblRelatedStructure ON Process.tblRelatedValuesArchive.SubjectTag_S1 = Process.tblRelatedStructure.SubjectTag_S1 AND Process.tblRelatedValuesArchive.SubjectTag_S2 = Process.tblRelatedStructure.SubjectTag_S2 
#   WHERE Process.tblRelatedStructure.RelationshipPath IN ", relationshipPathsString, " 
#       AND (
#         Process.tblRelatedValuesArchive.AlgorithmVersion IN (
#           SELECT TOP (2) AlgorithmVersion 
#           FROM Process.tblRelatedValuesArchive AS tblRelatedValuesArchive_1 
#           GROUP BY AlgorithmVersion 
#           ORDER BY AlgorithmVersion DESC
#         )
#       )"
# )

# @knitr load_data ---------------------------------------------------------------
dsPair <- NlsyLinks::Links79PairExpanded
dsOutcomes <- NlsyLinks::ExtraOutcomes79[, c("SubjectTag", oName)]
dsDetails <- NlsyLinks::SubjectDetails79[, c("SubjectTag", "Gender", "RaceCohort")]

# @knitr tweak_data --------------------------------------------------------------
dsPair <- dsPair[as.integer(dsPair$RelationshipPath) %in% relationshipPaths, ]

oName_1 <- paste0(oName, "_S1")
oName_2 <- paste0(oName, "_S2")
relationshipPathsString <- paste0("(", paste(relationshipPaths, collapse=","), ")")
relationshipPathsPretty <- paste0("(", paste(levels(dsPair$RelationshipPath)[relationshipPaths], collapse=", "), ")")

# dsOutcomes$RandomFakeOutcome <- rnorm(n=nrow(dsOutcomes))
dsSubject <- dsOutcomes %>%
  dplyr::left_join(dsDetails, by="SubjectTag")
rm(dsOutcomes, dsDetails)

if( dropIfHousematesAreNotSameGeneration ) {
  dsRaw <- dsRaw[dsRaw$SameGeneration==1L, ]
}

#They're called 'dirty' because the final cleaning stage hasn't occurred yet (ie, removing unwanted R groups)
dsDirty <- NlsyLinks::CreatePairLinksSingleEntered(
  outcomeDataset   = dsSubject, 
  linksPairDataset = dsPair, 
  linksNames       = rVersions, 
  outcomeNames     = c(oName, "Gender", "RaceCohort")
)
rm(dsSubject, dsPair)

# dsDirty <- dsDirty[dsDirty$Gender_S1 == "Male"   & dsDirty$Gender_S2 == "Male", ]
# dsDirty <- dsDirty[dsDirty$Gender_S1 == "Female" & dsDirty$Gender_S2 == "Female", ]
# dsDirty <- dsDirty[
#   (dsDirty$Gender_S1 == "Male"   & dsDirty$Gender_S2 == "Female")
#   | dsDirty$Gender_S1 == "Female"   & dsDirty$Gender_S2 == "Male", ]

# dsDirty <- dsDirty[dsDirty$RaceCohort_S1 == "Nbnh"   & dsDirty$RaceCohort_S1 == "Nbnh", ]
# dsDirty <- dsDirty[dsDirty$RaceCohort_S1 == "Black"   & dsDirty$RaceCohort_S1 == "Black", ]
# dsDirty <- dsDirty[dsDirty$RaceCohort_S1 == "Hispanic"   & dsDirty$RaceCohort_S1 == "Hispanic", ]

# table(dsDirty$RFull)
# mean(!is.na(dsDirty$RFull)) 
# mean(!is.na(dsDirty[!is.na(dsDirty[, oName_1]) & !is.na(dsDirty[, oName_2]), "RFull"])) 

# @knitr evaluate_groups ---------------------------------------------------------
groupDatasets <- list() 
dsAce <- data.frame(
  Version=rVersions, 
  ASq=NA_real_, CSq=NA_real_, ESq=NA_real_, 
  ASqSE=NA_real_, CSqSE=NA_real_, ESqSE=NA_real_, 
  N=NA_integer_
)
for( i in seq_along(rVersions) ) {
  rVersion <-  rVersions[i] # rVersion <- "RFull"
#  print(rVersion)
  dsGroupSummary <- NlsyLinks::RGroupSummary(dsDirty, oName_1, oName_2, rName=rVersion)#, determinantThreshold=determinantThreshold)
  dsGroupSummary[dsGroupSummary[, rVersion] %in% rGroupsToDrop, "Included"] <- FALSE
  
  groupDatasets[[(i-1)*2 + 1]] <- dsGroupSummary
       
  dsClean <- NlsyLinks::CleanSemAceDataset(dsDirty=dsDirty, dsGroupSummary, oName_1, oName_2, rName=rVersion)
  
  ace <- NlsyLinks::AceLavaanGroup(dsClean)
  est <- lavaan::parameterEstimates(ace@Details$lavaan)
  aSqSE <- est[est$label=="a2", "se"]
  cSqSE <- est[est$label=="c2", "se"]
  eSqSE <- est[est$label=="e2", "se"]
  
  dsAce[i, 2:8] <- c(
    ace@ASquared, ace@CSquared, ace@ESquared, 
    aSqSE, cSqSE, eSqSE, 
    ace@CaseCount
  )
}
# dsAce
# groupDatasets[[1]]

PrintGroupSummary <- function( dSummary, title="Group Summary"  ) {
  dSummary <- dplyr::rename_(dSummary,
    "Included in SEM" = "Included",
    "$N_{Pairs}$"     = "PairCount",
    "$\\bar{x}_1$"    = "O1Mean",
    "$\\bar{x}_2$"    = "O2Mean",
    "$s_1^2$"         = "O1Variance",
    "$s_2^2$"         = "O2Variance",
    "$s_{1,2}$"       = "O1O2Covariance",
    "$r$"             = "Correlation"
  )
  
  digitsFormat <- c(0,3,0,0,2,2,2,2,2,2,1,0) #Include a dummy at the beginning, for the row.names.
  textTable <- xtable::xtable(dSummary, caption=title, digits=digitsFormat)
  xtable:::print.xtable(textTable, include.rownames=F, sanitize.text.function = function(x) {x}) # size="large", size="small",
}
# PrintGroupSummary(dSummary=dsGroupSummary)

for( i in seq_along(rVersions) ) {
  rVersion <- rVersions[i]
  # if( startNewPage[i] ) cat('\\newpage \n')
  cat('\\section{Subgroups -- ', rVersion, '}')
  
  PrintGroupSummary(groupDatasets[[(i-1)*2 + 1]], title=rVersion)
}

# @knitr evaluate_ace ------------------------------------------------------------
PrintAces <- function( ) {
  dAcePretty <- dsAce
  dAcePretty <- dplyr::rename_(dAcePretty,
    "$R$ Variant" = "Version", 
    "$a^2$"       = "ASq", 
    "$c^2$"       = "CSq", 
    "$e^2$"       = "ESq",  
    "$se_{a^2}$"   = "ASqSE", 
    "$se_{c^2}$"   = "CSqSE", 
    "$se_{e^2}$"   = "ESqSE", 
    "$N$"         = "N"
  )
  
  dAcePretty <- cbind(dAcePretty[, 1], plyr::numcolwise(round)(dAcePretty,digits=2))
  dAcePretty <- cbind(dAcePretty[, 1], plyr::numcolwise(scales::comma)(dAcePretty))

  dAcePretty <- t(apply(dAcePretty, 1, function(x) gsub("^0.", ".", x)))
  
  digitsFormat <- 2# c(0,1,3,3,3,3,3,3) #Include a dummy at the beginning, for the row.names.
  alignnment <- "ll|rrr|rrr|r" #Include an initial dummy for the (suppressed) row names.
  textTable <- xtable::xtable(dAcePretty, caption="Comparison of R Variants (by rows) and of Links Versions (left vs right side).", digits=digitsFormat, align=alignnment)
  xtable:::print.xtable(textTable, include.rownames=F, size="large", sanitize.text.function = function(x) {x}, floating=T)
}
PrintAces()
# summary(dsAce)
