## ----preliminaries, echo=FALSE, message=FALSE---------------------------------
options(width=80, continue=" ")
library(NlsyLinks)

## ----echo=TRUE----------------------------------------------------------------
subset(Links79Pair, RelationshipPath=='Gen2Siblings' & R==.75)

## ----eval=FALSE, echo=TRUE, tidy=FALSE----------------------------------------
#  dsLinks <- Links79PairExpanded
#  isGen2Sib <- dsLinks$RelationshipPath=='Gen2Siblings'
#  
#  olderFullYoungerHalf <- (dsLinks$RExplicitOlderSibVersion==.5 &
#                             dsLinks$RExplicitYoungerSibVersion==.25)
#  
#  olderHalfYoungerFull <- (dsLinks$RExplicitOlderSibVersion==.25 &
#                             dsLinks$RExplicitYoungerSibVersion==.5)
#  
#  
#  dsLinks[isGen2Sib & (olderFullYoungerHalf | olderHalfYoungerFull), ]

## ----eval=FALSE, echo=TRUE, tidy=FALSE----------------------------------------
#  dsLinks[ isGen2Sib & (dsLinks$RExplicitOlderSibVersion==.375 |
#       dsLinks$RExplicitYoungerSibVersion==.375), ]

## ----echo=TRUE----------------------------------------------------------------
library(NlsyLinks)
filePathOutcomes <- file.path(path.package("NlsyLinks"), "extdata", "gen1-life-course.csv")

## ----echo=TRUE----------------------------------------------------------------
dsDemographics <- ReadCsvNlsy79Gen1(filePathOutcomes)
summary(dsDemographics)

## ----echo=TRUE----------------------------------------------------------------
dsDemographics <- RenameNlsyColumn(dsDemographics, "R0214700", "Race")
dsDemographics <- RenameNlsyColumn(dsDemographics, "R0214800", "Gender")

## ----echo=TRUE, tidy=FALSE----------------------------------------------------
#The official documentation calls this last level "NON-BLACK, NON-HISPANIC"
dsDemographics$Race <- factor(x=dsDemographics$Race,
                          levels=1:3,
                          labels=c("Hispanic", "Black", "NBNH"))
dsDemographics$Gender <- factor(x=dsDemographics$Gender,
                          levels=1:2,
                          labels=c("Male", "Female"))

## ----echo=TRUE, tidy=FALSE----------------------------------------------------
library(NlsyLinks)
filePathOutcomes <- file.path(path.package("NlsyLinks"), "extdata", "gen2-birth.csv")
dsDemographics <- ReadCsvNlsy79Gen2(filePathOutcomes)  #Notice this function is for Gen2.
# summary(dsDemographics) #Uncomment to see the summary

dsDemographics <- RenameNlsyColumn(dsDemographics, "C0005300", "Race")
dsDemographics <- RenameNlsyColumn(dsDemographics, "C0005400", "Gender")

dsDemographics$Race <- factor(x=dsDemographics$Race,
                          levels=1:3,
                          labels=c("Hispanic", "Black", "NBNH"))
dsDemographics$Gender <- factor(x=dsDemographics$Gender,
                          levels=1:2,
                          labels=c("Male", "Female"))

## ----echo=TRUE, eval=FALSE, tidy=FALSE----------------------------------------
#  dsCombined <- merge(x=dsDemographics, y=dsOutcomes, by="SubjectTag", all=TRUE)

