## ------------------------------------------------------------------------
library(radarchart)
skills

## ------------------------------------------------------------------------
skillsByName

## ---- tidyR--------------------------------------------------------------
library(tidyr)
skillsByLabel <- gather(skillsByName, key=Label, value=Score, -Name) %>%
                   spread(key=Name, value=Score)
skillsByLabel

## ---- baseR--------------------------------------------------------------
skillsByLabel <- as.data.frame(t(skillsByName[,-1]))
names(skillsByLabel) <- skillsByName$Name
skillsByLabel <- cbind(Label=row.names(skillsByLabel), skillsByLabel)
row.names(skillsByLabel) <- NULL

## ---- eval=FALSE---------------------------------------------------------
#  chartJSRadar(scores = skillsByLabel, maxScale = 10)

