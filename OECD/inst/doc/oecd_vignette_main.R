## ----options, echo=FALSE-------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, warning = FALSE, error = FALSE, eval = TRUE)
library(OECD)

## ----loadLibrary, eval=FALSE---------------------------------------------
#  # from CRAN
#  install.packages("OECD")
#  
#  # from Github
#  library(devtools)
#  install_github("expersso/OECD")
#  
#  library(OECD)

## ----get_datasets, eval=FALSE--------------------------------------------
#  dataset_list <- get_datasets()
#  search_dataset("unemployment", data = dataset_list)

## ----dataset-------------------------------------------------------------
dataset <- "DUR_D"

## ----get_data_structure--------------------------------------------------
dstruc <- get_data_structure(dataset)
str(dstruc, max.level = 1)

## ----show_var_desc-------------------------------------------------------
dstruc$VAR_DESC

## ----breakdowns----------------------------------------------------------
dstruc$SEX
dstruc$AGE

## ----get_dataset---------------------------------------------------------
filter_list <- list(c("DEU", "FRA"), "MW", "2024")
df <- get_dataset(dataset = dataset, filter = filter_list)
head(df)

## ----duration------------------------------------------------------------
unique(df$DURATION)
dstruc$DURATION

## ----plot, cache = FALSE, fig.width = 7, fig.height = 4, fig.cap = "Unemployment rates of foreign- and native-born populations"----
df_plot <- df[df$DURATION == "UN5", ]
df_plot$obsTime <- as.numeric(df_plot$obsTime)

library(ggplot2)

qplot(data = df_plot, x = obsTime, y = obsValue, color = COUNTRY, geom = "line") +
  labs(x = NULL, y = "Persons, thousands", color = NULL,
       title = "Number of long-term unemployed men, aged 20-24")

## ----browse_metadata, eval=FALSE-----------------------------------------
#  browse_metadata(dataset)

## ----example-------------------------------------------------------------
df <- get_dataset("PATS_REGION",
                  filter = "PCT_A.INVENTORS.BEL+BE10.TOTAL+BIOTECH", 
                  pre_formatted = TRUE)
head(df)

