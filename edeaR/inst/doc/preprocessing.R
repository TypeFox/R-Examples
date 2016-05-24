## ----echo = F, message=F, warning=FALSE----------------------------------
library(edeaR)
library(dplyr)

## ----eval = F------------------------------------------------------------
#  data <- eventlog_from_xes()

## ----echo = F------------------------------------------------------------
data("BPIC15_1_imported")
data <- BPIC15_1_imported

## ------------------------------------------------------------------------
data

## ------------------------------------------------------------------------
table(data$event_lifecycle.transition)

## ------------------------------------------------------------------------
n_events(data)
n_activity_instances(data)

## ----results = "hold"----------------------------------------------------
case_id(data)
activity_id(data)
activity_instance_id(data)
lifecycle_id(data)
timestamp(data)

## ------------------------------------------------------------------------
library(lubridate)
data[1:4,timestamp(data)]

data$event_time.timestamp <- ymd_hms(data$event_time.timestamp)

## ----cache = F-----------------------------------------------------------
data("csv_example", package = "edeaR")

## ------------------------------------------------------------------------
head(csv_example)

## ------------------------------------------------------------------------
csv_example$activity_instance <- 1:nrow(csv_example)

## ----eval = F------------------------------------------------------------
#  library(tidyr)
#  csv_example <- gather(csv_example, LIFECYCLE, TIMESTAMP, -CASE, - ACTIVITY, -ACTIVITY_INSTANCE)
#  head(csv_example)

## ----echo = F------------------------------------------------------------
data(example_log)
example_log <- as.data.frame(example_log)
example_log$LIFE_CYCLE <- factor(example_log$LIFE_CYCLE, labels = c("START","COMPLETE"))
csv_example <- example_log %>% rename(LIFECYCLE = LIFE_CYCLE)
head(csv_example)

## ------------------------------------------------------------------------
csv_example$LIFECYCLE <- factor(csv_example$LIFECYCLE, labels = c("start","complete"))
head(csv_example)

## ------------------------------------------------------------------------
csv_example$TIMESTAMP <- ymd_hms(csv_example$TIMESTAMP)

## ------------------------------------------------------------------------
log <- eventlog(eventlog = csv_example, 
				case_id = "CASE",
				activity_id = "ACTIVITY", 
				activity_instance_id = "ACTIVITY_INSTANCE", 
				lifecycle_id = "LIFECYCLE", 
				timestamp = "TIMESTAMP")

log

