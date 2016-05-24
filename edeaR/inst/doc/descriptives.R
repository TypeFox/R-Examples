## ----message = F---------------------------------------------------------
library(edeaR)
data("BPIC15_1")

## ------------------------------------------------------------------------
summary(BPIC15_1)

## ------------------------------------------------------------------------
case_information <- cases(BPIC15_1)
case_information

## ------------------------------------------------------------------------
library(dplyr)
summary(select(case_information, first_activity, last_activity))

## ----fig.width=7---------------------------------------------------------
library(ggplot2)
ggplot(case_information) + 
	geom_bar(aes(duration_in_days), binwidth = 30, fill = "#0072B2") + 
	scale_x_continuous(limits = c(0,500)) +
	xlab("Duration (in days)") + 
	ylab("Number of cases") 

## ------------------------------------------------------------------------
activity_information <- activities(BPIC15_1)
activity_information

## ----fig.width = 7-------------------------------------------------------
ggplot(activity_information) +
	stat_ecdf(aes(absolute_frequency), lwd = 1, col = "#0072B2") + 
	scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
	xlab("Absolute activity frequencies") +
	ylab("Cumulative percentage")

## ------------------------------------------------------------------------
activity_selfloops <- number_of_selfloops(BPIC15_1, level_of_analysis = "activity")
activity_selfloops

## ----fig.width=7, fig.height=4-------------------------------------------
ggplot(activity_selfloops) + 
	geom_bar(aes(reorder(event_concept.name, -absolute), absolute), stat = "identity", fill = "#0072B2") + 
	theme(axis.text.x = element_text(angle = 90)) + 
	xlab("Activity") + 
	ylab("Number of selfloops")

## ------------------------------------------------------------------------
activity_repetitions <- repetitions(BPIC15_1, level_of_analysis = "activity")
activity_repetitions

## ----fig.width=7, fig.height=4-------------------------------------------
ggplot(activity_repetitions) + 
	geom_bar(aes(reorder(event_concept.name, -absolute), absolute), stat = "identity", fill = "#0072B2") + 
	theme(axis.text.x = element_text(angle = 90)) + 
	xlab("Activity") + 
	ylab("Number of repetitions")

## ----fig.width=7, fig.height=7-------------------------------------------
data <- bind_rows(mutate(activity_selfloops, type = "selfloops"),
			  mutate(select(activity_repetitions, event_concept.name, absolute), type = "repetitions"))

ggplot(data) + 
	geom_bar(aes(reorder(event_concept.name, -absolute), absolute), stat = "identity", fill = "#0072B2") + 
	facet_grid(type ~ .) +
	theme(axis.text.x = element_text(angle = 90)) + 
	xlab("Activity") + 
	ylab("Number of selfloops and repetitions")

