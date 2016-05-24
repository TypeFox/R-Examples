## ----echo = F, message=F-------------------------------------------------
library(edeaR)
library(ggplot2)
library(dplyr)

## ------------------------------------------------------------------------
data(BPIC15_1)

## ----echo = F, fig.width = 7---------------------------------------------
activity_information <- activities(BPIC15_1)
ggplot(activity_information) +
	stat_ecdf(aes(absolute_frequency), lwd = 1, col = "#0072B2") + 
	scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
	xlab("Absolute activity frequencies") +
	ylab("Cumulative percentage") +
	geom_hline(y = 0.75)

## ------------------------------------------------------------------------
filtered_log <- filter_activity_frequency(BPIC15_1, percentile_cut_off = 0.25, reverse = T)
activities(filtered_log) %>% select(absolute_frequency) %>% summary

## ----fig.width = 7-------------------------------------------------------
case_throughput <- throughput_time(BPIC15_1, "case")
ggplot(case_throughput) +
	geom_histogram(aes(throughput_time), fill = "#0072B2", binwidth = 10) +
	xlab("Duration (in days)") +
	ylab("Number of cases")

## ----fig.width = 7-------------------------------------------------------
filtered_log <- filter_throughput_time(BPIC15_1, lower_threshold = 0, upper_threshold = 500)
case_throughput <- throughput_time(filtered_log, "case")
ggplot(case_throughput) +
	geom_histogram(aes(throughput_time), fill = "#0072B2", binwidth = 10) +
	xlab("Duration (in days)") +
	ylab("Number of cases")

## ----fig.width = 7-------------------------------------------------------
filtered_log <- filter_throughput_time(BPIC15_1, percentile_cut_off = 0.99, reverse = T)
throughput_time(filtered_log, "case")

## ----fig.width = 7-------------------------------------------------------
start <- BPIC15_1 %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
complete <- BPIC15_1 %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
bind_rows(start, complete) %>% 
	ggplot() +
	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*30) +
	facet_grid(type ~ .) +
	scale_fill_brewer(palette = "Dark2") +
	theme(legend.position = "none")

## ----fig.width = 7-------------------------------------------------------
library(lubridate)
a <- ymd("20120101")
b <- ymd("20120331")
filtered_log <- filter_time_period(BPIC15_1, a, b, "contained")
start <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
complete <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
bind_rows(start, complete) %>% 	
	ggplot() +
	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*7) +
	facet_grid(type ~ .) +
	scale_fill_brewer(palette = "Dark2") +
	theme(legend.position = "none")

## ----fig.width = 7-------------------------------------------------------
filtered_log <- filter_time_period(BPIC15_1, a, b, "start")
start <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
complete <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
bind_rows(start, complete) %>% 		ggplot() +
	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*7) +
	facet_grid(type ~ .) +
	scale_fill_brewer(palette = "Dark2") +
	theme(legend.position = "none")

## ----fig.width = 7-------------------------------------------------------
filtered_log <- filter_time_period(BPIC15_1, a, b, "complete")
start <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
complete <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
bind_rows(start, complete) %>% 	
	ggplot() +
	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*7) +
	facet_grid(type ~ .) +
	scale_fill_brewer(palette = "Dark2") +
	theme(legend.position = "none")

## ----fig.width = 7-------------------------------------------------------
filtered_log <- filter_time_period(BPIC15_1, a, b, "intersecting")
start <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
complete <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
bind_rows(start, complete) %>%
	ggplot() +
	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*7) +
	facet_grid(type ~ .) +
	scale_fill_brewer(palette = "Dark2") +
	theme(legend.position = "none")

## ----eval = F, fig.width = 7---------------------------------------------
#  filtered_log <- filter_time_period(BPIC15_1, a, b, "trim")
#  start <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = min(event_time.timestamp)) %>% mutate(type = "start")
#  complete <- filtered_log %>% group_by(case_concept.name) %>% summarize(timestamp = max(event_time.timestamp)) %>% mutate(type = "end")
#  bind_rows(start, complete) %>%
#  	ggplot() +
#  	geom_histogram(aes(timestamp, fill = type), binwidth = 60*60*24*7) +
#  	facet_grid(type ~ .) +
#  	scale_fill_brewer(palette = "Dark2") +
#  	theme(legend.position = "none")

