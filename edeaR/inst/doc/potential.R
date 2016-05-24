## ----message = F, warning = F--------------------------------------------
library(edeaR)
library(dplyr)
library(ggplot2)
data("BPIC14_incident_log")
data("BPIC14_incident_case_attributes")
BPIC14_incident_log %>% print
BPIC14_incident_case_attributes %>% print

## ------------------------------------------------------------------------
trace_coverage(BPIC14_incident_log, level_of_analysis = "trace") %>% print(width = Inf)

## ----fig.width = 7-------------------------------------------------------
trace_length(BPIC14_incident_log, "trace") %>% ggplot(aes(relative_trace_frequency, absolute)) + geom_jitter() + scale_x_continuous(limits = c(0,0.01)) + ylab("Trace length") + xlab("Relative trace frequency")

## ------------------------------------------------------------------------
case_performance <- BPIC14_incident_log %>% throughput_time("case") %>% 
	left_join(BPIC14_incident_log %>% trace_length("case")) %>% 
	left_join(BPIC14_incident_log %>% repetitions("case") %>% select(incident_id, absolute) %>% rename(repetitions = absolute)) %>% 
	left_join(BPIC14_incident_log %>% number_of_selfloops("case") %>% select(incident_id, absolute) %>% rename(number_of_selfloops = absolute)) 

case_performance %>% summary

## ----fig.width = 7, fig.align='center'-----------------------------------
input <- scale(case_performance[,2:5]) %>% as.data.frame()

clusters <- data.frame(i = 1:15)
for(i in 1:nrow(clusters)) {
	for(j in 1:100) {
		cl <- kmeans(input,i, iter.max = 20)
		min_sse <- min(Inf, cl$totss - cl$betweenss)
	}
	clusters$sse[i] <- min_sse
}
clusters %>% ggplot(aes(i, sse)) + geom_line() +
	xlab("Number of clusters") + ylab("SSE") + 
	scale_x_continuous(breaks = 1:15) + scale_y_continuous(breaks = seq(0,16000,2000))

## ------------------------------------------------------------------------
set.seed(4)
cl <- kmeans(input,5, iter.max = 20)
case_performance <- case_performance %>% bind_cols(data.frame(cluster = factor(cl$cluster)))
cl %>% str

## ----fig.align="center", fig.width= 7, echo = F--------------------------
case_performance %>%
	group_by(cluster) %>% 
	summarize(freq = n(),
			  mean_nr_of_repetitions = mean(repetitions),
			  mean_nr_of_selfloops = mean(number_of_selfloops),
			  mean_trace_length = mean(trace_length),
			  mean_throughput_time = mean(throughput_time)) %>%
	print(width = Inf)

case_performance %>% ggplot(aes(cluster, repetitions)) + geom_boxplot(aes(fill = cluster))
case_performance %>% ggplot(aes(cluster, number_of_selfloops)) + geom_boxplot(aes(fill = cluster))
case_performance %>% ggplot(aes(cluster, trace_length)) + geom_boxplot(aes(fill = cluster))
case_performance %>% ggplot(aes(cluster, throughput_time)) + geom_boxplot(aes(fill = cluster))

## ------------------------------------------------------------------------
BPIC14_incident_case_attributes <- BPIC14_incident_case_attributes %>%
	merge(case_performance) 

## ----fig.width = 7-------------------------------------------------------
BPIC14_incident_case_attributes %>% ggplot(aes(reorder(ci_type_aff, as.numeric(cluster) == 4, FUN = "mean"), fill = cluster)) + geom_bar(position = "fill")  +
	scale_fill_brewer() + coord_flip() + xlab("ci_type_aff") + scale_y_continuous(breaks = seq(0,1,0.1))

## ----fig.width = 7-------------------------------------------------------
BPIC14_incident_case_attributes %>% ggplot(aes(cluster, x_reassignments)) + geom_boxplot()  +
	scale_fill_brewer() + coord_flip() 

## ----fig.width = 7-------------------------------------------------------
BPIC14_incident_case_attributes %>% ggplot(aes(reorder(closure_code, as.numeric(cluster) == 4, FUN = "mean"), fill = cluster)) + geom_bar(position = "fill")  +
	scale_fill_brewer() + coord_flip() + xlab("closure_code") + scale_y_continuous(breaks = seq(0,1,0.1))

## ------------------------------------------------------------------------
cluster_5 <- BPIC14_incident_log %>% 
	merge(case_performance) %>%
	filter(cluster == 5) %>%
	eventlog(case_id(BPIC14_incident_log), 
			 activity_id(BPIC14_incident_log),
			 activity_instance_id(BPIC14_incident_log),
			 lifecycle_id(BPIC14_incident_log),
			 timestamp(BPIC14_incident_log))

## ------------------------------------------------------------------------
cluster_5 
cluster_5 %>% traces

