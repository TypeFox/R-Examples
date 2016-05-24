create_analytics <- function(container,classification_results,b=1) {	
	create_documentSummary <- function(container, scores) {
		testing_codes <- as.numeric(as.vector(container@testing_codes))
		return(cbind(MANUAL_CODE=testing_codes,CONSENSUS_CODE=scores$BEST_LABEL,CONSENSUS_AGREE=scores$NUM_AGREE,CONSENSUS_INCORRECT=testing_codes!=scores$BEST_LABEL,PROBABILITY_CODE=as.numeric(as.vector(scores$BEST_PROB)),PROBABILITY_INCORRECT=testing_codes!=scores$BEST_PROB))
	}
	
	create_topicSummary <- function(container, scores) {
		topic_codes <- unique(c(levels(container@training_codes),levels(container@testing_codes)))
		testing_codes <- as.numeric(as.vector(container@testing_codes))

		manually_coded <- table(factor(container@testing_codes,levels=topic_codes))
		automatically_coded_label <- table(factor(scores$BEST_LABEL,levels=topic_codes))#sapply(topic_codes,compare,scores$BEST_LABEL)
		automatically_coded_prob <- table(factor(scores$BEST_PROB,levels=topic_codes))#sapply(topic_codes,compare,scores$BEST_PROB)
		correctly_coded_label <- table(factor(container@testing_codes[which(testing_codes==scores$BEST_LABEL)],levels=topic_codes))/manually_coded*100

		correctly_coded_prob <- table(factor(testing_codes[which(testing_codes==scores$BEST_PROB)],levels=topic_codes))/manually_coded*100
		over_coded_label <- automatically_coded_label/manually_coded*100
		over_coded_prob <- automatically_coded_prob/manually_coded*100

		return(cbind(TOPIC_CODE=as.numeric(as.vector(topic_codes)),NUM_MANUALLY_CODED=manually_coded,NUM_CONSENSUS_CODED=automatically_coded_label,NUM_PROBABILITY_CODED=automatically_coded_prob,PCT_CONSENSUS_CODED=over_coded_label,PCT_PROBABILITY_CODED=over_coded_prob,PCT_CORRECTLY_CODED_CONSENSUS=correctly_coded_label,PCT_CORRECTLY_CODED_PROBABILITY=correctly_coded_prob))
	}
	
	if (container@virgin == FALSE) {
		score_summary <- create_scoreSummary(container, classification_results)
		document_summary <- create_documentSummary(container, score_summary)
		topic_summary <- as.data.frame(create_topicSummary(container, score_summary))
		algorithm_summary <- as.data.frame(create_precisionRecallSummary(container, classification_results, b))
		
		topic_summary <- topic_summary[with(topic_summary, order(TOPIC_CODE)),]
		row.names(topic_summary) <- topic_summary[,1]
		
		raw_summary <- cbind(classification_results,document_summary)
		
		ensemble_summary <- create_ensembleSummary(as.data.frame(raw_summary))
		
		container <- new("analytics", label_summary=as.data.frame(topic_summary)[,-1], document_summary=as.data.frame(raw_summary), algorithm_summary=as.data.frame(algorithm_summary), ensemble_summary=ensemble_summary)
    } else {
		score_summary <- create_scoreSummary(container, classification_results)
		
		document_summary <- create_documentSummary(container, score_summary)
		document_summary <- document_summary[,c(2,3,5)]

		raw_summary <- cbind(classification_results, document_summary)
		
		topic_summary <- create_topicSummary(container, score_summary)
		topic_summary <- as.data.frame(topic_summary[,c(1,3,4)])
		topic_summary <- topic_summary[with(topic_summary, order(TOPIC_CODE)),]
		row.names(topic_summary) <- topic_summary[,1]
		
		container <- new("analytics_virgin", label_summary=as.data.frame(topic_summary)[,-1], document_summary=as.data.frame(raw_summary))
	}
	
    return(container)   
}