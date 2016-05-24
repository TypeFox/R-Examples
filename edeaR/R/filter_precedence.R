#' @title Filter: precedence relations
#'
#' @description Filters cases based on the precedence relations between two sets of activities: antecedents and consequent.
#' The filter can detect directly following activities as well as eventually following activites.
#'
#' @param eventlog The event log to be used. An object of class
#' \code{eventlog}.
#'
#' @param antecedents,consequents The set of antecendent and consequent activities. All pairs of antecedents and consequents are checked for.
#'
#' @param precedence_type When \code{directly_follows}, the consequent activity should happen immediately after the antecedent activities.
#' When \code{eventually_follows}, other events are allowed to happen in between.
#'
#' @param filter_method When \code{each}, only cases where all the relations are valid are preserved. When \code{one_of}, all the cases where
#' at least one of the conditions hold are preserved.
#'
#' @param reverse A logical parameter depicting whether the selection should be reversed.
#'
#' @export filter_precedence


filter_precedence <- function(eventlog,
							  antecedents,
							  consequents,
							  precedence_type,
							  filter_method,
							  reverse = F) {
	stop_eventlog(eventlog)

	if(!(precedence_type %in% c("directly_follows","eventually_follows")))
		stop("Precedence type should be one of the following: directly_follows, eventually_follows")
	if(precedence_type == "directly_follows")
		interleavings_allowed = FALSE
	else
		interleavings_allowed = TRUE

	if(!(filter_method %in% c("each","one_of")))
		stop("filter_method should be one of the following: each, one_of")

	sequences <- paste(rep(antecedents, each = length(consequents)),
					   rep(consequents, times = length(antecedents)), sep = ",")
	number_of_conditions <- length(sequences)

	patterns <- data.frame(pattern = sequences)
	tr <- traces(eventlog)

	dummies <- generate_pattern_dummies(patterns, eventlog, interleavings_allowed = interleavings_allowed)
	colnames(dummies)[colnames(dummies) == case_id(eventlog)] <- "case_classifier"
	colnames(eventlog)[colnames(eventlog) == case_id(eventlog)] <- "case_classifier"

	dummies$conditions_valid <- rowSums(select(dummies, starts_with("X")))
	if(filter_method == "one_of")
		case_selection <- filter(dummies, conditions_valid > 0)$case_classifier
	else
		case_selection <- filter(dummies, conditions_valid == number_of_conditions)$case_classifier

	if(reverse == FALSE)
		f_eventlog <- filter(eventlog, case_classifier %in% case_selection)
	else
		f_eventlog <- filter(eventlog, !(case_classifier %in% case_selection))

	colnames(f_eventlog)[colnames(f_eventlog)=="case_classifier"] <- case_id(eventlog)

	output <- eventlog(f_eventlog,
					   activity_id = activity_id(eventlog),
					   case_id = case_id(eventlog),
					   timestamp =timestamp(eventlog),
					   lifecycle_id = lifecycle_id(eventlog),
					   activity_instance_id = activity_instance_id(eventlog))
	return(output)

}
