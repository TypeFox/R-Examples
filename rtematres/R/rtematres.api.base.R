#' Basic tematres server api
#'
#' @param task The api task you like to execute. Use the the "availableTasks"
#' 	  to get an overview about the base api. It returns a data frame with
#'        descriptions and the arguments for the tasks.
#'
#' @param argument Is the argument for the api task. You find the information
#'        about the arguments when you call the task "availableTasks".
#'
#' @return The function returns either a dataframe for information or a list
#'         of keywords and ids
#' @examples \dontrun{
#'     rtematres.api(task = "availableTasks")
#'     rtematres.api(task = "fetchVocabularyData")
#'     rtematres.api(task = "fetchTopTerms")
#'     rtematres.api(task = "search", argument = "measurement")
#'     rtematres.api(task = "letter", argument = "t")
#'     rtematres.api(task = "fetchTerm", argument = 12)
#'     rtematres.api(task = "fetchDown", argument = 4 )
#'     rtematres.api(task = "fetchUp", argument = 4)
#'     rtematres.api(task = "fetchRelated", argument = 4)
#'     rtematres.api(task = "fetchAlt", argument = 12 )
#'     rtematres.api(task = "fetchCode", argument = "tree")
#'     rtematres.api(task = "fetchNotes", argument = 5 )
#'     rtematres.api(task = "fetchDirectTerms", argument = 12)
#'     rtematres.api(task = "fetchURI", argument = 12)
#'     rtematres.api(task = "fetchTargetTerms", argument = 12 )
#'     rtematres.api(task = "fetchSourceTerm", argument = 12)
#'     rtematres.api(task = "fetchTerms", argument = '12,13' )
#'     rtematres.api(task = "fetchRelatedTerms", argument = '12,13' )
#'     rtematres.api(task = "fetchSimilar", argument = 12)
#'     rtematres.api(task = "fetchLast")
#'     rtematres.api.conversion.keyword_id(given = "Measurement")
#'     rtematres.api.conversion.keyword_id(given = 8)
#'  }
#'
#' @import RCurl
#' @import XML
#' @import gdata
#' @export

rtematres.api <- function(task = "availableTasks", argument) {
	if(task == "availableTasks") {
		base_service_url = getURL("tematres.befdata.biow.uni-leipzig.de/vocab/services.php")
		base_service_info = xmlTreeParse(base_service_url, useInternalNodes=T)
		tasks_available = trim(xpathSApply(base_service_info, "//task", xmlValue)[-1])
		tasks_description = trim(xpathSApply(base_service_info, "//action", xmlValue))
		tasks_argument =  trim(xpathSApply(base_service_info, "//arg", xmlValue)[-1])
		return(data.frame(tasks_available, tasks_description, tasks_argument))
	}

	# select task
	task = match.arg(task, c("fetchVocabularyData", "suggest", "suggestDetails", "fetchTopTerms", "search", "fetchCode", "letter", "fetchTerm", "fetchAlt", "fetchDown", "fetchUp", "fetchRelated", "fetchNotes", "fetchDirectTerms", "fetchURI", "fetchTargetTerms", "fetchSourceTerms", "fetchTerms", "fetchRelatedTerms", "letter", "fetchLast", "fetchSimilar"))

	tasks_need_no_argument = c("fetchLast", "fetchTopTerms", "fetchVocabularyData")

	if(any(task == tasks_need_no_argument)) {
	  url_param_sep = "?"
	  param_sep = "&"
	  task_trigger_name = "task"
	  assignment = "="
	  service_url = paste0(rtematres.options("tematres_service_url"), url_param_sep, task_trigger_name, assignment, task)
	} else {
	  url_param_sep = "?"
	  param_sep = "&"
	  task_trigger_name = "task"
	  arg_trigger_name = "arg"
	  assignment = "="
	  service_url = paste0(rtematres.options("tematres_service_url"), url_param_sep, task_trigger_name, assignment, task, param_sep, arg_trigger_name, assignment, argument )
	}

	response = xmlTreeParse(service_url, useInternalNodes = T)


	# scheme to use
	base_list = list(id = "//term/term_id",
		term = "//term/string")

	notes_list = list(id = "//term/term_id",
		term = "//term/string",
		language = "//term/note_lang",
		description = "//term/note_text")

	fetchVocabularyData_list = list(author = "//author",
		title = "//title",
		language = "//lang",
		uri = "//uri",
		lastMod = "//lastMod",
		count_terms = "//cant_terms",
		status = "//status")

	suggest_list = list(term = "//result/term")

	sheme = switch(task, fetchVocabularyData = fetchVocabularyData_list,
				fetchTopTerms = base_list,
				fetchCode = base_list,
				search = base_list,
				suggest = suggest_list,
				suggestDetails = base_list,
				letter = base_list,
				fetchAlt = base_list,
				fetchTerm = base_list,
				fetchAlt = base_list,
				fetchDown = base_list,
				fetchUp = base_list,
				fetchRelated = base_list,
				fetchNotes = notes_list,
				fetchDirectTerms = base_list,
				fetchURI = base_list,
				fetchTargetTerms = base_list,
				fetchSourceTerms = base_list,
				fetchTerms = base_list,
				fetchRelatedTerms = base_list,
				fetchSimilar = base_list,
				fetchLast = base_list
				)
	if(is.list(sheme)) {
	  rapply(sheme, function(x) xmlNodesValue(path=x, doc=response), how="replace")
	} else {
	  return(response)
	}
}

