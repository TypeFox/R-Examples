#' @name HydeNetSummaries
#' @title HydeNet Summary Objects
#' 
#' @description Summaries of \code{HydeNetwork}, compiled network, and
#'   compiled decision network objects.
#'   
#' @param object A \code{HydeNet} object to be summarized
#' @param ... Additional arguments.
#' 
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @method summary HydeNetwork
#' @export

summary.HydeNetwork <- function(object, ...){
  
  decision_nodes <- names(object$nodeDecision)[sapply(object$nodeDecision, identity)]
  utility_nodes <- names(object$nodeUtility)[sapply(object$nodeUtility, identity)]
  deterministic_nodes <- names(object$nodeType)[sapply(object$nodeType, 
                                                        function(x) x == "determ")]
  
  random_nodes <- object$nodes[!object$nodes %in% c(decision_nodes, 
                                          utility_nodes,
                                          deterministic_nodes)]
  cat("Decision Nodes: \n",
      decision_node_summary(object, decision_nodes),
      "\n\n",
      "Utility Nodes: \n",
      utility_node_summary(object, utility_nodes),
      "\n\n",
      "Deterministic Nodes: \n",
      utility_node_summary(object, deterministic_nodes),
      "\n\n",
      "Random Nodes: \n",
      random_node_summary(object, random_nodes),
      sep = ""
  )
}

decision_node_summary <- function(object, nodes)
{
  name_summary <- summarise_node_name(nodes)
  parent_summary <- summarise_parents(object, nodes, max(nchar(name_summary)))
  policy_summary <- summarise_policy(object, nodes, 
                                     max(nchar(name_summary)), 
                                     max(nchar(parent_summary)))
  
  paste0(name_summary, parent_summary, policy_summary, collapse = "\n")
}

utility_node_summary <- function(object, nodes)
{
  name_summary <- summarise_node_name(nodes)
  parent_summary <- summarise_parents(object, nodes, max(nchar(name_summary)), end_sep = "")
  
  paste0(name_summary, parent_summary, collapse = "\n")
}

random_node_summary <- function(object, nodes)
{
  name_summary <- summarise_node_name(nodes)
  parent_summary <- summarise_parents(object, nodes, max(nchar(name_summary)))
  type_summary <- summarise_type(object, nodes)
  
  paste0(name_summary, parent_summary, type_summary, collapse = "\n")
}


summarise_node_name <- function(nodes, max.width = 20)
{
  max.width <- min(c(max.width, 
                     max(nchar(nodes)) + 3))
  
  ifelse(nchar(nodes) > (max.width - 2),
         paste0(substr(nodes, 1, 14), "...  |  "),
         paste0(stringr::str_pad(nodes, max.width - 2, "right"), "  |  "))
}

summarise_parents <- function(object, nodes, name_width, end_sep = "  |  ")
{
  max.width <- floor((getOption("width") - name_width) / 2)
  parents <- vapply(object$parents[nodes],
                    paste0,
                    character(1),
                    collapse = ", ")
  nparents <- vapply(object$parents[nodes],
                     length,
                     numeric(1))
  
  parents <- 
    ifelse(nchar(parents) > (max.width - 2),
           ifelse(nparents == 1,
                  "1 parent",
                  paste0(nparents, " parents  ")),
           parents)
  
  paste0(stringr::str_pad(parents, max(nchar(parents)), "right"), end_sep)
}

summarise_policy <- function(object, nodes, name_width, parent_width)
{
  max.width <- getOption("width") - name_width - parent_width
  decision_policy <-
    vapply(object$nodePolicyValues[nodes],
           paste0,
           character(1),
           collapse = ", ")
  
  decision_policy <- ifelse(decision_policy == "",
                            "(no policies defined)",
                            decision_policy)
  
  ifelse(nchar(decision_policy) > max.width,
         paste0(substr(decision_policy, 1, max.width - 3), "..."),
         decision_policy)
}

summarise_type <- function(object, nodes)
{
  unlist(object$nodeType[nodes])
}
