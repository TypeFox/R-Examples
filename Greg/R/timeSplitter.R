#' A function for splitting a time according to time periods
#' 
#' If we have a violation of the cox proprtional hazards assumption we need to
#' split an individual's followup time into several.See vignette("timeSplitter", package="Greg")
#' for a detailed description.
#' 
#' \emph{Important note:} The time variables must have the same time unit. I.e. function can not dedu
#' if all variables are in years or if one happens to be in days.
#' 
#' @section The time_offset - details:
#' 
#' Both time_var and other variables will be adjusted by the time_offset, 
#' e.g. if we the time scale is in years and we want to skip the first 4 years 
#' we set the \code{time_offset = 4}. In the outputted dataset the smallest 
#' \code{time_var} will be 0. \emph{Note:} 0 will not be included as we 
#' generally want to look at those that survived the start date, e.g. if a 
#' patient dies on the 4-year mark we would not include him/her in our study.
#'  
#' @param data The dataset that you want to split according to the \code{time_var}
#'  option.
#' @param by The time period that you want to split the dataset by. The size of the variable
#'  must be in proportion to the the \code{time_var}.
#' @param event_var The event variable
#' @param event_start_status The start status of the event status, e.g. "Alive"
#' @param time_var The name of the main time variable in the dataset. This variable
#'  must be a numeric variable.
#' @param time_offset If you want to skip the initial years you can offset the 
#'  entire dataset by setting this variable. See detailed description below.
#' @param time_related_vars A dataset often contains other variabels that you want
#'  to update during the split, most commonly these are age or calendar time. 
#' @return \code{data.frame} with the split data. The starting time for each period
#'  is named \code{Start_time} and the ending time is called \code{Stop_time}. Note
#'  that the resulting event_var will now contain the time-splitted eventvar.
#'  
#' @example inst/examples/timeSplitter_example.R
#' @export
#' @importFrom Epi cal.yr Lexis splitLexis
timeSplitter <- function (data, by, time_var, 
                          event_var, 
                          event_start_status, 
                          time_related_vars, time_offset) {
  # Save the original order of the names for restoring at the end
  # Also save the labels if any
  org_names <- names(data)
  attr(org_names, "labels") <- sapply(data, label)
  # There is a bug in Lexis that doesn't allow variables of non-base type
  if (any(attr(org_names, "labels") != "")){
    for (i in 1:ncol(data)){
      if (inherits(data[[i]], "labelled")){
        class(data[[i]]) <- class(data[[i]])[class(data[[i]]) != "labelled"]
        attr(data[[i]], "label") <- NULL
      }
    }
  }
  
  if (missing(time_var)){
    if ("time" %in% names(data) &&
        is.numeric(data$time)){
      time_var <- "time"
    }else{
      stop("You need to provide the main time variable name.",
           " If not provided the variable may be named 'time' in the dataset")
    }
  } else if (length(time_var) > 1) {
    stop("You can only have one time variable in the dataset")
  } else if (!time_var %in% names(data)) {
    stop("Could not find the time variable among the variables in the dataset",
         ", i.e. '", time_var, "' not in '", paste(names(data), 
                                                   collapse = "', '"), "'")
  } else if (!is.numeric(data[[time_var]])){
    stop("The time variable must be numeric.")
  }
  
  if (!missing(time_offset)){
    if (!is.numeric(time_offset))
      stop("If you want the time to be offsetted by a certain amount of time",
           " then you need to provide a numeric time_offset. You've provided",
           " '", time_offset, "'")
    if (time_offset < 0)
      stop("The time_offset must be more than 0")
    
    data <- subset(data,
                   parse(text = 
                           paste(time_var, '>', time_offset)))
    data[[time_var]] <- data[[time_var]] - time_offset
    if (!missing(time_related_vars)){
      for (var in time_related_vars)
        data[[var]] <- data[[var]] - time_offset
    }
  }

  # Change name to temporary names so that we can revert after
  if (!missing(time_related_vars)){

    if (any(!time_related_vars %in% names(data)))
      stop("You have specified variables not in the dataset: ",
           "'", paste(time_related_vars[!time_related_vars %in% names(data)],
                 collapse = "', '"),
           "'")
    
    for (var in time_related_vars){
      if (is.numeric(data[[var]]) &&
          max(data[[var]])*10 < max(data[[time_var]]))
        warning("Your time variable seems much larger than the time related variable '", var, "'",
                " - Are you sure that they are both the same time unit, i.e. both are days/years?")
      data[[sprintf("u__temp_time__%s", var)]] <- data[[var]]
      data[[var]] <- NULL
    }
  }
  
  if (!event_var %in% names(data))
    stop("Could not find the event varible '", event_var ,"'",
         " among the dataset variabels: ",
         "'", paste(names(data), 
                    collapse = "', '"),
         "'")

  # Setup varibles for the Lexis object
  data <- cal.yr(data)
  # The cal.yr adds a "cal.yr" class that needs to be removed or we get
  # Error: Duration is not the same on all time scales
  rmv_class <- sapply(data, 
                      function(v)
                        "cal.yr" %in% class(v))
  if (any(rmv_class)){
    for (i in which(rmv_class))
      class(data[[i]]) <- class(data[[i]])[!class(data[[i]]) %in% "cal.yr"]
  }
  if (missing(event_start_status)){
    if (is.factor(data[[event_var]])){
      event_start_status <- levels(data[[event_var]])[1]
    }else if(is.logical(data[[event_var]])){
      event_start_status <- FALSE
    }else if(all(data[[event_var]] %in% 0:1)){
      event_start_status <- 0
    }else{
      stop("No initial status provided for the event variable.",
           " This is usually the status corresponding to the",
           " alive status. The function tries to deduce this",
           " from the factors if available alternatively from",
           " the 0/FALSE if the event variable happens to be",
           " boolean or numeric with only 0 and 1 values.")
    }
  }
  
  if (is.character(data[[event_var]]))
    data[[event_var]] <- factor(data[[event_var]])
  
  if (!event_start_status %in% unique(data[[event_var]])){
    stop("The event status should be among the levels in the event variable.",
         " Currently '", event_start_status, "' is not found among",
         " '", paste(unique(data[[event_var]]),
                     collapse = "', '"),
         "'.")
  }
  if (is.factor(data[[event_var]])){
    data$Entry <- factor(rep(event_start_status, nrow(data)),
                         levels = levels(data[[event_var]]))
  }else{
    data$Entry <- event_start_status
  }

  if (is.null(data[["Start_time"]]))
    data$Start_time <- 0
  
  # Use strange names in order to avoid conflict
  exit_list <- sprintf("list(u___timeband__ = `%s`)", 
                       time_var)
  if (!missing(time_related_vars)){
    for (var in time_related_vars){
      exit_list <- gsub(")$",
                        sprintf(", `%s` = `u__temp_time__%s` + `%s`)", var, var, time_var),
                        exit_list)
    }
  }
  # This is necessary as the Lexis uses substitute 
  # with the data as a co-argument and thus very
  # complex to handle with parse()/deparse()/expression()
  lxs_data <- 
    paste0("Lexis(entry = list(u___timeband__ = Start_time),
                    exit = ", exit_list, ",
                    exit.status = ", event_var, ",
                    entry.status = Entry,
                    data = data)") %>% 
    parse(text = .) %>% 
    eval
  
  if (length(by) == 1)
    by <- seq(0, ceiling(max(lxs_data$lex.dur)), by = by)
  
  spl <- splitLexis(lxs_data,
                    breaks = by,
                    time.scale = "u___timeband__")
  
  spl$Start_time <- spl$u___timeband__
  spl$Stop_time <- spl$u___timeband__ + spl$lex.dur

  # Remove all extra Lexis variables since these will only be confusing
  spl[[event_var]] <- spl$lex.Xst
  spl$lex.Cst <- NULL
  spl$lex.Xst <- NULL
  spl$lex.id <- NULL
  spl$u___timeband__ <- NULL
  spl$lex.dur <- NULL
  spl$Entry <- NULL
  # Remove the time indicator as this shouldn't be used later on
  spl[[time_var]] <- NULL
  # Drop temp variables
  if (!missing(time_related_vars)){
    for (var in time_related_vars){
      spl[[sprintf("u__temp_time__%s", var)]] <- NULL
    }
  }
  
  re_order <- sapply(names(spl), function(x) {
    ordr <- which(x == org_names)
    if (length(ordr) == 0)
      ordr <- -1
    return(ordr)
  })
  for (i in 1:length(re_order)){
    if (re_order[i] == -1){
      re_order[i] <- max(re_order) + 1
    }
  }
  if (any(attr(org_names, "labels") != "")){
    lbl <- attr(org_names, "labels")
    for(var in names(lbl)[lbl != ""]){
      if (var %in% names(spl)){
        label(spl[[var]]) <- lbl[var]
      }
    }
  }
  
  return(spl[,order(re_order)])
}

globalVariables(c(".", "label<-"))