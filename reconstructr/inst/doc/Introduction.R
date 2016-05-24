## ---- eval=FALSE---------------------------------------------------------
#  library(reconstructr)
#  actions_by_user <- list(c(1417330230, 1417330250, 1417330295, 1417330324, 1417330416),
#                          c(1417401697, 1417401741, 1417401751, 1417403263))
#  sessions <- reconstruct_sessions(timestamps = actions_by_user, threshold = 1800)

## ---- eval=FALSE---------------------------------------------------------
#  length_of_sessions <- session_length(sessions, padding_value = 430, preserve_single_events = FALSE, strip_last = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  events_per_session <- session_events(sessions)

## ---- eval=FALSE---------------------------------------------------------
#  bounce_rate <- session_events(sessions, decimal_places = 2)

## ---- eval=FALSE---------------------------------------------------------
#  time_on_page <- event_time(sessions, as_vector = FALSE, fun, ...)

## ---- eval=FALSE---------------------------------------------------------
#  average_time_on_page <- event_time(sessions, as_vector = FALSE, mean)

