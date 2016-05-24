# Functions to quickly query certain information from the solution

#' Get summary information from all databases
#'
#' Get the list of phases, samples, timeslices and bands  that are available in each database.
#'
#' @inheritParams query_master
#'
#' @family special queries
#' 
#' @export
query_phase <- function(db) {
  sql <- "SELECT DISTINCT phase_id, period_type_id AS is_summary FROM key"
  query_sql(db, sql) %>%
    add_phase_names %>%
    select(scenario, position, phase_id, phase, is_summary) %>%
    arrange(position, phase_id, is_summary)
}

#' @rdname query_phase
#' @export
query_sample <- function(db) {
  sql <- "SELECT DISTINCT phase_id, period_type_id AS is_summary, sample FROM key"
  query_sql(db, sql) %>%
    add_phase_names %>%
    select(scenario, position, phase_id, phase, is_summary, sample) %>%
    arrange(position, phase_id, is_summary)
}

#' @rdname query_phase
#' @export
query_timeslice <- function(db) {
  sql <- "SELECT DISTINCT phase_id, period_type_id AS is_summary, timeslice FROM key"
  query_sql(db, sql) %>%
    add_phase_names %>%
    select(scenario, position, phase_id, phase, is_summary, timeslice) %>%
    arrange(position, phase_id, is_summary)
}

#' @rdname query_phase
#' @export
query_band <- function(db) {
  sql <- "SELECT DISTINCT phase_id, period_type_id AS is_summary, band FROM key"
  query_sql(db, sql) %>%
    add_phase_names %>%
    select(scenario, position, phase_id, phase, is_summary, band) %>%
    arrange(position, phase_id, is_summary)
}

#' @rdname query_class_member
#' @export
query_class <- function(db) {
  sql <- "SELECT DISTINCT class_group, class FROM key"
  query_sql(db, sql) %>%
    select(-filename) %>%
    arrange(position, class_group, class)
}

#' Get list of objects from all databases
#'
#' Get the list of objects and classes in each database. Shortcuts for generators, regions and zones
#' are provided for convenience.
#'
#' @inheritParams query_master
#' @param class Type of class to query
#'
#' @family special queries
#' 
#' @export
query_class_member <- function(db, class) {
  sql <- sprintf("SELECT DISTINCT name, parent, region, zone FROM key WHERE class = '%s'", class)
  query_sql(db, sql) %>%
    select(-filename) %>%
    arrange(position, name)
}

#' @rdname query_class_member
#' @export
query_generator <- function(db) {
  query_class_member(db, "Generator") %>%
    filter(parent == "System") %>%
    select(-parent)
}

#' @rdname query_class_member
#' @export
query_region <- function(db) {
  query_class_member(db, "Region") %>%
    select(-region, -zone)
}

#' @rdname query_class_member
#' @export
query_zone <- function(db) {
  query_class_member(db, "Zone") %>%
    select(-region, -zone)
}

#' Get time spans from all databases
#'
#' Get the time limits and time step lengths for each simulation phase.
#'
#' @inheritParams query_master
#'
#' @family special queries
#' 
#' @export
query_time <- function(db) {
  sql <- "SELECT phase_id, min(time) start, max(time) end, count(time) count
          FROM time GROUP BY phase_id"
  query_sql(db, sql) %>%
    add_phase_names %>%
    mutate(start = lubridate::ymd_hms(start, quiet = TRUE),
           end = lubridate::ymd_hms(end, quiet = TRUE),
           timestep = difftime(end, start, units = "mins") / (count - 1)) %>%
    select(scenario, position, phase_id, phase, start, end, count, timestep) %>%
    arrange(position, phase_id)
}

# Shortcut to add phase names to a result
add_phase_names <- function(x) {
  phases <- c("LT", "PASA", "MT", "ST")
  phases.df <- data.frame(phase_id = 1:4, phase = factor(phases, levels = phases))
  x %>% inner_join(phases.df, by = "phase_id")
}
