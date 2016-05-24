#' @importFrom R6 R6Class
simulator_class <- R6Class("simulator",
  private = list(
    name = "TEMPLATE",
    priority = 50
  ),
  public = list(
    simulate = function() stop("virtual method"),
    get_cmd = function() stop("virtual method"),
    get_name = function() private$name,
    get_priority = function() private$priority,
    get_info = function() c(name = private$name),
    initialize = function(priority) {
      assert_that(is.number(priority))
      private$priority <- priority
    }
  )
)

is_simulator <- function(simulator) inherits(simulator, "simulator")


# Keep a user modifiable list of available simulation programs in a private
# environment
if (!exists("simulators")) simulators <- new.env()

register_simulator <- function(simulator) {
  assert_that(is_simulator(simulator))
  simulators[[simulator$get_name()]] <- simulator
}

get_simulator <- function(name) {
  simulators[[name]]
}


fill_cmd_template <- function(template, model, parameters,
                              locus_group, eval_pars = TRUE) {

  locus_length <- get_locus_length(model, group = locus_group)
  total_locus_number <- get_locus_number(model, locus_group, TRUE)

  if (has_variation(model)) {
    locus_number <- rep(1, total_locus_number)
  } else {
    locus_number <- total_locus_number
  }
  locus_id <- seq(along = locus_number)

  args <- vapply(locus_id, function(l_id) {
    tmp_env <- create_par_env(model, parameters,
                              locus_length = locus_length,
                              locus_id = l_id,
                              locus_number = total_locus_number,
                              for_cmd = !eval_pars)

    paste(eval(parse(text = template), envir = tmp_env), collapse = " ")
  }, character(1L))

  sim_cmds <- data.frame(locus_number = locus_number,
                         command = args,
                         stringsAsFactors = FALSE)

  reduce_sim_commands(sim_cmds)
}


reduce_sim_commands <- function(sim_commands) {
  if (nrow(sim_commands) == 1) return(sim_commands)
  grouped_commands <- unique(sim_commands[ , 2])
  if (length(grouped_commands) == nrow(sim_commands)) return(sim_commands)
  grouped_locus_number <- vapply(grouped_commands, function(cmd) {
    sum(sim_commands[sim_commands[ , 2] == cmd, 1])
  }, numeric(1)) #nolint
  data.frame(locus_number = grouped_locus_number,
             command = grouped_commands,
             stringsAsFactors = FALSE)
}


#' Returns the available simulators
#'
#' This functions returns the usable simulators
#'
#' @export
#' @examples
#' list_simulators()
list_simulators <- function() {
  simulators <- do.call(rbind, lapply(ls(simulators), function(simulator) {
    info <- get_simulator(simulator)$get_info()
    name <- info[["name"]]
    info <- info[-1]
    pars <- paste(names(info), ":", info, collapse = ", ")
    data.frame(name = name,
               priority = get_simulator(simulator)$get_priority(),
               info = pars)
  }))
  simulators[order(simulators$priority, decreasing = TRUE), ]
}
