# Assertions used for checking function arguments

assert_param_presence <- function(param_name, given_names_vector) {
  # precondition: given_names vector is a vector
  if (!(param_name %in% given_names_vector))
    stop(paste(param_name, "is a required parameter"))
}

assert_param_class <- function(param, expected_class) {
  if (class(param) != expected_class)
    stop(paste("Param class should be:", expected_class))
}

assert_param_subset <- function(given_names_vector, allowed_names_vector,
                                error_message = "Invalid param:") {
  for (param_name in given_names_vector)
    if (is.null(param_name) || !(param_name %in% allowed_names_vector))
      stop(paste(error_message, param_name))
}

assert_param_not_null <- function(given_param) {
  if (is.null(given_param))
    stop("Given param cannot be NULL")
}

assert_param_integer <- function(given_param) {
  # precondition: given_param is a single element
  if (class(given_param) != 'numeric' || !(given_param%%1 == 0))
    stop(paste("Given param is not an integer:"), given_param)
}

assert_param_single <- function(given_param) {
  # precondition: given_param is not NULL
  if (length(given_param) != 1)
    stop(paste("Length of given param is not 1:", given_param))
}

assert_param_range <- function(given_param, min, max, param_name) {
  # precondition: given_param is an integer
  if (given_param < min || given_param > max)
    stop(paste("The ",
               param_name,
               " has to be an integer in range ",
               min, "-", max, sep = ""))
}

assert_param_nrow <- function(data_frame, minimum = NULL, maximum = NULL) {
  if (!is.null(minimum) && nrow(data_frame) < minimum)
    stop(paste("The number of rows in the data frame is below the minimum of", minimum))
  # maximum checking not used anywhere currently
  if (!is.null(maximum) && nrow(data_frame) > maximum)
    stop(paste("The number of rows in the data frame is above the maximum of", maximum))
}
