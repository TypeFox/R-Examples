
report <- list(
	# an object (well, list...) containing 
	# functions that report various errors and warnings
	missing_traits = function (pcall) {

		stopf (
			"%s: the parameter 'traits' was missing but is required\n", 
			pcall)
	},
	missing_value = function (pcall) {
		stopf (
			"%s: the parameter 'value' was missing but is required\n", 
			pcall)
	},
	traits_not_character = function (pcall, traits) {
		stopf (
			"%s: the parameter 'traits' must be a character vector\n
			actual class was %s
			", 
			pcall, paste0(class(traits), collapse = ", "))
	},
	inputs_not_list = function (pcall, inputs) {
		stopf(
			"%s: the parameter 'traits' must be a list\n
			actual class was %s
			", 
			pcall, paste0(class(inputs), collapse = ", "))
	},
	invalid_traits = function (pcall, invalid) {
		stopf(
			"%s: unrecognised trait(s): (%s)\n", 
			pcall, 
			paste0(invalid, collapse = ', '))		
	},
	non_boolean = 
		function (pcall, inputs, actual) {

			text <- '%s:
				the value %s returned a non true/false value when tested for the trait %s:\n
				actual value was %s\n'

			readable <- list(
				value = deparse_to_string(inputs$value),
				actual = deparse_to_string(actual)
			)

			stopf(text,
				pcall, 
				readable$value, inputs$trait,
				readable$actual)
		},
	no_match =
		function (pcall, value, traits) {
			
			text <- "%s: 
				the value %s didn't match any of the following compound traits:
				%s\n"

			and_collapse <- function (val) {
				paste0(val, collapse = ' and ')
			}
			or_collapse <- function (val) {
				paste0(unlist(val), collapse = ', or ')
			}

			readable <- list(
				value = deparse_to_string(value),
				expected = or_collapse(sapply(traits, and_collapse))
			)

			stopf(text,
				pcall,
				readable$value,
				readable$expected)
		},
	error_encountered = 
		function (pcall, error, inputs) {

			text <- '%s:\n
			an error was encountered while testing the value %s for the the trait "%s":\n
			%s\n'

			readable <- list(
				value = deparse_to_string(inputs$value)
			)

			stopf(text,
				pcall, readable$value,
				inputs$trait, error$message)
		},
	warning_encountered =
		function (pcall, warning, inputs) {

			text <- '%s:\n
			a warning was encountered while testing the value %s for the the trait "%s":\n
			%s\n'

			readable <- list(
				value = deparse_to_string(inputs$value)
			)

			stopf(text,
				pcall, readable$value,
				inputs$trait, warning$message)
		},
	trait_overridden = 
		function (pcall, name) {

			text <- "%s:\n the trait '%s' already exists: overwriting.\n"

			warningf(text,
				pcall, name)
		}
)
