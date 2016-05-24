
context(
	"error functions on report have to be shown to be fully functional
	when given reasonable inputs."
)

test_that(
	"check that each of the expected fields are being 
	displayed in the output (won't happen if function is broken", {

		pcall <- 'this is a call'
		traits <- c("trait one", "trait two")
		invalid <- "white-elephant"
		inputs <- list(value = 10000, trait = "string")
		actual <- "actual-value"
		value  <- "this is a value"
		name <- "cat"
		error <- list(message = "I am an error message")

		expect_error(report$missing_traits(pcall), pcall)

		expect_error(report$missing_value(pcall), pcall)

		expect_error(
			report$traits_not_character(pcall, traits),
			pcall)

		expect_error(
			report$invalid_traits(pcall, invalid),
			pcall)

		expect_error(
			report$non_boolean(pcall, inputs, actual),
			pcall)

		expect_error(
			report$no_match(pcall, value, traits),
			pcall)

		expect_error(
			report$error_encountered(pcall, error, inputs),
			pcall)

		expect_error(
			report$warning_encountered(pcall, error, inputs),
			pcall)
		expect_warning(
			report$trait_overridden(pcall, name),
			pcall)
	}
)
