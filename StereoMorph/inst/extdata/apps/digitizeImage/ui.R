require(bezier)
require(rjson)

# LOAD INITIAL PARAMETERS
source(paste0('initial_parameters.R'))

shinyUI(
	bootstrapPage(
		includeHTML("digitize_image.html"),
		title='StereoMorph | Digitizing App', theme = 'digitize_image.css',
		textInput("text_input", label = "", value = json_string),
		textOutput("text_output")
	)
)