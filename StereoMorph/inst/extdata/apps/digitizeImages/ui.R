require(bezier)
require(rjson)

# SESSION PARAMETERS SENT DIRECTLY TO BROWSER
text_input <- readLines('session_parameters.txt')

shinyUI(
	bootstrapPage(
		includeHTML("digitize_image.html"),
		title='StereoMorph | Digitizing App', theme = 'digitize_image.css',
		textInput("text_input", label = "", value = text_input),
		textOutput("text_output")
	)
)