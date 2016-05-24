digitizeImages <- function(image.file = image.file, shapes.file = NULL, 
	landmarks.file=NULL, control.points.file=NULL, curve.points.file=NULL, 
	cal.file = NULL, landmarks.ref = NULL, curves.ref = NULL, image.id=NULL, 
	landmark.color.blur = 'blue', landmark.color.focus = 'green', curve.color.blur = 'purple', 
	control.point.color.blur = 'purple', control.point.color.focus = 'red', landmark.radius = 4, 
	control.point.radius = 4, marker.stroke.width = 2,
	app.dir=NULL){

	# CHECK AND PROCESS INPUTS, ENSURE IMAGE.FILE AND SHAPES.FILE ARE MATRICES, WITH 
	#	COLUMNS CORRESPONDING TO VIEWS, SAVE INTO SESSION PARAMETERS LIST
	session_params <- process_digitize_images_input(image.file=image.file, 
		shapes.file=shapes.file, landmarks.file=landmarks.file, 
		control.points.file=control.points.file, curve.points.file=curve.points.file, 
		cal.file=cal.file, landmarks.ref=landmarks.ref, curves.ref=curves.ref)

	# GET STEREOMORPH SHINY APP DIRECTORY
	if(is.null(app.dir)){
		app_dir <- paste0(path.package("StereoMorph"), '/extdata/apps/digitizeImages')
	}else{
		app_dir <- app.dir
	}

	# REMOVE ANY IMAGE FILES IN WWW IMG FOLDER
	#if(length(list.files(paste0(app_dir, '/www/img/'))) > 0)
	#	file.remove(paste0(app_dir, '/www/img/', list.files(paste0(app_dir, '/www/img/'))))

	# SET IMAGE IDS
	if(!is.null(image.id)) warning("'image.id' is no longer supported. Images are labeled according to their file name. 'image.id' input has no effect.")

	# ADD ADDITIONAL PARAMETERS
	session_params$app_dir <- app_dir
	session_params$prev_wd <- getwd()
	session_params$landmark_color_blur <- landmark.color.blur
	session_params$landmark_color_focus <- landmark.color.focus
	session_params$curve_color_blur <- curve.color.blur
	session_params$control_point_color_blur <- control.point.color.blur
	session_params$control_point_color_focus <- control.point.color.focus
	session_params$landmark_radius <- landmark.radius
	session_params$control_point_radius <- control.point.radius
	session_params$marker_stroke_width <- marker.stroke.width

	# SAVE SESSION PARAMETERS TO JSON STRING FOR SERVER.R TO READ
	write(x=listToJSONStr(session_params), file=paste0(app_dir, "/session_parameters.txt"))
#	list2XML4R(list=session_params, file=paste0(app_dir, "/session_parameters.txt"))

	# COPY HTML FILE WITHOUT IMAGE TAG
	#file.copy(paste0(app_dir, "/digitize_image_pre.html"), paste0(app_dir, "/digitize_image.html"), overwrite=TRUE)

	# ADD IMAGE TAG TO HTML DOCUMENT (TO LOAD IMAGE AND GET SIZE)
	#img_tag <- paste0('\n<img style="display:none;" id="img1" src="img/', gsub(" ", "_", img_names[img_num]),'" ></img>')
	#write(img_tag, file=paste0(app_dir, "/digitize_image.html"), append=TRUE)

	# INITIATE SHINY APP
	#run_app <- runApp(app_dir, port = NULL, host = "127.0.0.1", launch.browser = TRUE, display.mode = "auto")
	run_app <- runApp(app_dir)

	# REMOVE ANY IMAGE FILES IN IMG FOLDER
	if(length(list.files(paste0(app_dir, '/www/img/'))) > 0)
		file.remove(paste0(app_dir, '/www/img/', list.files(paste0(app_dir, '/www/img/'))))

	return(NULL)
}