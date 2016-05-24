##
## Classes with Color-, Interaction- & Display-Settings
########################################################
## TODO: Eventually also save tk2d settings here

setClass(
		Class = "ColorSettings",
		representation = representation(
				background = "character",
				bullet = "character",
				bulletActive = "character",
				nodes = "character",
				nodesActive = "character",
				adjNodes = "character",
				adjNodesActive = "character",
				notVisitedEdge = "character",
				visitedEdge = "character",
				edgeActive = "character",
				labels = "character",
				labelsActive = "character",
				adjLabels = "character",
				adjLabelsActive = "character",
				path = "character"
		),
		prototype=list(
				background = "white",
				bullet = "#fdfd96",           ## Pastel Yellow
				bulletActive = "#ff6961",     ## Pastel Red
				nodes = "#444444",            ## Gray
				nodesActive = "#03C03C",      ## Pastel Dark Green
				adjNodes = "#FFB347",         ## Pastel Orange
				adjNodesActive =  "#ff6961",  ## Pastel Red
				notVisitedEdge = "#444444",   ## Gray
				visitedEdge ="#CFCFC4",       ## Pastel Gray
				edgeActive = "#03C03C",       ## Pastel Dark Green
				labels = "#444444",				    ## Gray
				labelsActive = "#03C03C",     ## Pastel Dark Green
				adjLabels = "#ff8f00",        ## Princeton Orange
				adjLabelsActive ="#ff6961",   ## Pastel Red
				path = "#c23b22"              ## Pastel Dark Red
		)
)

setClass(
		Class = "InteractionSettings",
		representation = representation(
				NSteps = "numeric",
				animationTime = "numeric",
				dragSelectRadius = "numeric",
				labelDistRadius = "numeric"
		),
		prototype = list(
				NSteps = 50,
				animationTime = 0.1,
				dragSelectRadius = 15,
				labelDistRadius = 30
		)
)

setClass(
		Class = "DisplaySettings",
		representation = representation(
				bulletRadius = "numeric",
				nodeRadius = "numeric",
				lineWidth = "numeric",
				highlightedLineWidth = "numeric"
		),
		prototype = list(
				bulletRadius = 15,
				nodeRadius = 10,
				lineWidth = 2,
				highlightedLineWidth = 4
		)
)

setClass(
		Class = "Tk2dDisplay",
		representation = representation(
				bg = "character",
				brush_colors = "character",
				brush_color = "character",
				linked = "logical"
		),
		prototype = list(
				bg = "white",
				## RColorBrewer qualitative, 9 classes, color scheme: Set1
				brush_colors = c('#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999', '#E41A1C'),
				brush_color = "magenta",
				linked = TRUE
		)
)


setClass(
		Class = "NG_Settings",
		representation = representation(
				color = "ColorSettings",
				interaction = "InteractionSettings",
				display = "DisplaySettings",
				tk2d = "Tk2dDisplay"
		),
		prototype = list(
				color = new('ColorSettings'),
				interaction = new("InteractionSettings"),
				display = new("DisplaySettings"),
				tk2d = new("Tk2dDisplay")
		)
)
