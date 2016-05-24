# nlarDialog <- function(series){
#     require(tcltk)||stop("tcltk package required for displaying the GUI")
# 
#     models <- availableModels()
# 	vModel <- tclVar(models[1])
# 
# 	frMain <- Frame(opts=list(anchor="w"))
# 	frTop <- Frame()
# 	add(frTop, frMain)
# 	add(frMain, Widget(opts=list(type="label", text="Model")))
# 
#     for(i in 1:length(models))
#         add(frMain, Widget(opts = list(type = "radiobutton", text=models[i], value=models[i], variable=vModel)))
# 
#     onNext <- function() {
#         model <- as.character(tclObj(vModel))
#         tkdestroy(frRoot$tkvar)
# 	    fun <- get(paste("showDialog",model,sep="."))
# 	    fun(series)
# 	    return(invisible(NULL))
# 	}
# 
#     onCancel <- function()
#         tkdestroy(frRoot$tkvar)
# 
# 	frBottom <- makeButtonsFrame(list(Next=onNext, Cancel=onCancel))
# 	frRoot <- Frame()
# 	add(frRoot, frTop, frBottom)
# 	buildDialog("nlar model fitting", frRoot)
# 	cat("\n")
#     invisible(NULL)
# 
# }
# 
# nlar.struct.Frame <- function(vM, vD, vSteps) {
# 	frMain <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4))
# 	add(frMain,
# 		Widget(opts=list(type="label", text="dimension")), 
# 		Widget(opts=list(type="spinbox", from=1, to=1000, increment=1, textvariable=vM, width=4)),
# 		Widget(opts=list(type="label", text="time delay")),
# 		Widget(opts=list(type="spinbox", from=1, to=1000, increment=1, textvariable=vD, width=4)),
# 		Widget(opts=list(type="label", text="forecast steps")),
# 		Widget(opts=list(type="spinbox", from=1, to=1000, increment=1, textvariable=vSteps, width=4))
# 	)
# 	return(frMain)
# }
# 
# namedEntry <- function(name, var, width=2, ...) {
# 	frMain <- Frame()
# 	add(frMain,
# 		Widget(opts=list(type="label",text=name)),
# 		Widget(opts=list(type="entry",textvariable=var, width=width, ...))
# 	)
# 	return(frMain)
# }
# 
# namedSpinbox <- function(name, var, from, to, increment=1, width=4, ...) {
# 	frMain <- Frame()
# 	add(frMain,
# 		Widget(opts=list(type="label",text=name)),
# 		Widget(opts=list(type="spinbox",from=from, to=to, increment=increment, 
# 			textvariable=var, width=width, ...))
# 	)
# 	return(frMain)
# }
# 
# outputDialog <- function(oC) {
#     series <- oC$series
#     m <- oC$m
#     d <- oC$d
#     steps <- oC$steps
#     model <- oC$model
#     cl <- oC$cl #nlar function call
#     res <- eval(cl)
#     res$call <- NULL
#     assign("nlarModel", res, envir=.GlobalEnv)
#     return( invisible( NULL ) )
# }
