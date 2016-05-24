`plotObservedEffects` <-
function(x, ...) {
	if (!is.null(attr(x, "class"))){
		UseMethod("plotObservedEffects")
	}
	else {
		plotObservedEffectsDefault(x, ...)
	}
}

