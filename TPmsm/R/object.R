AddClasses <- function(object, value) {
	invisible( .Call(Rf_SetClass, object, c(class(object)[1], value), PACKAGE="TPmsm") )
}

RemoveClasses <- function(object) {
	invisible( .Call(Rf_SetClass, object, class(object)[1], PACKAGE="TPmsm") )
}
