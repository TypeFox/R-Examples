TransMethod <- function(object, s, t, ...) {
	UseMethod("TransMethod")
}

TransBoot <- function(object, s, t, ...) {
	UseMethod("TransBoot")
}

TPBoot <- function(object, UT, nboot, ...) {
	UseMethod("TPBoot")
}

TPCBoot <- function(object, UT, UX, nboot, ...) {
	UseMethod("TPCBoot")
}

TransPROB <- function(object, UT, ...) {
	UseMethod("TransPROB")
}

toTPmsm <- function(lst, UT, s, t, statenames) {
	UseMethod("toTPmsm")
}

BtoTPmsm <- function(lst, UT, s, t, statenames, nboot, conflevel, methodboot) {
	UseMethod("BtoTPmsm")
}

TransMatrix <- function(x) {
	UseMethod("TransMatrix")
}
