openv <- new.env()

assign("N", 1e3, envir = openv)

openv.setN <- function(x) {
	assign("N", x, envir = openv)
}