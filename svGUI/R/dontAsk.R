dontAsk <- function (gui = .GUI)
	return(!interactive() || !guiAsk(gui))
