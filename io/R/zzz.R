.onLoad <- function(libname, pkgname) {

	# plotting options
	opts <- list(
		width = 5,
		height = 5,
		res = 300,
		units = "in"
	);
	options(plot=opts);

	if (!interactive() || identical(getOption("device"), pdf)) {
		# Either the current session is not interactive or a graphical
		# screen cannot be open; therefore, print directly to inferred device
		# (R sets the device to pdf if no graphical screen can be opened)
		options(plot.device=NULL);
	} else {
		options(plot.device=NA);
	}

}
