#
#	Rparallel_setEnable_pkg.R
#Mon Mar 25 17:06:14 CET 2013

# enable/disable parallelize.dynamic functionality during package use

parallelize_setEnable = function(state) {
	sourceFile = if (!state) {
		system.file('Rscripts/Rparallel_functions_std.R', package = 'parallelize.dynamic')
	} else {
		system.file('Rscripts/Rparallel_functions_parallel.R', package = 'parallelize.dynamic')
	}
	Log(sourceFile, 1);
	source(sourceFile);
}
