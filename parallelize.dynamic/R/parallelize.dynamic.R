#
#	parrallelize.dynamic.R
#

.onLoad = function(libname, pkgname) {
	Log.setLevel(4);
	perlPath = system.file('Perl', package = pkgname);
	Sys.setenv(PERL5LIB = sprintf('%s:%s', perlPath, Sys.getenv('PERL5LIB')));
	Sys.setenv(PATH = sprintf('%s:%s', perlPath, Sys.getenv('PATH')));
	#parallelize_setEnable(F);
}
