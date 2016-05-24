#
#	Rfunctions.R
#Tue 14 Aug 2007 01:39:42 PM CEST 

#
#	<§> abstract data functions
#

inverse = function(f, interval = c(-Inf, Inf)) {
	function(y, ...) {
		optimize(function(x, ...){ (y - f(x, ...))^2 }, interval = interval, ...)$minimum
	}
}

#
#	<p> meta functions
#

callWithArgs = function(fctName, args) {
	#arguments = paste(sapply(names(args), function(n)sprintf("%s = %s", n, args[[n]])), collapse = ", ");
	fhead = sprintf("%s(%s)", fctName, paste(names(args), collapse = ", "));
	eval(parse(text = fhead))
}

.do.call = function(f, args, restrictArgs = T) {
	if (restrictArgs) {
		fargs = names(as.list(args(f)));
		fargs = fargs[fargs != ''];
		args = args[which.indeces(fargs, names(args))];
	}
	do.call(f, args)
}

#
#	<p> benchmarking
#

benchmark.timed = function(.f, ..., N__ = 1e1) {
	t0 = Sys.time();
	for (i in 1:N__) {
		r = .f(...);
	}
	t1 = Sys.time();
	r = list(time = (t1 - t0)/N__, lastResult = r, t0 = t0, t1 = t1);
	print(r$time);
	print(r$t0);
	print(r$t1);
	r
}
