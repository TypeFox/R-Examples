
# x is name of a tcl variable holding list of character vectors
tclList2R <- function(x, convert = as.character) {
	.Tcl <- tcltk::.Tcl
	len <- as.integer(.Tcl(sprintf("llength $%s", x)))
	f <- function(i) convert(.Tcl(sprintf("lindex $%s %d", x, i)))
	lapply(seq(0, len-1), f)
}

# high performance strapply with hard coded FUN=c.  Guts in tcl.
strapplyc <- function(X, pattern, backref = NULL, ignore.case = FALSE, simplify = FALSE, USE.NAMES, engine) {
    if (missing(engine)) engine <- getOption("gsubfn.engine")
    if (missing(USE.NAMES)) USE.NAMES <- FALSE
    if (identical(engine, "R")) return(
		strapply(X = X, pattern = pattern, FUN = c, backref = backref, 
			ignore.case = ignore.case, simplify = simplify, 
			USE.NAMES = USE.NAMES, engine = engine)
	)
	.Tcl <- tcltk::.Tcl
	tcl <- tcltk::tcl
	tcl("set", "X", tcltk::as.tclObj(X))
	tcl("set", "pattern", pattern)
	tcl("set", "nocase", if (ignore.case) "-nocase" else "")
	if (missing(backref) || is.null(backref) || is.na(backref))  backref <- 999
	tcl("set", "backref", backref)
	.Tcl("set about [regexp -about -- $pattern]")
	.Tcl("set about [lindex $about 0]")
	.Tcl("if { min($about, $backref) <= 0 } { set mn 0 } else { set mn 1 }")
	.Tcl("set mx [expr min($about, abs($backref))]")
	s <- paste('set result {}
		 set k [expr $about + 1]
		 if { $about == 0 || $about <= -$backref} {
			 # this leg of the "if" returns everything from regexp so we 
			 # can avoid the extraction subloop of the "else" leg for speed
			 foreach item $X {
				# {*} is new feature in tcl 8.5 to add level of substitution
				set cmd [list regexp -all -inline {*}$nocase -- $pattern $item] 
				set res [{*}$cmd]
				lappend result $res
			}
		 } else {
			foreach item $X {
				# {*} is new feature in tcl 8.5 that adds level of substitution
				set cmd [list regexp -all -inline {*}$nocase -- $pattern $item] 
				set cmdout [{*}$cmd]
				set imin $mn
				set imax $mx
				set res {}
				while {$imax < [llength $cmdout]} {
					lappend res [lrange $cmdout $imin $imax]
					incr imin $k
					incr imax $k
				}
				lappend result [concat {*}$res]
			}
		}')
	.Tcl(s)
	out <- tclList2R("result")

    result <- sapply(out, identity, simplify = isTRUE(simplify), 
		USE.NAMES = USE.NAMES)
    if (is.logical(simplify)) result else {
		do.call(match.funfn(simplify), result)
	}
}
