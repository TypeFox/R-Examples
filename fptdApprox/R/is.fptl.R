is.fptl <-
function (obj) 
{
    if (inherits(obj, "fptl") & is.list(obj) & (length(obj) == 
        2)) 
        if (all(unlist(lapply(obj, is.numeric))) & identical(names(obj), 
            c("x", "y")) & all(is.element(names(attributes(obj)), 
            c("names", "dp", "Call", "vars", "class")))) 
            if ((length(obj$x) == length(obj$y)) & is.diffproc(attr(obj, 
                "dp")) & is.call(attr(obj, "Call"))) {
		    if (!is.null(attr(obj, "vars"))) if (!is.list(attr(obj, "vars"))) return(FALSE)
                args <- as.list(attr(obj, "Call"))

                if (any(!is.element(c("dp", "t0", "T", "x0", 
                  "S"), names(args)))) 
                  return(FALSE)
                if (any(!(unlist(lapply(args[3:5], is.numeric))))) 
                  return(FALSE)
                if (!is.character(args$S) & !is.numeric(args$S)) 
                  return(FALSE)
		    if (is.character(args$S)){
	  		#if (!is.element("t", all.vars(parse(text = args$S))))
			#	return(FALSE)
    		    } 
                if (any(unlist(lapply(args[3:6], length)) > 1)) 
                  return(FALSE)
                if (inherits(try(parse(text = args$S), silent = TRUE), 
                  "try-error")) 
                  return(FALSE)
                if (inherits(try(D(parse(text = args$S), "t"), silent = TRUE), 
                  "try-error")) 
                  return(FALSE)

                if (is.element("env", names(args))) {
                  env <- args$env
			if (!is.null(obj$vars)) if (!is.list(eval(env, obj$vars))) return(FALSE)                  
                }
                if (is.element("n", names(args))) 
			n <- args$n
                else n <- formals(FPTL)$n
                return((length(obj$x) == n) & all(obj$x >= args$t0) & 
                  all(obj$x <= args$T) & all(obj$y >= 0) & all(obj$y <= 
                  1))
            }

    return(FALSE)
}
