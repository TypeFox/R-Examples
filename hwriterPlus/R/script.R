script =
    local({
        warningCalls = vector("list", 50)
	warningMessages = character(50)
	nwarnings = 0
	renewwarnings = TRUE
	newwarnings = FALSE
        readLine =
	    function(prompt) {
	        cat(prompt)
	        flush(stdout())
	        readLines(n = 1)
	    }
	incompleteParse =
	    function(e)
	    (inherits(e, "error") &&
	     grepl("unexpected end of input", e$message))
	handleParseError =
	    function(e) {
	        msg = strsplit(conditionMessage(e), "\n")[[1]]
	        errortxt = msg[1]
	        msg = gsub("[0-9]+: ", "", msg[-c(1, length(msg))])
	        msg = msg[length(msg) - 1:0]
	        if (length(msg) == 1)
	            msg = paste(" in: \"", msg, "\"\n", sep = "")
	        else
	            msg = paste(" in:\n\"",
	                paste(msg, collapse = "\n"),
	                "\"\n", sep = "")
	        cat("Error",
	            gsub("\n.*", "",
	                 gsub("<text>:[0-9]+:[0-9]+", "",
	                      errortxt)),
	            msg, sep = "")
	    }
	handleError =
	    function(e) {
	        cat("Error in", deparse(conditionCall(e)),
	            ":", conditionMessage(e), "\n")
	    }
	handleValue =
	    function(e) {
	        if (e$visible) {
	            print(e$value)
	        }
	    }
	echoCommands =
	    function(cmd, transcon) {
	        cat(paste(c("> ",
	                    rep("+ ", max(length(cmd) - 1), 0)),
	                  cmd, "\n", sep = ""), sep = "",
	            file = transcon)
	    }
	echoOutput =
	    function(transcon, outcon) {
	        seek(outcon, 0)
	        lines = readLines(outcon, warn = FALSE)
	        writeLines(lines, transcon)
	        seek(outcon, 0)
	        truncate(outcon)
	    }
	warningHandler = function(w) {
	    newwarnings <<- TRUE
	    if (renewwarnings) {
	        renewwarnings <<- FALSE
	        nwarnings <<- 0
	    }
	    n = nwarnings + 1
	    if (n <= 50) {
	        warningCalls[[n]] <<- conditionCall(w)
	        warningMessages[n] <<- conditionMessage(w)
	        nwarnings <<- n
	    }
	    invokeRestart("muffleWarning")
	}
	tryCatchWithWarnings =
	    function(expr)
	    withCallingHandlers(tryCatch(expr,
	            error = function(e) e),
	        warning = warningHandler)
	displayWarnings =
	    function(n) {
	        if (n <= 10)
	            print(warnings())
	        else if (n < 50) {
	            cat("There were",
	                nwarnings,
	                "warnings (use warnings() to see them)\n")
	        }
	        else
	            cat("There were 50 or more warnings",
	                "(use warnings() to see the first 50)\n")
	    }
	isQuitCall =
	    function(e)
	    (!inherits(e, "error") &&
	     length(e) == 1 &&
	     deparse(e[[1]], nlines = 1) == "q()")
        repl =
	    function(env, transcon, outcon) {
	        sinkdepth = sink.number()
	        prompt = "script> "
	        cmd = character()
	        repeat {
	            ans = tryCatch(repeat {
		        repeat {
			    cmd = c(cmd, readLine(prompt))
			    ans = tryCatch(parse(text = cmd),
			        error = function(e) e)
			    if (inherits(ans, "error")) {
			        if (incompleteParse((ans))) {
			            prompt = "script+ "
			        }
			        else {
			            echoCommands(cmd, transcon)
			            sink(outcon, split = TRUE)
			            handleParseError(ans)
			            sink()
			            echoOutput(transcon, outcon)
			            prompt = "script> "
			            cmd = character()
			        }
			    }
			    else {
			        echoCommands(cmd, transcon)
			        if (length(ans) == 0) {
				    break
				}
				else if (isQuitCall(ans)) {
				    return()
				}
				else if (grepl("^script\\(",
				               deparse(ans[[1]], nlines = 1))) {
				    sink(outcon, split = TRUE)
				    cat("Error: You can't call \"script\" while scripting\n")
				    sink()
				    echoOutput(transcon, outcon)
				    break
				}
			        else {
				    renewwarnings <<- TRUE
				    newwarnings <<- FALSE
				    for(e in ans) {
				        sink(outcon, split = TRUE)
				        e = tryCatchWithWarnings(withVisible(eval(e,
				            envir = env)))
				        if (inherits(e, "error"))
				            handleError(e)
				        else
				            handleValue(e)
				        sink()
				        echoOutput(transcon, outcon)
				    }
				    if (newwarnings) {
				        warnings = warningCalls
				        names(warnings) = warningMessages
				        assign("last.warning",
				               warnings[1:nwarnings],
				               "package:base")
				        sink(outcon, split = TRUE)
				        displayWarnings(nwarnings)
				        sink()
				        echoOutput(transcon, outcon)
				    }
				}
			        prompt = "script> "
			        cmd = character()
			    }
			}
		    }, interrupt = function(x) x)
		    if (inherits(ans, "interrupt")) {
		        if (sink.number() > sinkdepth) {
		            sink()
		            echoOutput(transcon, outcon)
		        }
		        else
		            echoCommands(cmd, transcon)
		        cat("\nInterrupt!\n")
		        cat("Interrupt!\n", file = transcon)
		        prompt = "script> "
		        cmd = character()
		    }
		    else
		        stop("Interrupt catcher caught non-interrupt")
	        }
	    }

        function(file = "transcript.txt") {
	    ## if (!isatty(stdin()))
	    ##     stop("script can only be used interactively")
	    transcon = file(file, "w")
	    outcon = file("")
	    cat("Script started, file is \"", file, "\"\n", sep = "")
	    cat("Script started on", date(), "\n", file = transcon)
	    repl(sys.parent(), transcon, outcon)
	    cat("Script done on", date(), "\n", file = transcon)
	    cat("Script done, file is \"", file, "\"\n", sep = "")
	    close(outcon)
	    close(transcon)
	    invisible()
	}
    })
