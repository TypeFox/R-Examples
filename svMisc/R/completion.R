completion <- function (code, pos = nchar(code), min.length = 2,
print = FALSE, types = c("default", "scintilla"), addition = FALSE, sort = TRUE,
what = c("arguments", "functions", "packages"), description = FALSE,
max.fun = 100, skip.used.args = TRUE, sep = "\n", field.sep = "\t")
{
	finalize <- function (completions) {
		## Construct a data frame with completions
		ret <- data.frame(completion = completions,
			stringsAsFactors = FALSE)
		
		## Do we add types?
		if (isTRUE(add.types)) {
			tl <- numeric(length(completions))
			tl[grep(" = $", completions)] <- 4L
			tl[grep("::$", completions)] <- 3L
			tl[grep("<-$", completions)] <- 1L
			tl[completions %in% .reserved.words] <- 5L
			tl[!tl] <- ifelse(sapply(completions[!tl],
				function(x) existsFunction(x, where = .GlobalEnv)), 1L, 2L)
			tl <- factor(tl, levels = 1:5, labels = types)
			ret <- cbind(ret, data.frame(type = tl, stringsAsFactors = FALSE))
		}
		
		## Do we add descriptions?
		if (isTRUE(description)) {
			ret <- cbind(ret, data.frame(desc = rep("", nrow(ret)),
				context = rep("", nrow(ret)), stringsAsFactors = FALSE))
						
			## Deal with packages (completions ending with ::)
			if (length(test.pack <- grep("::$", completions))) {
				pkgDesc <- function (pkg) {
					## This is to deal with completion of :, ::, ::: in pkg base
					if (grepl(":$", pkg)) return("") else
						return(packageDescription(pkg, fields = "Description"))
				}
				ret[test.pack, "desc"] <- sapply(sub(":{2,3}$", "",
					completions[test.pack]), pkgDesc)
			}

			## Deal with argument completions (ending with " = ")
			if (length(test.arg <- grep(" = ", completions))) {
				fun <- getNamespace("utils")$.CompletionEnv[["fguess"]]
				ret[test.arg, "context"] <- fun
				ret[test.arg, "desc"] <- descArgs(fun,
					sub(" = $", "", completions[test.arg]))	
			}

			## Deal with completions with "$" (excluding things like base::$)
			if (length(test.dollar <- grep("[^:]\\$", completions))) {
				elements <- completions[test.dollar]
				object <- gsub("\\$.*$", "", completions)[1]
				items <- gsub("^.*\\$", "", completions)
				pack <- .find.multiple(object)
				ret[test.dollar, "context"] <- pack
				ret[test.dollar, "desc"] <- .descData(object, items,
					package = pack)
			}

			## Deal with completions with "@" (excluding things like base::$)
			if (length(test.slot <- grep("[^:]@", completions))) {
				elements <- completions[test.slot]
				object <- gsub("@.*$", "", completions)[1]
				slots <- gsub("^.*@", "", completions)
				pack <- .find.multiple(object)
				ret[test.slot, "context"] <- pack
				ret[test.slot, "desc"] <- .descSlots(object, slots,
					package = pack)
			}

			## Deal with completions with "["
			if (length(test.square <- grep("\\[", completions))) {
				ret[test.square, "desc"] <- .descSquare(completions[test.square],
					package = pack)
			}
		
			## TODO: do not know what to do with these?
			test.others <- grep(" ", completions)
			## TODO: are there other kind of completions I miss here?

			## Deal with function completions
			test.fun <- setdiff(1:length(completions), c(test.arg, test.pack,
				test.others, test.dollar, test.slot, test.square))
			if (length(test.fun)) {
				funs <- completions[test.fun]
				## If we have nmspace::fun, or nmspace:::fun, split it
				test.nms <- grep(".+::.+", funs)
				packs <- rep("", length(funs))
				if (length(test.nms)) {
					packs[test.nms] <- sub(":{2,3}[^:]+$", "", funs[test.nms])
					funs[test.nms] <- sub("^.+:{2,3}", "", funs[test.nms])
					packs[-test.nms] <- .find.multiple(funs[-test.nms])
				} else packs <- .find.multiple(funs)
				desc.fun <- rep("", length(packs))
				## Do not try to find description for functions in those envs
				isPack <- !packs %in% c("", ".GlobalEnv", "SciViews:TempEnv",
					"Autoloads", "tools:RGUI")
				## The following code is too slow for many function
				## (it takes 6-7sec for the 1210 base::XXXX functions)
				## So, do it only if less than max.fun
				## Note, without descriptions, it takes 0.3sec on my MacBook Pro
				if (length(isPack) < max.fun)
					desc.fun[isPack] <- descFun(funs[isPack], packs[isPack])
				ret[test.fun, "context"] <- packs
				ret[test.fun, "desc"] <- desc.fun
			}
		}
		
		## Do we sort results alphabetically?
		if (isTRUE(sort)) ret <- ret[order(completions), ]
		
		## Add metadata as attributes
		attr(ret, "token") <- token
		attr(ret, "triggerPos") <- triggerPos
		attr(ret, "fguess") <- fguess
		attr(ret, "funargs") <- funargs
		attr(ret, "isFirstArg") <- isFirstArg

		if (isTRUE(print)) {
			if (is.null(ret$desc)) {
				cat(triggerPos, paste(ret$completion, ret$type, sep = field.sep),
					sep = sep)
			} else {
				cat(triggerPos, paste(ret$completion, ret$type, ret$desc,
					ret$context, sep = field.sep), sep = sep)
			}
			if (sep != "\n") cat("\n")
			invisible(ret)
		} else ret
	}

	## Do we return the type of the entry, and if yes, in which format?
	if (is.character(types[1L])) {
		types <- switch(match.arg(types),
			default = .default.completion.types,
			scintilla = .scintilla.completion.types,
			.default.completion.types)
	}
	add.types <- as.logical(!is.na(types[1L]))

	## Default values for completion context
	token <- ""
	triggerPos <- 0L
	fguess <- ""
	funargs <- list()
	isFirstArg <- FALSE

	## Is there some code provided?
	code <- paste(as.character(code), collapse = "\n")
	if (is.null(code) || !length(code) || code == "" ||
		nchar(code, type = "chars") < min.length) {
		## Just return a list of objects in .GlobalEnv
		## TODO: look if we are inside a function and list
		## local variables (code analysis is required!)
		return(finalize(ls(envir = .GlobalEnv)))
	}

	## If code ends with a single [, then look for names in the object
	if (regexpr("[^[][[]$", code) > 0) {
		## TODO: look for object names... currently, return nothing
		return(invisible(""))
	}

	## If code ends with a double [[, then, substitute $ instead and indicate
	## to quote returned arguments (otherwise, [[ is not correctly handled)!
	if (regexpr("[[][[]$", code) > 0) {
		code <- sub("[[][[]$", "$", code)
		dblBrackets <- TRUE
	} else dblBrackets <- FALSE
	
	## Save funarg.suffix and use " = " locally
	utilsNS <- getNamespace("utils")
	ComplEnv <- utilsNS$.CompletionEnv
	opts <- ComplEnv$options
	funarg.suffix <- opts$funarg.suffix
	on.exit({
		opts$funarg.suffix <- funarg.suffix
		ComplEnv$options <- opts
	})
	opts$funarg.suffix <- " = "
	ComplEnv$options <- opts

	## Calculate completion with standard R completion tools
	utilsNS$.assignLinebuffer(code)
	utilsNS$.assignEnd(pos)
	utilsNS$.guessTokenFromLine()
	## The standard utils:::.completeToken() is replaced by our own version:
	.completeTokenExt()
	completions <- utilsNS$.retrieveCompletions()
	triggerPos <- pos - ComplEnv[["start"]]
	token <- ComplEnv[["token"]]

	## If token is empty, we complete by using objects in .GlobalEnv by default
	if (!length(completions) && token == "") {
		triggerPos <- nchar(code, type = "chars")
		## TODO: look if we are inside a function and list
		## local variables (code analysis is required!)
		return(finalize(ls(envir = .GlobalEnv)))
	}

	## For tokens like "a[m", the actual token should be "m"
    ## completions are modified accordingly
    rx <- regexpr("[[]+", ComplEnv$token)
    if (rx > 0) {
    	## Then we need to trim out whatever is before the [ in the completion
    	## and the token
    	start <- rx + attr(rx, "match.length")
    	ComplEnv$token <- substring(ComplEnv$token, start)
    	completions <- substring(completions, start)
    }
	if (!length(completions)) return(invisible(""))

	## Remove weird object names (useful when the token starts with ".")
    i <- grep("^[.]__[[:alpha:]]__", completions)
	if (length(i) > 0) completions <- completions[-i]
    if (!length(completions)) return(invisible(""))

    ## Restrict completion for which information is gathered (speed things up)
    if (!"arguments" %in% what)
		completions <- completions[regexpr("=$", completions) < 0]
    if (!length(completions)) return(invisible(""))

    if (!"packages" %in% what)
		completions <- completions[regexpr("::$", completions) < 0]
    if (!length(completions)) return(invisible(""))

    if (!"functions" %in% what)
		completions <- completions[regexpr("(::|=)$", completions) > 0]
    if (!length(completions)) return(invisible(""))

	## Eliminate function arguments that are already used
	fguess <- ComplEnv$fguess
	if (skip.used.args && length(fguess) && nchar(fguess))
		completions <- completions[!(completions %in% ComplEnv$funargs)]
	if (!length(completions)) return(invisible(""))

	## Eliminate function names like `names<-`
	i <- grep("<-.+$", completions)
	if (length(i) > 0) completions <- completions[-i]

	## Do we return only additional strings for the completion?
	if (isTRUE(addition) && triggerPos > 0L)
		completions <- substring(completions, triggerPos + 1)

	## In case of [[, restore original code
	if (dblBrackets) {  # Substitute var$name by var[["name"
		completions <- sub("[$](.+)$", '[["\\1"', completions)
		token <- sub("[$]$", "[[", token)
		triggerPos <- triggerPos + 1
	}

	## Finalize processing of the completion list
	funargs <- ComplEnv$funargs
	isFirstArg <- ComplEnv$isFirstArg
	finalize(completions)
}

.reserved.words <- c("if", "else", "repeat", "while", "function", "for", "in",
	"next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_",
	"NA_real_", "NA_complex_", "NA_character_")

.default.completion.types <- list(fun = "function", var = "variable",
	env = "environment", args = "arg", keyword = "keyword")

.scintilla.completion.types <- list(fun = "1", var = "3",
	env = "8", args = "11", keyword = "13")

## Modified utils:::inFunction()
## (checked equivalent with R 2.11.1)
## Only difference: it also gets current arguments list (if applicable).
## They are assigned to utils:::.CompletionEnv$funargs
.inFunctionExt <-
function (line, cursor)
{
	utilsNS <- getNamespace("utils")
	if (missing(line)) line <- utilsNS$.CompletionEnv[["linebuffer"]]
	if (missing(cursor)) cursor <- utilsNS$.CompletionEnv[["start"]]
		
	parens <- sapply(c("(", ")"), function(s)
		gregexpr(s, substr(line, 1L, cursor), fixed = TRUE)[[1L]],
		simplify = FALSE)
	parens <- lapply(parens, function(x) x[x > 0])
	temp <- data.frame(i = c(parens[["("]], parens[[")"]]),
		c = rep(c(1, -1), sapply(parens, length)))
	if (nrow(temp) == 0)
		return(character(0L))
	temp <- temp[order(-temp$i), , drop = FALSE]
	wp <- which(cumsum(temp$c) > 0)
	if (length(wp)) {
		index <- temp$i[wp[1L]]
		prefix <- substr(line, 1L, index - 1L)
		suffix <- substr(line, index + 1L, cursor + 1L)
		if ((length(grep("=", suffix, fixed = TRUE)) == 0L) &&
			(length(grep(",", suffix, fixed = TRUE)) == 0L))
			utilsNS$setIsFirstArg(v = TRUE)
		if ((length(grep("=", suffix, fixed = TRUE))) && (length(grep(",",
			substr(suffix, tail(gregexpr("=", suffix, fixed = TRUE)[[1L]],
			1L), 1000000L), fixed = TRUE)) == 0L)) {
			return(character(0L))
		} else {
			## This is the code added to utils:::inFunction()
			wp2 <- rev(cumsum(temp$c[-(wp[1L]:nrow(temp))]))
			suffix <- sub("^\\s+", "", suffix, perl = TRUE)
			## TODO: simplify this:
			if (length(wp2)) {
				funargs <- strsplit(suffix,	"\\s*[\\(\\)][\\s,]*",
					perl = TRUE)[[1]]
				funargs <- paste(funargs[wp2 == 0], collapse = ",")
			} else {
				funargs <- suffix
			}
			funargs <- strsplit(funargs, "\\s*,\\s*", perl=TRUE)[[1]]
			funargs <- unname(sapply(funargs, sub, pattern = "\\s*=.*$",
				replacement = utilsNS$.CompletionEnv$options$funarg.suffix,
					perl=TRUE))
			assign("funargs", funargs, utilsNS$.CompletionEnv)
			## TODO: how to take non named arguments into account too?
			## ... addition ends here

			possible <- suppressWarnings(strsplit(prefix, utilsNS$breakRE,
				perl = TRUE))[[1L]]
			possible <- possible[possible != ""]
			if (length(possible)) {
				return(tail(possible, 1))
			} else {
				return(character(0L))
			}
		}
	} else {
		return(character(0L))
	}
}

## Modified utils:::.completeToken()
## (checked equivalent with R 2.11.1)
## Main difference is that calls .inFunctionExt instead of utils:::inFunction
## and it also makes sure completion is for Complete in 'Complete("anova(", )'!
.completeTokenExt <- function () {
	utilsNS <- getNamespace("utils")
	ComplEnv <- utilsNS$.CompletionEnv
	text <- ComplEnv$token
	linebuffer <- ComplEnv$linebuffer
	st <- ComplEnv$start

	if (utilsNS$isInsideQuotes()) {
		probablyNotFilename <- (st > 2L &&
			(substr(linebuffer, st - 1L, st - 1L) %in% c("[", ":", "$")))
		if (ComplEnv$settings[["files"]]) {
			if (probablyNotFilename) {
				ComplEnv[["comps"]] <- character(0L)
			} else {
				ComplEnv[["comps"]] <- utilsNS$fileCompletions(text)
			}
			utilsNS$.setFileComp(FALSE)
		} else {
			ComplEnv[["comps"]] <- character(0L)
			utilsNS$.setFileComp(TRUE)
		}
	} else {

		## Completion does not a good job when there are quoted strings,
		## e.g for linebuffer = "Complete("anova(", )" would give arguments for
		## anova rather than for Complete.
		# Replace quoted strings with sequences of "_" of the same length.
		# This is a temporary solution though, there should be a better way...
		mt <- gregexpr('(?<!\\\\)(["\']).*?((?<!\\\\)\\1|$)', linebuffer,
			perl = TRUE)[[1]]
		if (mt[1L] != -1) {
			ml <- attr(mt, "match.length")
			y <- sapply(lapply(ml, rep, x = "a"), paste, collapse = "")
			for (i in seq_along(mt))
				substr(linebuffer, mt[i], mt[i] + ml[i]) <- y[i]
		}
		## ... additions until here

		utilsNS$.setFileComp(FALSE)
		utilsNS$setIsFirstArg(FALSE)
		guessedFunction <- ""
		if (ComplEnv$settings[["args"]]) {
			## Call of .inFunctionExt() instead of utils:::inFunction()
			guessedFunction <- .inFunctionExt(linebuffer, st)
		} else {
			guessedFunction <- ""
		}

		assign("fguess", guessedFunction, ComplEnv)
		fargComps <- utilsNS$functionArgs(guessedFunction, text)

		if (utilsNS$getIsFirstArg() && length(guessedFunction) &&
			guessedFunction %in% c("library", "require", "data")) {
			assign("comps", fargComps, ComplEnv)
			return()
		}
		lastArithOp <- tail(gregexpr("[\"'^/*+-]", text)[[1L]], 1)
		if (haveArithOp <- (lastArithOp > 0)) {
			prefix <- substr(text, 1L, lastArithOp)
			text <- substr(text, lastArithOp + 1L, 1000000L)
		}
		spl <- utilsNS$specialOpLocs(text)
		if (length(spl)) {
			comps <- utilsNS$specialCompletions(text, spl)
		} else {
			appendFunctionSuffix <- !any(guessedFunction %in%
				c("help", "args", "formals", "example", "do.call",
				"environment", "page", "apply", "sapply", "lapply",
				"tapply", "mapply", "methods", "fix", "edit"))
			comps <- utilsNS$normalCompletions(text,
				check.mode = appendFunctionSuffix)
		}
		if (haveArithOp && length(comps))
			comps <- paste(prefix, comps, sep = "")
		comps <- c(comps, fargComps)
		assign("comps", comps,  ComplEnv)
	}
}

## Similar to "find" but `what` can be a vector
## also, this one only searches in packages (position of the search path
## matching '^package:') and only gives one result per what
.find.multiple <- function (what)
{
    stopifnot(is.character(what))
    sp <- grep( "^package:", search(), value = TRUE)
    out <- rep( "" , length(what))
    for (i in sp) {
        ok <- what %in% ls(i, all.names = TRUE) & out == ""
        out[ok] <- i
        if (all(out != "")) break
    }
    names(out) <- what
    sub("^package:", "", out)
}
