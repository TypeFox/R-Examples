#'  Remove Redundant Spaces
#'
#' 	This function removes redundant spaces from character vectors
#'
#' 	@param x character vector
#'
#'  @return  character vector without trailing or multiple spaces
#'  @examples stopifnot(confidence:::strip_spaces(" Hello  World  ") == "Hello World")
strip_spaces <- 
function(x) {
	x <- gsub(pattern = " {2,}",   replacement = " ", x = x)
	x <- gsub(pattern = "^ +| +$", replacement = "",  x = x)
	x
}


#' Sanitize Text to Give Proper Filenames
#'
#' @param x character vector to sanitize
#' 
#' @return sanitized character vector
sanitize <- 
function(x) {
    gsub(
        pattern = "[^a-zA-Z0-9_-]", 
        replacement = "_", 
        x = x,
        perl = TRUE
    ) 
}


#'  Transformations
#'  
#'  Performs log or logit transformations.
#'
#'  @param x value to transform
#'  @param type type of transform (log, logit).
#'
#'  @return transformed value
transform <- 
function(x, type = c("identity", "log", "logit", "none", NA_character_)) {
    switch(
        match.arg(type),
        log = log(x),
        logit = log(x/(1-x)),
        x
    )
}


#'  Back-transformations
#'  
#'  Performs inverse log or logit transformations.
#'
#'  @param x value to back-transform
#'  @param type type of transform (log, logit).
#'
#'  @return backtransformed value
backtransform <- 
function(x, type = c("identity", "log", "logit", "none", NA_character_)) {
    switch(
        match.arg(type),
        log = exp(x),
        logit = exp(x) / (1 + exp(x)),
        x
    )
}


#'  Check Confidence data
#'
#'     This function checks \code{data.frames} to be used by the confidence 
#'     package. The format has been specified in Van Loon (2014) and should
#'     contain the following columns:
#'    \itemize{
#' 		    \item{OBJECTID: water body code, e.g., NL89_os;}
#'  	    \item{PAR: parameter, e.g., Cadmium;}
#'  	    \item{DATE: date according to ISO 8601 (YYYY-mm-dd) for point values
#'              or year YYYY for annual means;}
#'          \item{VALUE: numerical value.}
#'          \item{TARGET: target value for the European Water Framework 
#'              Directive;}
#'          \item{UNIT: measurement unit of PAR. This unit should be the same for all
#'              records with the same PAR and is the same for both VALUE 
#'              and TARGET;}
#'          \item{transform: data transformation, i.e., log, logit, NA.}
#'      }
#'
#'  @details The function performs the following tasks:
#' 	    \itemize{
#'  	    \item{checks availablitity of required columns (case insensitive);}
#'          \item{make column names case-insensitive;}
#'  	    \item{removes redundant spaces;}
#'          \item{checks on missing values in required columns;}
#'          \item{checks if DATE-field adheres to ISO 8601 (YYYY-mm-dd) or YYYY;}
#'          \item{checks mixtures of annual averages and point values for a each year;}
#'          \item{checks if measurement units are the same for a specific 
#'              OBJECTID-PAR-pair;}
#'          \item{checks if TARGET-value is the same for a specific 
#'              OBJECTID-PAR-pair;}
#'          \item{checks if transform is one of log, logit, NA in transform column;}
#'          \item{checks that the EQR-column contains identical values fo each
#'              OBJECTID-PAR combination.}
#'          
#'      }
#'      
#'  @param x \code{data.frame} to be checked
#'      
#'  @return \code{data.frame} that has passed all checks
#'  
#'  @export
#'  
conf_input <- function(x) {

    # check column names (case insensitive)
	requiredColumns <- c("OBJECTID", "PAR", "YEAR", "VALUE", "TARGET",  "UNIT")
    names(x) <- toupper(names(x))
	missingColumns <- setdiff(requiredColumns, names(x))
	if (length(missingColumns) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing\n: %s", 
                toString(missingColumns)
            ),
            call. = FALSE
        )
	}

    # check on missing values in required columns (except UNIT)
    isNA <- apply(
        X = is.na(x[, requiredColumns[requiredColumns != "UNIT"]]), 
        MARGIN = 1L, 
        FUN = any
    )
    if (any(isNA)) {
        message <- sprintf(
            paste0(
                "missing values found in record(s): %s\n",
                "These record(s) will be removed\n"
            ),
            toString(which(isNA))
        )
        warning(message, call. = FALSE)
        x <- x[!isNA, ]
    }

	# remove redundant spaces
	x <- as.data.frame(
		lapply(X = x, FUN = function(x) {
			if (is.character(x)) {
				x <- strip_spaces(x)
			}
			x
		}),
        stringsAsFactors = FALSE
	)

    # check UNIT field
    tmp <- unique(x[, c("OBJECTID", "PAR", "UNIT")])
    isDuplicated <- duplicated(tmp[, c("OBJECTID", "PAR")])
    if (any(isDuplicated)) {
        stop(
            "All values in column 'UNIT' should be identical",
            " for 'OBJECTID'-'PAR' pairs", 
            call. = FALSE)
    }

    # check TARGET field
    tmp <- unique(x[, c("OBJECTID", "PAR", "TARGET")])
    isDuplicated <- duplicated(tmp[, c("OBJECTID", "PAR")])
    if (any(isDuplicated)) {
        stop(
            "All values in column 'TARGET' should be identical",
            " for 'OBJECTID'-'PAR' pairs", 
            call. = FALSE)
    }

    # check column 'YEAR'
    x$YEAR <- suppressWarnings(as.integer(x$YEAR))
    isValidYear <- (x$YEAR > 1000) & (x$YEAR < 9999)
    if (!all(isValidYear)) {
    	stop(
			"Invalid years found. Years should consist of 4 digits (YYYY)",
			call. = FALSE
		)
    }

    # check that only one annual average is given for each year
    if (anyDuplicated(x[, c("OBJECTID", "PAR", "YEAR")])) {
        stop(
            "Only one annual average allowed in each year", 
            call. = FALSE
        )
    }

    # check data transformation
    if (is.null(x$TRANSFORM)) {
        x$TRANSFORM <- "none"
    } else {
        validTransforms <- c("log", "logit", "identity", "none", NA_character_)
        if (!(all(x$TRANSFORM %in% validTransforms))) {
        	stop(
    			"Invalid transformations found in column 'transform'.",
                " Only 'log', 'logit' or NA are allowed.",
    			call. = FALSE
    		)
        }
        x$TRANSFORM[is.na(x$TRANSFORM)] <- "identity"
    }
    
    # check type of COLOR column
    if (is.null(x$COLOR)) {
        x$COLOR <- "green/orange"
    } else {
        x$COLOR <- tolower(x$COLOR)
        if (!all(x$COLOR %in% c("green/orange", "orange/green"))) {
            stop(
    			"'COLOR' column should be either",
                " 'green/orange' or 'orange/green'. ",
    			call. = FALSE
    		)
        }
    }
    
    # check that COLOR-column contains identical values fo each OBJECTID-PAR
    if (
        nrow(unique(x[, c("OBJECTID", "PAR", "COLOR")])) >
        nrow(unique(x[, c("OBJECTID", "PAR")]))) {
        stop(
    		"values in COLOR column should be identical",
            " for each OBJECTID-PAR combination.",
			call. = FALSE
		)
    }
   
    # return result
    class(x) <- c("conf_input", "data.frame")
    x
}




#' Multi-Year Average
#' 
#' Estimates the multi-year average of environmental properties and
#' associated confidence intervals.
#' 
#' @param x object of class \code{\link{conf_input}} or a \code{\link{data.frame}}
#'  that can be coerced to an instance of class \code{\link{conf_input}}.
#' @param \dots further arguments to be passed to other methods
#' 
#' @return a \code{\link{data.frame}} with the following columns:
#' \itemize{
#'     	\item{\code{MYA}: }{the multi-year arithmetic average;}
#'  	\item{\code{PROB_LTT}: }{the probability that \code{MYA} is less than the target value specified;}
#'  	\item{\code{PROB_GTT}: }{the probability that \code{MYA} is greater than the target value specified;}
#'      \item{\code{q05}: }{the lowerbound of the 90\% confidence interval of \code{MYA}}
#'      \item{\code{q95}: }{the upperbound of the 90\% confidence interval of \code{MYA}}
#'  }
#' 
#' @seealso \code{\link{conf}}
#'
#'  @export
mya <- function(x, ...) {
    UseMethod("mya")
}

#' @export
mya.character <- function(x, ...) {
    mya(conf_input(x))
}

#' @export
mya.data.frame <- function(x, ...) {
    mya(conf_input(x))
}

#' @export
mya.conf_input <- function(x, ...) {
        
    # optionally transform data
    x$VALUE_t  <- mapply(FUN = transform, x = x$VALUE,  type = x$TRANSFORM)
    x$TARGET_t <- mapply(FUN = transform, x = x$TARGET, type = x$TRANSFORM)

    # estimate multi-year averages
    x <- ddply(
        .data = x, 
        .variables = c("OBJECTID", "PAR"), 
        .fun = function(x) {
            N <- nrow(x) # Note: NAs have already been removed by check_input
            x <- data.frame(
                COLOR = x$COLOR[1],
                PERIOD = paste(min(x$YEAR), max(x$YEAR), sep = "-"),
                N = N,
                df = N-1L,
                MYA_t = mean(x$VALUE_t),
                STERR_t = sd(x$VALUE_t) / sqrt(N),
                TARGET_t = x$TARGET_t[1],
                UNIT = x$UNIT[1],
                TRANSFORM = x$TRANSFORM[1],
                stringsAsFactors = FALSE
            )
            x$q05_t <- x$MYA_t + qt(p = 0.05, df = x$df) * x$STERR_t
            x$q95_t <- x$MYA_t + qt(p = 0.95, df = x$df) * x$STERR_t
            x$PROB <- pt(
                q = (x$MYA_t - x$TARGET_t) / x$STERR_t,
                df = x$df, 
                lower.tail = TRUE
            )
            x
        }
    )

    # backtransform data
    x$MYA    <- mapply(FUN = backtransform, x = x$MYA_t,    type = x$TRANSFORM)
    x$TARGET <- mapply(FUN = backtransform, x = x$TARGET_t, type = x$TRANSFORM)
    x$q05    <- mapply(FUN = backtransform, x = x$q05_t,    type = x$TRANSFORM)
    x$q95    <- mapply(FUN = backtransform, x = x$q95_t,    type = x$TRANSFORM)

    # add Prob(VALUE > TARGET) and Prob(VALUE <= TARGET)
    x$PROB_GTT <- x$PROB
    x$PROB_LTT <- 1.0 - x$PROB
    x$PROB <- NULL
    
    # return result
    class(x) <- c("mya", "data.frame")
    x
}

#' @export
as.data.frame.mya <- function(x, ...) {
    class(x) <- "data.frame"
    x[c("OBJECTID", "PAR", "PERIOD", "MYA", "TARGET", 
         "PROB_LTT", "PROB_GTT", "q05", "q95")]
}



#' @export
print.mya <- function(x, ...) {
    print(as.data.frame(x))
}


#' @export
plot.mya <- function(x, which, ...) {
    
    # check 'which' argument
    if (missing(which)) {
        if (nrow(x) == 1L) {
            which <- 1L
        } else {
            stop("argument 'which' is missing", call. = FALSE)
        }
    }
    if (length(which) != 1L) {
        stop("'which' should be of length 1", call. = FALSE)
    }
    if (!(which %in% 1:nrow(x))) {
        stop(
            "'which' should be an integer >=1 and <=", nrow(x),
            call. = TRUE
        )
    }
    x <- x[which, ]

    # standardize target value (transformed scale)
    target0_t <- with(x, (TARGET_t - MYA_t) / STERR_t)

    # t-values (t(0,1, N-1)) including TARGET
    t01 <- seq(
        from = min(-5, 1.1 * target0_t), 
        to   = max( 5, 1.1 * target0_t),
        length.out = 1001
    )

    # t-values left and right of target (Note: target should be included)
    t01_l <- t01[t01 <= target0_t]
    t01_r <- t01[t01 >= target0_t]

    # number of values right and left of target
    n_l <- length(t01_l)
    n_r <- length(t01_r)

    # convert t01 to MYA on the orignal scale
    MYA_l <- backtransform(t01_l * x$STERR_t + x$MYA_t, type = x$TRANSFORM)
    MYA_r <- backtransform(t01_r * x$STERR_t + x$MYA_t, type = x$TRANSFORM)
    
    # construct data.frame for plotting
    p <- rbind(
        data.frame(
            id = if (x$COLOR == "orange/green") {"right"} else {"left"},
            OBJECTID = x$OBJECTID,
            PAR = x$PAR, 
            MYA = c(MYA_l[1], MYA_l, MYA_l[n_l]), 
            density = c(0, dt(x = t01_l, df = x$df), 0),
            TARGET = x$TARGET,
            PROB = x$PROB_GTT,
            PERIOD = x$PERIOD
        ),
        data.frame(
            id = if (x$COLOR == "orange/green") {"left"} else {"right"},
            OBJECTID = x$OBJECTID,
            PAR = x$PAR, 
            MYA = c(MYA_r[1], MYA_r, MYA_r[n_r]), 
            density = c(0, dt(x = t01_r, df = x$df), 0),
            TARGET = x$TARGET,
            PROB = x$PROB_GTT,
            PERIOD = x$PERIOD
        )
    )
    
    # create plot
    ggplot(data = p) +
        geom_polygon(
            mapping = aes_string(
                x = "MYA", 
                y = "density",
                group = "id", 
                fill = "id"
            )
        ) +
        geom_vline(
            mapping = aes_string(xintercept = "TARGET"),
            colour = "red"
        ) +
        scale_fill_manual(
            name = "area",
            # http://jfly.iam.u-tokyo.ac.jp/color/#pallet
            values = c(
                left  = rgb(  0, 158, 115, maxColorValue = 255), # green
                right = rgb(230, 159,  0, maxColorValue = 255)   # orange
            )
        ) +
        geom_text(
            mapping = aes_string(
                x = "TARGET", 
                y = "0",
                label = 'paste0(round(100*PROB), "%")'
            ),
            hjust = -0.2,
            vjust = -0.1
        ) +
        geom_text(
            mapping = aes_string(
                x = "TARGET", 
                y = "0",
                label = 'paste0(round(100-100*PROB), "%")'
            ), 
            hjust = 1.2,
            vjust = -0.1
        ) +
        scale_x_continuous(name = "multiyear average") +
        scale_y_continuous(name = "confidence density") +
        theme(legend.position = "none") +
        ggtitle(label = paste(x$OBJECTID, x$PAR, x$PERIOD, paste("target", x$TARGET), sep = "; "))
}




write_html <- function(x, ...) {
    UseMethod("write_html")
}

write_html.character <- function(x, outputDir, browse = TRUE, ...) {
    write_html(x = conf_input(x), outputDir = outputDir, browse = browse)
}

write_html.data.frame <- function(x, outputDir, browse = TRUE, ...) {
    write_html(x = conf_input(x), outputDir = outputDir, browse = browse)
}

write_html.conf_input <- function(x, outputDir, browse = TRUE, ...) {

    # create temporary directory
    tmpdir <- tempfile(pattern = "confidence")
	dir.create(path = tmpdir)

	# copy template of report to temporary directory
	templates <- list.files(
		path = system.file("Rmd", package = "confidence"),
		pattern = "\\.Rmd$", full.names = TRUE)
	file.copy(from = templates, to = tmpdir)
    
    # create Markdown document 
    # (code below works better than knit2html, see BEQI2)
    owd <- setwd(tmpdir)
    on.exit(setwd(owd), add = TRUE)
	suppressMessages(
        res <- try(knit(input = "confidence.Rmd", quiet = TRUE), silent = TRUE)
    )
    if (inherits(res, "try-error")) {
        stop(toString(attr(res, "condition")$message), call. = FALSE)
    }
	markdownToHTML(
		file  = "confidence.md", 
		output = "output.html",
        options = c("use_xhtml", "smartypants", "mathjax", "highlight_code"),
        extensions = getOption("markdown.extensions"),
    	title = "Confidence Report",
        stylesheet = system.file("css", "confidence.css", 
                                 package = "confidence")
	)

    # copy all files from the temporary directory to the output directory
    if (!file.exists(outputDir)) {
        dir.create(path = outputDir, recursive = TRUE)
    }
    tmpfiles <- list.files(path = tmpdir, pattern = "\\.html$|\\.png$", 
                           recursive = TRUE)
    file.copy(from = tmpfiles, to = outputDir)
    
	# view result
	if (isTRUE(browse)) {
		browseURL(file.path(outputDir, "output.html"))
	}
}