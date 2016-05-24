### tk2Fonts.R - Manage Tk fonts
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2009-04-23: fix of a bug in fntl$family, $ operator not allowed (thanks Brian Ripley)
### - 2007-01-11: fisrt version (for tcltk2_1.0-0)
###
### To do:
### -

tk2font.get <- function (font, what = c("family", "size", "bold", "italic"))
{
	## font is the TkFont name to use, in case of several items, other ones
	## are secondary, tertiary, ... options
	## what indicate what characteristic of the font to report in the list
	## 'family', 'size', 'bold', 'italic', 'underline', 'overstrike' (last two rarely used)
	if (!is.tk()) return("")
	allTkFonts <- as.character(tkfont.names())
	for (fnt in font) {
		if (fnt %in% allTkFonts) break
	}
	if (!fnt %in% allTkFonts) {
		return("")
		##if (length(font) == 1) {
		##	stop("'", font, "' is not currently defined in Tk")
		##} else {
		##	stop("'", paste(font, collapse = "', '"), "' are not currently defined in Tk")
		##}
	}
	fontspec <- as.character(tkfont.configure(fnt))
	res <- list()
	if (length(fontspec) != 12) return(res)	# There is a problem here!
	if ("family" %in% what) res$family <- fontspec[2]
	if ("size" %in% what) res$size <- as.numeric(fontspec[4])
	if ("bold" %in% what) res$bold <- (fontspec[6] == "bold")
	if ("italic" %in% what) res$italic <- (fontspec[8] == "italic")
	if ("underline" %in% what) res$underline <- (fontspec[10] == "1")
	if ("overstrike" %in% what) res$overstrike <- (fontspec[12] == "1")
	return(res)
}

tk2font.set  <- function (font, settings)
{
### TODO: allow for multiple fonts specifications => take first one available
	## font is the name of the TkFont to create/change
	## settings is a list with font characteristics
	if (!is.tk()) return(NULL)
	font <- as.character(font)
	l <- length(font)
	if (!is.list(settings) && !is.character(settings))
		stop("'settings' must be a list or a character string")
	## If settings is a character string,
	## it is assumed to be a text description of a Tk font
	if (is.character(settings)) {
		## Do not recycle... make sure that lengths match
		if (length(settings) != l)
			stop("length of 'font' and 'settings' do not match")
		for (i in 1:l) {
			.Tcl(paste("catch {font create ", font[i], "}", sep = ""))
			.Tcl(paste("catch {font configure ", font[i], " ", settings[i], "}",
				sep = ""))
		}
	} else {  # This is a list of font characteristics
		## Do not recycle... make sure that lengths match
		if (l > 1) {
			if (length(settings) != l)
				stop("length of 'font' and 'settings' do not match")
		} else {  # Is it the list of characteristics, or a lit containing it?
			if (any(names(settings) %in%
				c("family", "size", "bold", "italic", "underline", "overstrike")))
				settings <- list(settings)
		}
		fntfamilies <- as.character(tkfont.families())
		for (i in 1:l) {
			## Construct the font descriptor
			fntl <- as.list(settings[[i]])
			fnt <- " "
			if (!is.null(fntl$family)) {
				## Look for the first font family provided that is available
				fntfamily <- fntl$family
				if (length(fntfamily) > 1) {
					fntexists <- fntfamily %in% fntfamilies
					if (any(fntexists)) fntfamily <- fntfamily[fntexists][1] else
						fntfamily <- fntfamily[1]  # No fonts found... take first one
				}
			 	fnt <- paste(fnt, "-family {", fntfamily, "}", sep = "")
			}
			if (!is.null(fntl$size)) fnt <- paste(fnt, "-size", fntl$size)
			if (!is.null(fntl$bold)) {
				value <- if(fntl$bold == TRUE) "bold" else "normal"
				fnt <- paste(fnt, "-weight", value)
			}
			if (!is.null(fntl$italic)) {
				value <- if(fntl$italic == TRUE) "italic" else "roman"
				fnt <- paste(fnt, "-slant", value)
			}
			if (!is.null(fntl$underline)) fnt <- paste(fnt, "-underline",
				as.numeric(fntl$underline == TRUE))
			if (!is.null(fntl$overstrike)) fnt <- paste(fnt, "-overstrike",
				as.numeric(fntl$overstrike == TRUE))
			## Possibly create the font in Tk
			.Tcl(paste("catch {font create ", font[i], "}", sep = ""))
			if (fnt != " ")
				.Tcl(paste("catch {font configure ", font[i], fnt, "}", sep = ""))
		}
	}
	res <- font %in% as.character(tkfont.names())
	names(res) <- font
	return(res)
}

tk2font.setstyle <- function (text = TRUE, system = FALSE, default.styles = FALSE)
{
	## Set default fonts according to currently defined style
	## .SystemFonts and .Fonts must be defined in SciViews:TempEnv!

	if (!is.tk()) {
		warning("Package Tk is required but not loaded")
		return(NULL)
	}

	## This is a copy of assignTemp(), getTemp() and existsTemp() functions from
	## svMisc, so that we do not link to this package
	TempEnv <- function () {
	    pos <-  match("SciViews:TempEnv", search())
	    if (is.na(pos)) { # Must create it
	        `SciViews:TempEnv` <- list()
	        Attach <- function (...) get("attach", mode = "function")(...)
			Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
	        rm(`SciViews:TempEnv`)
	        pos <- match("SciViews:TempEnv", search())
	    }
	    return(pos.to.env(pos))
	}
	
	assignTemp <- function (x, value, replace.existing = TRUE)
	    if (replace.existing || !exists(x, envir = TempEnv(), mode = "any", inherits = FALSE))
	        assign(x, value, envir = TempEnv())
	
	existsTemp <- function (x, mode = "any")
	    exists(x, envir = TempEnv(), mode = mode, inherits = FALSE)
	
	getTemp <- function (x, default = NULL, mode="any") {
	    if  (exists(x, envir = TempEnv(), mode = mode, inherits = FALSE)) {
	        return(get(x, envir = TempEnv(), mode = mode, inherits = FALSE))
	    } else {  # Variable not found, return the default value
	        return(default)
	    }
	}

	if (system) {  # Set system fonts
		## We collect back system fonts settings (other values may be imposed by Tk)
		sysfonts <- list(
			defaultclassic = tk2font.get("TkClassicDefaultFont"),
			default = tk2font.get("TkDefaultFont"),
			caption = tk2font.get("TkCaptionFont"),
			smallcaption = tk2font.get(c("TkSmallCaptionFont", "TkCaptionFont")),
			menu = tk2font.get(c("TkMenuFont", "TkDefaultFont")),
			status = tk2font.get(c("TkStatusFont", "TkTooltipFont")),
			tooltip = tk2font.get("TkTooltipFont"),
			heading = tk2font.get("TkHeadingFont"),
			icon = tk2font.get(c("TkIconFont", "TkDefaultFont"))
		)
		## Make sure these are correctly defined
		assignTemp(".SystemFonts", sysfonts)
		res <- TRUE
	} else res <- character(0)

	if (default.styles) {  # Define default styles
		## These are the four default Font themes one can use
		assignTemp(".FontsStyleClassic", list(
			Text = list(family = c("Times New Roman", "Times"), size = -12),
			Title = list(family = c("Arial", "Helvetica"), size = -14, bold = TRUE),
			BigTitle = list(family = c("Arial", "Helvetica"), size = -16, bold = TRUE),
			Fixed = list(family = c("Courier New", "Courier"), size = -12)
		))

		assignTemp(".FontsStyleAlternate", list(
			Text = list(family = "Georgia", alt.family = "Times", size = -12),
			Title = list(family = c("Trebuchet MS", "Trebuchet"),
				alt.family = "Helvetica", size = -14, bold = TRUE),
			BigTitle = list(family = c("Trebuchet MS", "Trebuchet"),
				alt.family = "Helvetica", size = -16, bold = TRUE),
			Fixed = list(family = "Andale Mono", alt.family = "Courier", size = -12)
		))

		assignTemp(".FontsStylePresentation", list(
			Text = list(family = "Verdana", alt.family = "Helvetica", size = -12),
			Title = list(family = "Verdana", alt.family = "Helvetica", size = -14,
				bold = TRUE),
			BigTitle = list(family = "Verdana", alt.family = "Helvetica", size = -16,
				bold = TRUE),
			Fixed = list(family = "Lucida Console", alt.family = "Courier",
				size = -12)
		))

		assignTemp(".FontsStyleFancy", list(
			Text = list(family = c("Trebuchet MS", "Trebuchet"),
				alt.family = "Helvetica", size = -12),
			Title = list(family = c("Comic Sans MS", "Comic Sans"),
				alt.family = "Helvetica", size = -14, bold = TRUE),
			BigTitle = list(family = c("Comic Sans MS", "Comic Sans"),
				alt.family = "Helvetica", size = -16, bold = TRUE),
			Fixed = list(family = "Lucida Console", alt.family = "Courier",
				size = -12)
		))
	}

	if (text) {  # Set text, titles and fixed fonts
		## Determine which font style we currently use
		curStyle <- getTemp(".FontsStyle", default = "Classic", mode = "character")
		curSFonts <- getTemp(paste(".FontsStyle", curStyle, sep = ""),
			default = getTemp(".FontsStyleClassic"))
		assignTemp(".Fonts", curSFonts)

		## Create corresponding fonts in Tk (note, we create bold, italic, and
		## bolditalic equivalents for TkTextFont and TkFixedFont
		Fonts <- list()
		Fonts$Text <- curSFonts$Text
		Fonts$Text$bold <- FALSE
		Fonts$Text$italic <- FALSE
		Fonts$TextBold <- Fonts$Text
		Fonts$TextBold$bold <- TRUE
		Fonts$TextItalic <- Fonts$Text
		Fonts$TextItalic$italic <- TRUE
		Fonts$TextBoldItalic <- Fonts$TextBold
		Fonts$TextBoldItalic$italic <- TRUE
		Fonts$Title <- curSFonts$Title
		Fonts$BigTitle <- curSFonts$BigTitle
		Fonts$Fixed <- curSFonts$Fixed
		Fonts$Fixed$bold <- FALSE
		Fonts$Fixed$italic <- FALSE
		Fonts$FixedBold <- Fonts$Fixed
		Fonts$FixedBold$bold <- TRUE
		Fonts$FixedItalic <- Fonts$Fixed
		Fonts$FixedItalic$italic <- TRUE
		Fonts$FixedBoldItalic <- Fonts$FixedBold
		Fonts$FixedBoldItalic$italic <- TRUE

		FNames <- c("TkTextFont", "TkTextBoldFont", "TkTextItalicFont",
			"TkTextBoldItalicFont", "TkTitleFont", "TkBigTitleFont", "TkFixedFont",
			"TkFixedBoldFont", "TkFixedItalicFont", "TkFixedBoldItalicFont")
		res <- c(res, tk2font.set(FNames, Fonts))
	}
	## Check the results
	if (system && any(!res))
		warning("One or several Tk fonts not set: '",
			paste(names(res)[!res], collapse = "', '", "'"))
	return(res)
}
