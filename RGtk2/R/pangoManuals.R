# reason: we need to hide the text length parameter
pangoGetLogAttrs <-
function(text, level, language)
{
        text <- as.character(text)
        level <- as.integer(level)
        checkPtrType(language, "PangoLanguage")

        w <- .RGtkCall("S_pango_get_log_attrs", text, level, language)

        return(invisible(w))
}
pangoGlyphStringGetLogicalWidths <-
function(object, text, embedding.level)
{
        checkPtrType(object, "PangoGlyphString")
        text <- as.character(text)
        embedding.level <- as.integer(embedding.level)

        w <- .RGtkCall("S_pango_glyph_string_get_logical_widths", object, text, embedding.level)

        return(w)
}

# reason: removed text length and attrs length parameters
pangoBreak <-
function(text, analysis)
{
        text <- as.character(text)
        checkPtrType(analysis, "PangoAnalysis")

        w <- .RGtkCall("S_pango_break", text, analysis)

        return(invisible(w))
}

# reason: PangoMatrix is a boxed type that needs to be treated as a reference,
# but there is no way to create one from an API function
# This mimics the aliased static initializer in the API.

pangoMatrixInit <-
function()
{
	w <- .RGtkCall("S_pango_matrix_init")
	
	return(w)
}

# Here are some Pango constants
PANGO_SCALE <- 1024
PANGO_SCALE_XX_SMALL <-  (1 / (1.2 * 1.2 * 1.2))
PANGO_SCALE_X_SMALL <- (1 / (1.2 * 1.2))
PANGO_SCALE_SMALL <- (1 / 1.2)
PANGO_SCALE_MEDIUM <- 1
PANGO_SCALE_LARGE <- 1 * 1.2
PANGO_SCALE_X_LARGE <- 1 * 1.2 * 1.2
PANGO_SCALE_XX_LARGE <- 1 * 1.2 * 1.2 * 1.2

# reason: now let's handle the functions where the text length is omitted
pangoLayoutSetMarkupWithAccel <-
function(object, markup, accel.marker)
{
        checkPtrType(object, "PangoLayout")
        markup <- as.character(markup)
        length <- -1
        accel.marker <- as.character(accel.marker)

        w <- .RGtkCall("S_pango_layout_set_markup_with_accel", object, markup, length, accel.marker)

        return(w)
}
pangoParseMarkup <-
function(markup.text, accel.marker, .errwarn = TRUE)
{
        markup.text <- as.character(markup.text)
        length <- -1
        accel.marker <- as.character(accel.marker)

        w <- .RGtkCall("S_pango_parse_markup", markup.text, length, accel.marker)

        w <- handleError(w, .errwarn)
        
        return(w)
}
pangoGlyphStringIndexToX <-
function(object, text, analysis, index, trailing)
{
        checkPtrType(object, "PangoGlyphString")
        text <- as.character(text)
        length <- nchar(text)
        checkPtrType(analysis, "PangoAnalysis")
        index <- as.integer(index)
        trailing <- as.logical(trailing)

        w <- .RGtkCall("S_pango_glyph_string_index_to_x", object, text, length, analysis, index, trailing)

        return(w)
}
pangoGlyphStringXToIndex <-
function(object, text, analysis, x.pos)
{
        checkPtrType(object, "PangoGlyphString")
        text <- as.character(text)
        length <- nchar(text)
        checkPtrType(analysis, "PangoAnalysis")
        x.pos <- as.integer(x.pos)

        w <- .RGtkCall("S_pango_glyph_string_x_to_index", object, text, length, analysis, x.pos)

        return(invisible(w))
}
pangoShape <-
function(text, analysis, glyphs)
{
        text <- as.character(text)
        length <- nchar(text)
        checkPtrType(analysis, "PangoAnalysis")
        checkPtrType(glyphs, "PangoGlyphString")

        w <- .RGtkCall("S_pango_shape", text, length, analysis, glyphs)

        return(w)
}
# reason: var-args, just set aligns/positions for each one
pangoTabArrayNewWithPositions <-
function(size, positions.in.pixels, ...)
{
        size <- as.integer(size)
        positions.in.pixels <- as.logical(positions.in.pixels)
        args <- list(...)
		alignments <- args[seq(1,length(args),by=2)]
		positions <- args[seq(2,length(args),by=2)]
		
		array <- pangoTabArrayNew(size, positions.in.pixels)
		for (i in 1:length(args))
			pangoTabArraySetTab(array, i-1, alignments[[i]], positions[[i]])
        
		return(array)
}

# these are macros that seem simple enough to reimplement in R
# of course, if they change this will also need to be updated

pangoAscent <- function(x)
{
	rect <- as.PangoRectangle(x)
	-rect$y
}
pangoDescent <- function(x)
{
	rect <- as.PangoRectangle(x)
	rect$y + rect$height
}
pangoLbearing <- function(x)
{
	rect <- as.PangoRectangle(x)
	rect$x
}
pangoRbearing <- function(x)
{
	rect <- as.PangoRectangle(x)
	rect$x + rect$width
}

pangoReorderItems <-
function(logical.items)
{
	.notimplemented("wouldn't work automatically due to memory management details and the API says its useless")
}

## version checking ##

boundPangoVersion <- function() {
  as.numeric_version(paste(.RGtkCall("boundPangoVersion"), collapse="."))
}

checkPango <- function(version) {
  .Deprecated("boundPangoVersion() >= version")
  boundPangoVersion() >= version
}
