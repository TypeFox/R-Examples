subst <- function(strings, ..., digits=7) {
  substitutions <- list(...)
  names <- names(substitutions)
  if (is.null(names)) names <- rep("", length(substitutions))
  for (i in seq_along(names)) {
    if ((n <- names[i]) == "")
      n <- as.character(sys.call()[[i+2]])
    value <- substitutions[[i]]
    if (is.numeric(value))
      value <- formatC(value, digits=digits, width=1)
    strings <- gsub(paste("%", n, "%", sep=""), value, strings)
  }
  strings
}

.writeWebGL <- function(dir="webGL", filename=file.path(dir, "index.html"),
                       template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
                       prefix = "",
                       snapshot = TRUE, commonParts = TRUE, reuse = NULL,
		       font="Arial",
                       width = NULL, height = NULL) {

  # Lots of utility functions and constants defined first; execution starts way down there...

  header <- function()
  	if (commonParts)
  	  c(
        as.character(includeScript(system.file("htmlwidgets/lib/CanvasMatrix/CanvasMatrix.js", package = "rglwidget"))),
        as.character(includeScript(system.file("htmlwidgets/lib/rglClass/rglClass.src.js", package = "rglwidget")))
      )

  scriptheader <- function() subst(
  '
<div id="%elementId%" class="rglWebGL"></div>
<script type="text/javascript">
	var %prefix%div = document.getElementById("%elementId%"),
      %prefix%rgl = new rglwidgetClass();
  %prefix%div.width = %width%;
  %prefix%div.height = %height%;
  %prefix%rgl.initialize(%prefix%div,
                         %json%);
  %prefix%rgl.prefix = "%prefix%";
</script>', prefix, elementId, json, width, height)

  footer <- function() subst('
	<p id="%prefix%debug">
	You must enable Javascript to view this page properly.</p>',
    prefix)

  #  Execution starts here!

  # Do a few checks first

  elementId <- paste0(prefix, "div")

  if (!file.exists(dir))
    dir.create(dir)
  if (!file.info(dir)$isdir)
    stop(gettextf("'%s' is not a directory", dir), domain = NA)

  if (!is.null(template)) {
    templatelines <- readLines(template)
    templatelines <- subst(templatelines, rglVersion = packageVersion("rgl"), prefix = prefix)

    target <- paste("%", prefix, "WebGL%", sep="")
    replace <- grep( target, templatelines, fixed=TRUE)
    if (length(replace) != 1)
      stop(gettextf("template '%s' does not contain '%s'", template, target),
           domain = NA)

    result <- c(templatelines[seq_len(replace-1)], header())
  } else
    result <- header()

  scene <- convertScene(width = width, height = height,
                        elementId = elementId, reuse = reuse,
                        snapshot = snapshot)
  if (is.null(width)) width <- scene$width
  if (is.null(height)) height <- scene$height

  reuse <- attr(scene, "reuse")
  json <- toJSON(I(scene),
                 dataframe = "columns", null = "null", na = "null",
                 auto_unbox = TRUE, digits = getOption("shiny.json.digits",
                                                       7),
                 use_signif = TRUE, force = TRUE, POSIXt = "ISO8601",
                 UTC = TRUE, rownames = FALSE, keep_vec_names = TRUE)

  result <- c(result,
              scriptheader(),
              footer(),
              if (!is.null(template))
              	templatelines[replace + seq_len(length(templatelines)-replace)]
              else
              	subst("<script>%prefix%rgl.start();</script>", prefix = prefix)
             )

  cat(result, file=filename, sep="\n")
#   if (!is.null(reuse)) {
#     prefixes <- prefixes[!duplicated(prefixes$id),]
#     attr(filename, "reuse") <- prefixes
#   }
  invisible(structure(filename, reuse = reuse))
}
