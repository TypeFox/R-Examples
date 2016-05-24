##' asciiOptions
##'
##' @param select select
##' @param .backends .backends
##' @param .outputs .outputs
##' @param .extensions .extensions
##' @param .cli .cli
##' @param .args .args
##' @param .O .O
##' @param .f .f
##' @param .e .e
##' @param .preambule .preambule
##' @return options
##' @author David.hajage
##' @keywords internal
asciiOpts <- function(select = "all", .backends = NULL, .outputs = NULL, .extensions = NULL, .cli = NULL, .args = NULL, .O = NULL, .f = NULL, .e = NULL, .preambule = NULL) {

  if (is.null(.backends)) {
    .backends <- c("asciidoc", "a2x", "t2t", "pandoc", "markdown2pdf")
  }

  if (is.null(.outputs)) {
    .outputs = list(
      asciidoc = c("html", "docbook", "slidy"),
      a2x = c("xhtml", "chunked", "htmlhelp", "epub", "pdf", "text", "dvi", "ps", "tex", "asciidoc"),
      t2t = c("aap", "aas", "aat", "adoc", "bbcode", "creole", "csv", "dbk", "doku", "gwiki", "html", "html5", "lout", "man", "md", "mgp", "moin", "ods", "pm6", "pmw", "red", "rtf", "sgml", "spip", "tex", "txt", "txt2t", "wiki", "xhtml", "xhtmls"),
      pandoc = c("native", "json", "html", "html+lhs", "s5", "slidy", "docbook", "opendocument", "latex", "latex+lhs", "context", "texinfo", "man", "markdown", "markdown+lhs", "plain", "rst", "rst+lhs", "mediawiki", "textile", "rtf", "org", "odt", "epub"),
      markdown2pdf = "")
  }

  if (is.null(.extensions)) {
    .extensions = list(
      docbook = "xml",
      slidy = "html",
      xhtml = "html",
      chunked = "html",
      htmlhelp = "html",
      asciidoc = "adoc",
      html5 = "html",
      txt2t = "t2t",
      xhtml = "html",
      xhtmls = "html",
      "html+lhs" = "html",
      s5 = "html",
      slidy = "html",
      opendocument = "odt",
      latex = "tex",
      "latex+lhs" = "tex",
      markdown = "md",
      "markdown+lhs" = "md",
      "rst+lhs" = "rst")
  }

  if (is.null(.cli)) {
    .cli = list(
      asciidoc = c("asciidoc %options", paste(Sys.getenv("COMSPEC"), "/c", "asciidoc.py %options"), "bash -c \"asciidoc %options \""),
      a2x = c("a2x %options", paste(Sys.getenv("COMSPEC"), "/c", "a2x.py %options"), "bash -c \"a2x %options \""),
      t2t = c("txt2tags %options", paste(Sys.getenv("COMSPEC"), "/c", "txt2tags.py %options"), "bash -c \"txt2tags %options \""),
      pandoc = c("pandoc %options", paste(Sys.getenv("COMSPEC"), "/c", "pandoc %options"), "bash -c \"pandoc %options \""),
      markdown2pdf = c("markdown2pdf %options", paste(Sys.getenv("COMSPEC"), "/c", "markdown2pdf %options"), "bash -c \"markdown2pdf %options \""))
  }

  if (is.null(.args)) {
    .args = list(
      asciidoc = "-a encoding=%e -b %f %O -o %d/%o %i",
      a2x = "-a encoding=%e -D %d -f %f %O %i",
      t2t = "--encoding=%e -t %f %O -o %d/%o %i",
      pandoc = "-t %f -o %d/%o %O %i",
      markdown2pdf = "-o %d/%o %O %i")
  }

  if (is.null(.O)) {
    .O = list(
      asciidoc = "-a toc",
      a2x = "-a toc",
      t2t = "",
      pandoc = "-s",
      markdown2pdf = "")
  }

  if (is.null(.f)) {
    .f = list(
      asciidoc = "html",
      a2x = "xhtml",
      t2t = "html",
      pandoc = "html",
      markdown2pdf = "pdf")
  }

  if (is.null(.e)) {
    .e = list(
      asciidoc = "UTF-8",
      a2x = "UTF-8",
      t2t = "UTF-8",
      pandoc = "",
      markdown2pdf = "")
  }

  if (is.null(.preambule)) {
    .preambule = list(
      asciidoc =
"= %title
:author:    %author
:email:     %email
:revdate:   %date

",
      t2t =
"%title
%author
%date

",
      pandoc =
"% %title
% %author %email
% %date

")
  }

  opts <- list(".backends" = .backends, ".outputs" = .outputs, ".extensions" = .extensions, ".cli" = .cli, ".args" = .args, ".O" = .O, ".f" = .f, ".e" = .e, ".preambule" = .preambule)

  if (select == "all")
    return(opts)
  else
    return(opts[[select]])
}

##' replace
##'
##' @param backend backend
##' @param cygwin cygwin
##' @param i i
##' @param f f
##' @param d d
##' @param e e
##' @param O O
##' @keywords internal
##' @author David Hajage
replace <- function(backend = getOption("asciiBackend"), cygwin = FALSE, i, f = NULL, d = NULL, e = NULL, O = NULL) {
  if (is.null(f))
    f <- asciiOpts(".f")[[backend]]
  if (is.null(d))
    d <- dirname(i)
  if (is.null(e))
    e <- asciiOpts(".e")[[backend]]
  if (is.null(O))
    O <- asciiOpts(".O")[[backend]]

  extension <- ifelse(f %in% names(asciiOpts(".extensions")), asciiOpts(".extensions")[[f]], f)
  basefile <- sub("(.+)(\\..+$)", "\\1", basename(i))
  file <- paste(basefile, extension, sep = ".")

  if (grepl("w|W", .Platform$OS.type)) {
    if (cygwin) {
      cli <- asciiOpts(".cli")[[backend]][3]
    } else {
      cli <- asciiOpts(".cli")[[backend]][2]
    }
  } else {
    cli <- asciiOpts(".cli")[[backend]][1]
  }

  args <- asciiOpts(".args")[[backend]]
  args <- sub("%i", i, args)
  args <- sub("%f", f, args)
  args <- sub("%d", d, args)
  args <- sub("%o", file, args)
  args <- sub("%e", e, args)
  args <- sub("%O", O, args)

  results <- sub("%options", args, cli)
  attr(results, "file") <- file
  attr(results, "directory") <- d
  attr(results, "f") <- f
  attr(results, "cygwin") <- cygwin
  results
}

##' Convert a file with specified backend
##'
##' This function convert a file with asciidoc, txt2tags or pandoc backend
##' @title Convert a file with specified backend
##' @param i input file
##' @param d output directory
##' @param f format
##' @param e encoding
##' @param O other options
##' @param backend backend (\code{"asciidoc"}, \code{"t2t"} or \code{"pandoc"})
##' @param cygwin use cygwin?
##' @param open open resulting file?
##' @return Nothing
##' @export
##' @author David Hajage
convert <- function(i, d = NULL, f = NULL, e = NULL, O = NULL, backend = getOption("asciiBackend"), cygwin = FALSE, open = FALSE) {

  if (!(backend %in% asciiOpts(".backends")))
    stop(paste("Wrong backend. Please choose: ", paste(asciiOpts(".backends"), collapse = ", "), ".", sep = ""))

  cmd <- replace(backend, cygwin = cygwin, i = i, d = d, f = f, e = e, O = O)
  windows <- grepl("w|W", .Platform$OS.type)
  if (windows) { # because pandoc doesn't like path with "/"
    cmd <- gsub("/", "\\\\", cmd)
    cmd <- sub("cmd.exe \\\\c ", "cmd.exe /c ", cmd)
  }
  err <- system(cmd, wait = TRUE)

  file <- attr(cmd, "file")
  f <- attr(cmd, "f")
  directory <- attr(cmd, "d")

  if (f == "chunked") {
    dfile <- paste(directory, paste(sub("(.+)(\\..+$)", "\\1", file), "chunked", sep = "."), "index.html", sep = "/")
  } else {
    dfile <- paste(directory, file, sep = "/")
  }

  if (open) {
    cat("Trying to open ", dfile, sep = "")
    if (windows) {
      cat(" with shell.exec...\n")
      if (!grepl(":", dfile)) {
        dfile <- paste("\"", sub("(^.+)(/$)", "\\", getwd()), "/", dfile, "\"", sep = "")
      }
      shell.exec(dfile)
    } else if (grepl("darwin", version$os)) {
      cat(" with open...\n")
      system(paste(shQuote("open"), shQuote(dfile)), wait = FALSE, ignore.stderr = TRUE)
    } else {
      cat(" with xdg-open...\n")
      system(paste(shQuote("/usr/bin/xdg-open"), shQuote(dfile)), wait = FALSE, ignore.stderr = TRUE)
    }
  }
  invisible(cmd)
}

##' Create a section
##'
##' \code{section} can be used with \code{export} function to add\dots a section
##' @param caption a string
##' @param caption.level caption level
##' @return A section object.
##' @export
##' @author David Hajage
section <- function(caption, caption.level = 1) {
  results <- list(caption = caption, caption.level = caption.level)
  class(results) <- "section"
  results
}

##' Print a section object
##'
##' Print a section object
##' @param x a section object
##' @param backend ascii backend
##' @param ... not used
##' @method print section
##' @export
##' @author David Hajage
print.section <- function(x, backend = getOption("asciiBackend"), ...) {
  caption <- x$caption
  caption.level <- x$caption.level
  if (backend == "asciidoc" | backend == "a2x")
    results <- header.asciidoc(caption, caption.level+1)
  if (backend == "t2t")
    results <- header.t2t(caption, caption.level)
  if (backend == "pandoc" | backend == "markdown2pdf")
    results <- header.pandoc(caption, caption.level)
  cat("\n", results, sep = "")
}

##' Create a paragraph
##'
##' \code{paragraph} can be used with \code{export} function to add\dots a paragraph
##' @param ... strings composing the paragraph
##' @param new whether to create a new paragraph or to continue a preceding one
##' @return A paragraph object.
##' @export
##' @author David Hajage
paragraph <- function(..., new = TRUE) {
  results <- list(text = list(...), new = new)
  class(results) <- "paragraph"
  results
}

##' Print a paragraph object
##'
##' Print a paragraph object
##' @param x a paragraph object
##' @param ... not used
##' @method print paragraph
##' @export
##' @author David Hajage
print.paragraph <- function(x, ...) {
  newline <- "\n"
  text <- x$text
  new <- x$new

  if (new) {
    cat(newline)
  }
  for (i in 1:length(text)) {
    if (class(text[[i]]) == "sexpr") {
      print(text[[i]])
    }
    else {
      cat(text[[i]], " ", sep = "")
    }
  }
  cat("\n")
}

##' Create a verbatim paragraph
##'
##' \code{verbatim} can be used with \code{export} function to add a verbatim paragraph
##' @param ... strings composing the paragraph (line by line)
##' @return A verbatim object.
##' @export
##' @author David Hajage
verbatim <- function(...) {
  results <- c(...)
  class(results) <- "verbatim"
  results
}

##' Print a verbatim object
##'
##' Print a verbatim object
##' @param x a verbatim object
##' @param backend ascii backend
##' @param ... not used
##' @method print verbatim
##' @export
##' @author David Hajage
print.verbatim <- function(x, backend = getOption("asciiBackend"), ...) {
  cat("\n")
  if (backend == "asciidoc" | backend == "a2x")
    cat("----\n")
  if (backend == "t2t")
    cat("```\n")
  if (backend == "pandoc" | backend == "markdown2pdf")
    cat("\n~~~~~~~{.R}\n")

  cat(x, sep = "\n", ...)

  if (backend == "asciidoc" | backend == "a2x")
    cat("----\n\n")
  if (backend == "t2t")
    cat("```\n\n")
  if (backend == "pandoc" | backend == "markdown2pdf")
    cat("~~~~~~~~~~~\n\n")
}

##' Insert an inline R result
##'
##' \code{sexpr} can be used with \code{export} function to insert an inline R results
##' @param x an R results (of length one)
##' @return A sexpr object.
##' @export
##' @author David Hajage
sexpr <- function(x) {
  results <- x
  class(results) <- "sexpr"
  results
}

##' Print a sexpr object
##'
##' Print a sexpr object
##' @param x a sexpr object
##' @param ... not used
##' @method print sexpr
##' @export
##' @author David Hajage
print.sexpr <- function(x, ...) {
  cat(x, sep = "")
}

##' Export R objects
##'
##' \code{out} can be used with \code{export} function to insert an R results
##' @param x an R object
##' @param results if \code{'verbatim'}, the output is included in a verbatim environment. If \code{'ascii'}, the output is taken to be already proper markup and included as is.
##' @return An out object
##' @export
##' @author David Hajage
out <- function(x, results = "verbatim") {
  results <- list(x, results)
  class(results) <- "out"
  results
}

##' Print an out object
##'
##' Print an out object
##' @param x an out object
##' @param backend ascii backend
##' @param ... not used
##' @method print out
##' @export
##' @author David Hajage
print.out <- function(x, backend = getOption("asciiBackend"), ...) {
  results <- x[[2]]
  cat("\n")
  if (results == "verbatim") {
    if (backend == "asciidoc" | backend == "a2x")
      cat("----\n")
    if (backend == "t2t")
      cat("```\n")
    if (backend == "pandoc" | backend == "markdown2pdf")
      cat("\n~~~~~~~{.R}\n")
  }
  print(x[[1]], ...)
  if (results == "verbatim") {
    if (backend == "asciidoc" | backend == "a2x")
      cat("----\n\n")
    if (backend == "t2t")
      cat("```\n\n")
    if (backend == "pandoc" | backend == "markdown2pdf")
      cat("~~~~~~~~~~~\n\n")
  } else {
    cat("\n")
  }
}

##' Insert figure
##'
##' \code{graph} can be used with \code{export} function to insert an R graphic.
##' @aliases graph
##' @param file character string (
##' @param graph a recordedplot, a lattice plot, a ggplot, or an expression producing a plot (optional if the file already exists)
##' @param format jpg, png or pdf (or guessed with the file name)
##' @param ... additional arguments (passed to format options)
##' @return A fig object
##' @export
##' @author David Hajage
fig <- function(file = NULL, graph = NULL, format = NULL, ...) {

  if (is.null(file) & is.null(graph)) {
    stop("Please provide a graph or a link to an existing graph.")
  }

  if (is.null(format)) {
    if (!is.null(file)) {
      ext <- tolower(sub("(^.+\\.)(.+$)", "\\2", file))
      if (grepl(ext, "jpg|jpeg"))
        format <- "jpg"
      if (grepl(ext, "png"))
        format <- "png"
      if (grepl(ext, "pdf"))
        format <- "pdf"
    } else {
      format <- "png"
    }
  }

  if (is.null(file)) {
    file <- paste(tempfile("graph"), format, sep = ".")
  }

  if (!is.null(graph)) {
    if (format == "jpg") {
      jpeg(file, ...)
    }
    if (format == "png") {
      png(file, ...)
    }
    if (format == "pdf") {
      pdf(file, ...)
    }
    if (is.expression(graph)) {
      eval(graph)
    }
    else {
      print(graph)
    }
    dev.off()
  }

  results <- file
  class(results) <- "fig"
  results
}

graph <- fig

##' Print an graph object
##'
##' Print an graph object
##' @param x an graph object
##' @param backend ascii backend
##' @param ... not used
##' @method print fig
##' @export
##' @author David Hajage
print.fig <- function(x, backend = getOption("asciiBackend"), ...) {
  if (backend == "asciidoc" | backend == "a2x")
    results <- paste("image::", x, "[]", sep = "")
  if (backend == "t2t")
    results <- paste("[", x, "]", sep = "")
  if (backend == "pandoc" | backend == "markdown2pdf")
    results <- paste("![](", x, ")", sep = "")
  cat("\n", results, "\n\n", sep = "")
}


##' Produce a report
##'
##' Produce a report from a list of R objects. This function can be
##' used directly, or through a \code{Report} object (see
##' examples). \code{Report$new()} creates a new object,
##' \code{Report$create()} produce a report. Exportation options can
##' be specified with \code{Report$nameoftheoption <- option} or
##' directly in \code{Report$create(nameoftheoption = option)}.
##'
##' Special objects can be used to create sections (see
##' \code{?section}), paragraphs (see \code{?paragraph}), verbatim
##' environment (see \code{?verbatim} and to insert figures (see
##' \code{?fig}) or inline results (see \code{?sexpr}). Helpers exist:
##' \code{Report$addSection()}, \code{Report$addParagraph()},
##' \code{Report$addVerbatim()}, \code{Report$addFig()}.
##'
##' It needs a working installation of asciidoc, a2x tool chain,
##' txt2tags, pandoc and/or markdown2pdf.
##'
##' @aliases Report
##' @title Report creation
##' @param ... R objects (not used if \code{"list"} is not NULL)
##' @param list list of R objects
##' @param file name of the output file (without extension)
##' @param format format of the output file
##' @param open open resulting file?
##' @param backend backend
##' @param encoding encoding
##' @param options other options
##' @param cygwin use cygwin?
##' @param title title of the report
##' @param author author of the report
##' @param email email of the author
##' @param date date
##' @return Nothing
##' @export
##' @rdname createreport
##' @author David Hajage
##' @examples
##' \dontrun{
##' options(asciiType = "asciidoc")
##' createreport(head(esoph))
##'
##' r <- Report$new(author = "David Hajage", email = "dhajage at gmail dot com")
##' r$add(section("First section"))
##' r$addSection("First subsection", 2)
##' r$add(paragraph("The data set has", sexpr(nrow(esoph)), " lines. See yourself:"), esoph)
##' r$addSection("Second subsection: age and alc group", 2)
##' tab <- with(esoph, table(alcgp, agegp))
##' r$add(ascii(tab), ascii(summary(tab), format = "nice"))
##' r$create()
##' r$format <- "slidy"
##' r$createt()
##'
##' r$title <- "R report example"
##' r$author <- "David Hajage"
##' r$email <- "dhajage at gmail dot com"
##' options(asciiType = "pandoc")
##' r$backend <- "pandoc"
##' r$format <- "odt"
##' r$create()
##'
##' r$create(backend = "markdown2pdf", format = "pdf")
##' }
createreport <- function(..., list = NULL, file = NULL, format = NULL, open = TRUE, backend = getOption("asciiBackend"), encoding = NULL, options = NULL, cygwin = FALSE, title = NULL, author = NULL, email = NULL, date = NULL) {

  if (is.null(file)) {
    file <- tempfile("R-report")
  }

  wd <- dirname(file)
  file <- paste(wd, basename(file), sep = "/")

  if (is.null(title)) {
    title <- sub("\\~", "\\\\~", paste(wd, basename(file), sep = "/"))
    if (backend == "asciidoc")
      title <- beauty.asciidoc(title, "m")
    if (backend == "text2tags")
      title <- beauty.t2t(title, "m")
    if (backend == "pandoc")
      title <- beauty.pandoc(title, "m")
  }

  if (backend == "a2x") {
    preambule <- asciiOpts(".preambule")[["asciidoc"]]
  } else if (backend == "markdown2pdf") {
    preambule <- asciiOpts(".preambule")[["pandoc"]]
  } else {
    preambule <- asciiOpts(".preambule")[[backend]]
  }

  if (is.null(author)) {
    author <- paste(version$language, " ", paste(version$major, version$minor, sep = "."), " ascii ", packageDescription("ascii")$Version, sep = "")
  }

  if (is.null(email)) {
    email <- "cran@r-project.org"
  }

  if (is.null(date)) {
    date <- format(Sys.time(), "%Y/%m/%d %X")
  }

  preambule <- sub("%title", title, preambule)
  preambule <- sub("%author", author, preambule)
  preambule <- sub("%email", email, preambule)
  preambule <- sub("%date", date, preambule)

  if (is.null(list)) {
    args <- list(...)
  } else {
    args <- list
  }

  lines <- capture.output({
    cat(preambule)

    for (i in seq_along(args)) {
      arg <- args[[i]]
      if (!is.null(names(args))) {
        if (names(args)[i] != "") {
          cat(".", names(args)[i], "\n", sep = "")
        }
      }
      if ("asciiTable" %in% class(arg) | "asciiList" %in% class(arg) | "asciiMixed" %in% class(arg)) {
        cat("\n")
        print(arg)
        cat("\n")
      } else if ("out" %in% class(arg) | "section" %in% class(arg) | "paragraph" %in% class(arg) | "verbatim" %in% class(arg) | "fig" %in% class(arg)) {
        print(arg, backend = backend)
      } else {
        print(out(arg, "verbatim"), backend = backend)
      }
    }})
  textfile <- paste(file, "txt", sep = ".")
  f <- file(textfile, "w")
  writeLines(lines, f)
  close(f)

  convert(i = textfile, d = NULL, f = format, e = encoding, O = options, backend = backend, cygwin = cygwin, open = open)
}

##' Report generator
##'
##' @author David Hajage
##' @rdname createreport
##' @export
Report <- setRefClass("Report",
                      fields = c("file", "format", "open", "backend", "encoding", "options", "cygwin", "title", "author", "email", "date", "objects"),

                      methods = list(
                        initialize = function(...) {
                          initFields(file = NULL)
                          initFields(format = "html")
                          initFields(open = TRUE)
                          initFields(backend = getOption("asciiBackend"))
                          initFields(encoding = NULL)
                          initFields(options = NULL)
                          initFields(cygwin = FALSE)
                          initFields(title = NULL)
                          initFields(author = NULL)
                          initFields(email = NULL)
                          initFields(date = NULL)
                          initFields(objects = list())

                          callSuper(...)
                        },

                        add = function(...) {
                          obj <- list(...)
                          .self$objects <- c(.self$objects, obj)
                        },

                        addSection = function(x, caption.level = 1) {
                          .self$objects <- c(.self$objects, list(section(x, caption.level)))
                        },

                        addParagraphs = function(..., new = TRUE) {
                          .self$objects <- c(.self$objects, list(paragraph(..., new = new)))
                        },

                        addVerbatim = function(...) {
                          .self$objects <- c(.self$objects, list(verbatim(...)))
                        },

                        addFig = function(file = NULL, graph = NULL, format = NULL, ...) {
                          .self$objects <- c(.self$objects, list(fig(file, graph, format, ...)))
                        },

                        show.Report = function(help = FALSE) {
                          cat("Report object:\n\n")
                          cat("  title:   ", ifelse(is.null(.self$title), "None", .self$title), "\n")
                          cat("  author:  ", ifelse(is.null(.self$author), "None", .self$author), "\n")
                          cat("  email:   ", ifelse(is.null(.self$email), "None", .self$email), "\n")
                          cat("  date:    ", ifelse(is.null(.self$date), format(Sys.time(), "%Y/%m/%d %X"), .self$date), "\n")
                          cat("  file:    ", ifelse(is.null(.self$file), "Temporary file", .self$file), "\n")
                          cat("  open:    ", .self$open, "\n")
                          cat("  backend: ", .self$backend, "\n")
                          cat("  format:  ", .self$format, "\n")
                          cat("  encoding:", ifelse(is.null(.self$encoding), asciiOpts(".e")[[.self$backend]], .self$encoding), "\n")
                          cat("  options: ", ifelse(is.null(.self$options), asciiOpts(".O")[[.self$backend]], .self$options), "\n")
                          cat("  cygwin:  ", .self$cygwin, "\n")

                          if(help) {
                            cat("\nTo change a slot:\n")
                            cat("\tyourreport$slot <- 'value'\n\n")
                            cat("To create the report:\n")
                            cat("\tyourreport$create()\n")
                          }
                        },

                        create = function(list = .self$objects, file = .self$file, format = .self$format, open = .self$open, backend = .self$backend, encoding = .self$encoding, options = .self$options, cygwin = .self$cygwin, title = .self$title, author = .self$author, email = .self$email, date = .self$date) {
                          createreport(list = list, file = file, format = format, open = open, backend = backend, encoding = encoding, options = options, cygwin = cygwin, title = title, author = author, email = email, date = date)
                        }
                        )
                      )
