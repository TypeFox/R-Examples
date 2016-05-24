## NOTE: needs commands "convert" from ImageMagick (http://www.imagemagick.org/)
tex2image <- function(tex, format = "png", width = 6, 
  pt = 12, density = 350, dir = NULL, tdir = NULL, idir = NULL,
  width.border = 0L, col.border = "white", resize = 650, shave = 2,
  packages = c("amsmath", "amssymb", "amsfonts"),
  header = c("\\setlength{\\parindent}{0pt}",
    "\\renewcommand{\\sfdefault}{phv}",
    "\\IfFileExists{sfmath.sty}{\n\\RequirePackage{sfmath}\n\\renewcommand{\\rmdefault}{phv}}{}"),
  header2 = NULL, Sweave = TRUE, show = TRUE, name = "tex2image")
{
  ## directory handling
  if(is.null(tdir)) {
    tdir <- tempfile()
    on.exit(unlink(tdir))
  }
  tdir <- file.path(path.expand(tdir), "tex2image")
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

  if(!is.list(tex) && (length(text) < 2L) && file.exists(tex)) {
    texfile <- file_path_as_absolute(tex)
    tex <- readLines(con = texfile)
    texdir <- dirname(texfile)
    cfiles <- list.files(texdir)
    cfiles <- cfiles[cfiles != "tex2image"]
    cfiles <- cfiles[cfiles != basename(texfile)]
    if(length(pdfs <- grep("pdf", file_ext(cfiles))))
      cfiles <- cfiles[-pdfs]
    file.copy(file.path(texdir, cfiles), file.path(tdir, cfiles))
    texfile <- paste("tex2image-", basename(texfile), sep = "")
    name <- file_path_sans_ext(texfile)
  } else texdir <- tempdir()

  if(any(grepl("\\documentclass", tex, fixed = TRUE))) {
    begin <- grep("\\begin{document}", tex, fixed = TRUE)
    end <- grep("\\end{document}", tex, fixed = TRUE)
    tex <- tex[(begin + 1):(end - 1)]
  }

  if(is.null(dir)) dir <- texdir
  if(dir == ".") dir <- getwd()

  owd <- getwd()
  setwd(tdir)
  on.exit(setwd(owd), add = TRUE)

  ## output formats
  format <- tolower(format)

  ## LaTeX packages  
  packages <- unique(c(packages, c("a4wide", "graphicx", "url", "color", "realboxes")))

  if(length(graphics <- grep("includegraphics", unlist(tex), fixed = TRUE, value = TRUE))) {
    if(is.null(idir))
      idir <- texdir
    idir <- path.expand(idir)
    files <- list.files(idir)
    cp <- NULL
    for(k in 1L:length(graphics)) {
      graphics[k] <- extract_command(graphics[k], "includegraphics")
      cp <- c(cp, grep(graphics[k], files, fixed = TRUE, value = TRUE))
    }
    if(length(cp)) {
      for(f in cp) 
        file.copy(from = file.path(idir, f), to = file.path(tdir, f), overwrite = TRUE)
    } else stop(paste("graphic is missing in ", texdir, "!", sep = ""))
  }
  
  ## auxiliary LaTeX file
  texlines <- paste("\\documentclass[a4paper,", pt, "pt]{article}", sep = "")
  for(i in packages) {
    brackets <- if(grepl("{", i, fixed = TRUE)) NULL else c("{", "}")
    texlines <- c(texlines, paste("\\usepackage", brackets[1], i, brackets[2], sep = ""))
  }
  if(Sweave) texlines <- c(texlines, paste("\\usepackage{",
    file.path(R.home("share"), "texmf", "tex", "latex", "Sweave"), "}", sep = ""))
  texlines <- c(
    texlines,
    "\\pagestyle{empty}"
  )
  texlines <- c(texlines, header)
  texlines <- c(texlines, paste("\\setlength{\\textwidth}{", width, "in}", sep = ""))
  texlines <- c(texlines, "\\begin{document}")
  texlines <- c(texlines, header2)
  tex <- if(!is.list(tex)) list(tex) else tex
  nt <- length(tex)
  pic_names <- if(is.null(names(tex))) {
    if(nt > 1) paste(name, 1:length(tex), sep = "_") else name
  } else {
    paste(name, names(tex), sep = "_")
  }
  pic_names <- paste(pic_names, format, sep = ".")
  for(i in 1:nt) {
    if(!any(grepl("begin{figure}", tex[[i]], fixed = TRUE)) &&
      !any(grepl("caption{", tex[[i]], fixed = TRUE))) {
      texlines <- c(texlines, paste("\\Fbox{\\begin{minipage}[t]{", width, "in}", sep = ""))
    }
    texlines <- c(texlines, tex[[i]])
    if(!any(grepl("begin{figure}", tex[[i]], fixed = TRUE)) &&
      !any(grepl("caption{", tex[[i]], fixed = TRUE))) {
      texlines <- c(texlines, "\\end{minipage}}")
    }
    if(nt > 1)
      texlines <- c(texlines, "\\newpage")
  }
  texlines <- c(texlines, "\\end{document}")
  file.create(paste(tdir, "/", name, ".log", sep = ""))
  writeLines(text = texlines, con = paste(tdir, "/", name, ".tex", sep = ""))

  ## compile LaTeX into PDF
  tools::texi2dvi(file = paste(name, ".tex", sep = ""), pdf = TRUE, clean = TRUE, quiet = TRUE)

  ## shell command on Windows
  shcmd <- Sys.getenv("COMSPEC")
  if(shcmd != "") shcmd <- paste(shQuote(shcmd), "/c")

  ## convert to images
  image <- paste(name, if(nt > 1) 1:nt else NULL, ".", format, sep = "")
  dirout <- rep(NA, length(name))
  for(i in 1:nt) {
     if(!(format %in% c("pdf", "svg"))) {
      if(format == "png") {
        cmd <- paste("convert -trim -shave ", shave, "x", shave," -density ", density, " ",
          name, ".pdf[", i - 1, "] -transparent white ", image[i], " > ", name, i, ".log", sep = "")
      } else {
        cmd <- paste("convert -trim -shave ", shave, "x", shave," -density ", density, " ",
          name, ".pdf[", i - 1, "] ", image[i], " > ", name, i, ".log", sep = "")
      }
      system(paste(shcmd, cmd))
      if(!is.null(resize)) {
        cmd <- paste("convert -resize ", resize, "x ", image[i], " ", image[i], " > ", name[i], ".log", sep = "")
        system(paste(shcmd, cmd))
      } else resize <- 800
      width.border <- as.integer(width.border)
      if(width.border > 0L) {
        width.border <- paste(width.border, "x", width.border, sep = "")
        cmd <- paste("convert ", image[i], " -bordercolor ", col.border, " -border ", width.border, " ",
          image[i], " > ", name[i], ".log", sep = "")
        system(paste(shcmd, cmd))
      }
    } else {
      system(paste(shcmd, "pdfcrop --margins", -shave, "--clip", paste0(name, ".pdf"), paste0(name, ".pdf")), ignore.stdout = TRUE)
      if(format == "svg") system(paste(shcmd, "pdf2svg", paste0(name, ".pdf"), paste0(name, if(nt > 1) "_%d" else NULL, ".svg"), "all"), ignore.stdout = TRUE)
    }
    dirout[i] <- file.path(path.expand(dir), pic_names[i])
    file.copy(from = file.path(tdir, image[i]), 
      to = dirout[i], overwrite = TRUE)
    dirout[i] <- normalizePath(dirout[i])
    if(show) browseFile(dirout[i])
  }
  
  return(invisible(dirout))
}
