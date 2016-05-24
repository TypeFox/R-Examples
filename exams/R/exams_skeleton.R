exams_skeleton <- exams.skeleton <- function(dir = ".",
  type = c("num", "schoice", "mchoice", "string", "cloze"),
  writer = c("exams2html", "exams2pdf", "exams2moodle", "exams2qti12", "exams2qti21", "exams2arsnova", "exams2nops"),
  markup = "latex", absolute = FALSE, encoding = "")
{
  ## match available types/writers
  type <- as.vector(sapply(type, match.arg,
    c("num", "schoice", "mchoice", "cloze", "string")))
  writer <- as.vector(sapply(writer, match.arg,
    c("exams2html", "exams2pdf", "exams2moodle", "exams2qti12", "exams2qti21", "exams2arsnova", "exams2nops")))
  markup <- match.arg(markup, c("latex", "markdown"))

  ## create output directory (if necessary)
  create_dir <- function(path) {
    isdir <- file.info(path)$isdir
    if(identical(isdir, FALSE)) stop(sprintf("file (rather than directory) %s exists", path))
    if(is.na(isdir)) dir.create(path)
  }
  create_dir(dir)
  create_dir(edir <- file.path(dir, "exercises"))
  if(any(c("exams2html", "exams2pdf", "exams2qti12", "exams2qti21") %in% writer)) create_dir(templ <- file.path(dir, "templates"))
  if("exams2nops" %in% writer) create_dir(file.path(dir, "nops"))
  pdir <- find.package("exams")
  if(absolute) dir <- file_path_as_absolute(dir)
  
  ## select exercises for demo script and all available exercises
  exrc <- c(
    "num"     = "tstat",
    "schoice" = "tstat2",
    "mchoice" = "boxplots",
    "cloze"   = "boxhist",
    "string"  = "function"
  )
  exrc <- exrc[type]
  exrc <- paste0(exrc, if(markup == "latex") ".Rnw" else ".Rmd")
  axrc <- list.files(path = file.path(pdir, "exercises"), pattern = if(markup == "latex") "Rnw$" else "Rmd$")
  axrc <- axrc[axrc != "confint.Rnw"]
  file.copy(file.path(pdir, "exercises", axrc), file.path(edir, axrc))

  ## encoding
  enc <- gsub("-", "", tolower(encoding), fixed = TRUE)
  if(enc %in% c("iso8859", "iso88591")) enc <- "latin1"
  if(enc == "iso885915") enc <- "latin9"
  charset <- encoding

  if((markup == "markdown" || "exams2arsnova" %in% writer)) {
    if(!(enc %in% c("", "utf8"))) {
      warning("pandoc-based conversion needs UTF-8 encoding")
      encoding <- "UTF-8"
    }
    enc <- "utf8"
  }
  if(enc == "utf8") {
    exrc_enc <- "currency8.Rnw"
    if(markup == "latex" | encoding != "") exrc <- c(exrc, exrc_enc)
    charset <- "UTF-8"
  }
  if(enc == "latin1") {
    exrc_enc <- "currency1.Rnw"
    if(markup == "latex") exrc <- c(exrc, exrc_enc)
    charset <- "ISO-8859-1"
  }
  if(enc == "latin9") {
    exrc_enc <- "currency9.Rnw"
    if(markup == "latex") exrc <- c(exrc, exrc_enc)
    charset <- "ISO-8859-15"
  }
  
  ## copy templates
  if("exams2pdf" %in% writer) {
    file.copy(file.path(pdir, "tex", c("exam.tex", "solution.tex")), templ)
    if(encoding != "" | markup == "markdown") {
      x <- readLines(file.path(pdir, "tex", "plain.tex"))
      i <- grep("documentclass", x, fixed = TRUE)[1L]
      x <- c(x[1L:i],"", sprintf('\\usepackage[%s]{inputenc}', enc), "", x[-(1L:i)])
      if(markup == "markdown") x <- gsub(",amsmath", ",amsmath,hyperref", x, fixed = TRUE)
      writeLines(x, file.path(templ, "plain.tex"))      
      x <- readLines(file.path(pdir, "tex", "exam.tex"))
      i <- grep("documentclass", x, fixed = TRUE)[1L]
      x <- c(x[1L:i], "", sprintf('\\usepackage[%s]{inputenc}', enc), "", x[-(1L:i)])
      if(markup == "markdown") x <- gsub(",amsmath", ",amsmath,hyperref", x, fixed = TRUE)
      writeLines(x, file.path(templ, "exam.tex"))
      x <- readLines(file.path(pdir, "tex", "solution.tex"))
      i <- grep("documentclass", x, fixed = TRUE)[1L]
      x <- c(x[1L:i], "", sprintf('\\usepackage[%s]{inputenc}', enc), "", x[-(1L:i)])
      if(markup == "markdown") x <- gsub(",amsmath", ",amsmath,hyperref", x, fixed = TRUE)
      writeLines(x, file.path(templ, "solution.tex"))      
    }
  }
  if("exams2html" %in% writer) {
    file.copy(file.path(pdir, "xml", "plain.html"), templ)
    if(encoding != "" | markup == "markdown") {
      x <- readLines(file.path(pdir, "xml", "plain.html"))
      i <- grep("</head>", x, fixed = TRUE)[1L] - 1L
      x <- c(x[1L:i], "", sprintf('<meta charset="%s">', charset), "", x[-(1L:i)])
      writeLines(x, file.path(templ, "plain.html"))
    }
  }
  if("exams2qti12" %in% writer) {
    file.copy(file.path(pdir, "xml", "qti12.xml"), templ)
    if(encoding != "") {
      x <- readLines(file.path(pdir, "xml", "qti12.xml"))
      x[1L] <- gsub("UTF-8", charset, x[1L], fixed = TRUE)
      writeLines(x, file.path(templ, "qti12.xml"))
    }
  }

  demo_all <- c(
    '## exams ----------------------------------------------------------------------------',
    '',
    '## load package',
    'library("exams")',
    '',
    '## this script gives an overview of the example exercises provided',
    '## and basic usage of exams2html/exams2pdf - for more advanced usage',
    '## and further interfaces, see the other demo-*.R scripts',
    '',
    '## to get an overview of the available exercises in this demo,',
    '## switch to the "exercises" folder',
    sprintf('setwd("%s")', if(absolute) file.path(dir, "exercises") else "exercises"),
    'dir()',
    '',
    sprintf('## in the following the exercises in %s (.Rnw) format are discussed',
      if(markup == "latex") "R/LaTeX" else "R/Markdown"),
    '',
    '',
    '## inspect individual exercises -----------------------------------------------------',
    '',
    '## simply turn a single exercise into a HTML file (shown in browser)',
    sprintf('exams2html("%s")', exrc[1L]),
    '## or a PDF file (shown in PDF viewer)',
    sprintf('exams2pdf("%s")', exrc[1L]),
    '',
    '## extract the meta-information to check whether it is processed correctly',
    sprintf('exams_metainfo(exams2html("%s"))', exrc[1L]),
    '',
    '',
    '## types of exercises ---------------------------------------------------------------',
    '',
    if(!("num" %in% type)) NULL else c(
    '## numeric exercises',
    'exams2html(c(',
    '  "tstat.Rnw",      ## computation of t-statistic',
    '  "dist.Rnw",       ## computation of Euclidean distance',
    '  "regression.Rnw", ## prediction in simple linear regression',
    '  "lagrange.Rnw"    ## optimization under constraint',
    '))',
    ''),
    if(!("schoice" %in% type)) NULL else c(
    '## single choice exercises',
    'exams2html(c(',
    '  "tstat2.Rnw", ## single-choice list of numeric tstat exercise (by hand)',
    '  "dist3.Rnw"   ## single-choice list of numeric dist exercise (via num_to_schoice)',
    '))',
    ''),
    if(!("mchoice" %in% type)) NULL else c(
    '## multiple choice exercises',
    'exams2html(c(',
    '  "anova.Rnw",       ## interpret ANOVA results',
    '  "boxplots.Rnw",    ## interpret parallel boxplots',
    '  "scatterplot.Rnw", ## interpret scatterplot',
    '  "ttest.Rnw",       ## 2-sample t-test',
    '  "relfreq.Rnw",     ## table of relative frequencies',
    '  "cholesky.Rnw"     ## Cholesky factorization',
    '))',
    ''),
    if(!("string" %in% type)) NULL else c(
    '## string exercises',
    'exams2html(c(',
    '  "function.Rnw", ## names of R functions',
    '  "countrycodes.Rnw"   ## three-letter ISO country codes',
    '))',
    ''),
    if(!("cloze" %in% type)) NULL else c(
    '## cloze exercises (combining several num/schoice/mchoice/string parts)',
    'exams2html(c(',
    '  "boxhist.Rnw",  ## download and describe artificial data (schoice/num)',
    '  "confint2.Rnw", ## compute a confidence interval (num)',
    '  "dist2.Rnw",    ## compute several types of distances (num)',
    '  "fourfold.Rnw"  ## compute all elements of a fourfold table (num)',
    '))',
    '',
    ''),
    if(encoding == "") NULL else c(
    '## encoding -------------------------------------------------------------------------',
    '',
    sprintf('exams2html("%s",', exrc_enc),
    sprintf('  encoding = "%s",', encoding),
    sprintf('  template = "%s")', if(absolute) file.path(dir, "templates", "plain.html") else file.path("..", "templates", "plain.html")),
    sprintf('exams2pdf("%s",', exrc_enc),
    sprintf('  encoding = "%s",', encoding),
    sprintf('  template = "%s")', if(absolute) file.path(dir, "templates", "plain.tex") else file.path("..", "templates", "plain.tex")),
    '',
    ''),
    '## other interfaces -----------------------------------------------------------------',
    '',
    '## switch back to the original folder',
    sprintf('setwd("%s")', if(absolute) dir else ".."),
    '',
    '## other interfaces include:',
    '## - exams2pdf for customizable PDF output',
    '## - exams2nops for a fixed PDF format that can be automatically scanned and evaluated',
    '##',
    '## - exams2html for customizable HTML output',
    '## - exams2moodle for Moodle XML that can be imported into Moodle quizzes',
    '## - exams2qti12/exams2qti21 for QTI XML (version 1.2 or 2.1) that can be imported',
    '##   into various learning management systems (e.g., OLAT or OpenOLAT among others)',
    '##',
    '## - exams2arsnova for a JSON format that can be imported into ARSnova live quizzes ',
    '',
    '## see the demo-*.R scripts in this directory for more examples',
    'dir()',
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )
  if(markup == "markdown") demo_all <- gsub(".Rnw", ".Rmd", demo_all, fixed = TRUE)

  demo_exams <- c(
    '## exams ----------------------------------------------------------------------------',
    '',
    '## load package',
    'library("exams")',
    '',
    sprintf('## exam with a simple vector of exercises in %s format',
      if(markup == "latex") "R/LaTeX (.Rnw)" else "R/Markdown (.Rmd)"),
    '## -> alternatively try a list of vectors of more exercises',
    sprintf('myexam <- c("%s")', paste(exrc, collapse = '", "')),
    if("currency" %in% substr(exrc, 1, 8)) sprintf(
    '## note that the currency exercise is in %s encoding',
      paste(unique(c(encoding, charset)), collapse = "/")) else NULL,
    '',
    ''
  )

  demo_html <- c(demo_exams,
    '## exams2html -----------------------------------------------------------------------',
    '## HTML output (1 file per exam)',
    '## -> typically used for quickly checking if an exercise can be converted to HTML',
    '## -> or customized via suitable templates',
    '',
    '## generate the HTML version of a single exercise (shown in browser)',
    '## with default settings',
    sprintf('exams2html("%s")', if(absolute) file.path(dir, "exercises", exrc[1L]) else file.path("exercises", exrc[1L])),
    '',
    '## generate a single HTML exam (shown in browser)',
    '## with specification of a template (the default) %s encoding',
    'exams2html(myexam, n = 1,',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
    sprintf('  template = "%s")', if(absolute) file.path(dir, "templates", "plain.html") else file.path("templates", "plain.html")),
    '',
    '## generate three HTML exams without solutions in output directory',
    'exams2html(myexam, n = 3, name = "html-demo", solution = FALSE,',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
    sprintf('  template = "%s")', if(absolute) file.path(dir, "templates", "plain.html") else file.path("templates", "plain.html")),
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )

  demo_pdf <- c(demo_exams,
    '## exams2pdf ------------------------------------------------------------------------',
    '## PDF output (1 file per exam)',
    '## -> typically used for quickly checking if an exercise can be converted to PDF',
    '## -> or customized via suitable templates',
    '',
    '## generate the PDF version of a single exercise (shown in PDF viewer)',
    '## with default settings',
    sprintf('exams2pdf("%s")', if(absolute) file.path(dir, "exercises", exrc[1L]) else file.path("exercises", exrc[1L])),
    '',
    '## generate a single PDF exam (shown in PDF viewer)',
    '## with specification of a template (for an exam) %s encoding',
    'exams2pdf(myexam, n = 1,',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
    sprintf('  template = "%s")', if(absolute) file.path(dir, "templates", "exam.tex") else file.path("templates", "solution.tex")),
    '',
    '## generate three PDF exams and corresponding solutions in output directory',
    '## (with the header used to set a custom Date and ID for the exam)',
    'exams2pdf(myexam, n = 3, name = c("pdf-exam", "pdf-solution"),',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
    sprintf('  template = c("%s", "%s"),',
      if(absolute) file.path(dir, "templates", "exam.tex") else file.path("templates", "exam.tex"),
      if(absolute) file.path(dir, "templates", "exam.tex") else file.path("templates", "solution.tex")),
    '  header = list(',
    '    Date = "2015-01-01",',
    '    ID = function(i) formatC(i, width = 5, flag = "0")',
    '  ))',
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )

  demo_moodle <- c(demo_exams,
    '## exams2moodle ---------------------------------------------------------------------',
    '## Moodle XML output (1 file containing all exams)',
    '## -> for import into a Moodle system',
    '',
    '## generate Moodle exam with three replications per question',
    'exams2moodle(myexam, n = 3, name = "moodle-demo",',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),    
    sprintf('  edir = "%s")', if(absolute) file.path(dir, "exercises") else "exercises"),
    '',
    '## hint: to quickly check (prior to Moodle export) whether each exercise can be',
    '## converted to HTML, exams2html() can be used',
    sprintf('exams2html("%s")', if(absolute) file.path(dir, "exercises", exrc[1L]) else file.path("exercises", exrc[1L])),
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )

  demo_qti <- c(demo_exams,
    if(!("exams2qti12" %in% writer)) NULL else c(
    '## exams2qti12 ----------------------------------------------------------------------',
    '## XML output in QTI 1.2 format (1 file containing all exams, zipped by default)',
    '## -> for import into QTI-based learning management systems (e.g., OLAT/OpenOLAT)',
    '',
    '## generate QTI 1.2 exam with three replications per question',
    '## (showing correct solutions after failed attempts and passing if solving',
    '## at least two items)',
    'exams2qti12(myexam, n = 3, name = "qti12-demo",',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),    
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
            '  solutionswitch = TRUE, maxattempts = 1, cutvalue = 2)',
    '',
    '## hint: if customization of the QTI 1.2 template is necessary, modify',
    '## the file templates/qti12.xml and then set the argument:',
    sprintf('## template = "%s"', if(absolute) file.path(dir, "templates", "qti12.xml") else file.path("templates", "qti12.xml")),
    '',
    '## hint: to quickly check (prior to QTI export) whether each exercise can be',
    '## converted to HTML, exams2html() can be used',
    sprintf('exams2html("%s")', if(absolute) file.path(dir, "exercises", exrc[1L]) else file.path("exercises", exrc[1L])),
    '',
    ''),
    if(!("exams2qti21" %in% writer)) NULL else c(
    '## exams2qti21 ----------------------------------------------------------------------',
    '## XML output in QTI 2.1 format (1 file containing all exams, zipped by default)',
    '## -> for import into QTI-based learning management systems',
    '',
    '## generate QTI 2.1 exam with three replications per question',
    '## (showing correct solutions after failed attempts and passing if solving',
    '## at least two items)',
    'exams2qti21(myexam, n = 3, name = "qti21-demo",',
    if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
    sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),    
    sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
            '  solutionswitch = TRUE, maxattempts = 1, cutvalue = 2)',
    '',
    '## hint: if customization of the QTI 2.1 template is necessary, modify',
    '## the file templates/qti21.xml and then set the argument:',
    sprintf('## template = "%s"', if(absolute) file.path(dir, "templates", "qti21.xml") else file.path("templates", "qti21.xml")),
    '',
    '## hint: to quickly check (prior to QTI export) whether each exercise can be',
    '## converted to HTML, exams2html() can be used',
    sprintf('exams2html("%s")', if(absolute) file.path(dir, "exercises", exrc[1L]) else file.path("exercises", exrc[1L])),
    '',
    ''),
    '## ----------------------------------------------------------------------------------'
  )

  ## omit cloze exercises for exams2arsnova
  exrc <- exrc[substr(exrc, 1, 7) != "boxhist"]
  demo_exams[8L] <- sprintf('myexam <- c("%s")', paste(exrc, collapse = '", "'))
  
  demo_arsnova <- c(demo_exams,
    '## exams2arsnova --------------------------------------------------------------------',
    '## JSON output in ARSnova format (1 file per exam)',
    '## -> for import into an ARSnova session',
    '',
    '## generate ARSnova exam/quiz with one replication per question',
    '## (without allowing abstentions)',
    'exams2arsnova(myexam, n = 1, name = "arsnova-demo",',
      if(encoding != "") sprintf('  encoding = "%s",', encoding) else NULL,
      sprintf('  dir = "%s",', if(absolute) file.path(dir, "output") else "output"),    
      sprintf('  edir = "%s",', if(absolute) file.path(dir, "exercises") else "exercises"),
    '  abstention = FALSE)',
    '',
    '## hint: if the arguments url/sessionkey/jsessionid are specified the JSON output',
    '## can be directly sent to the running ARSnova session (without manually saving and',
    '## importing it)',
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )  

  demo_nops <- c(
    '## exams ----------------------------------------------------------------------------',
    '',
    '## load package',
    'library("exams")',
    '',
    '## switch to dedicated directory for nops output',
    sprintf('setwd("%s")', if(absolute) file.path(dir, "nops") else "nops"),
    '',
    '',
    '## exams2nops -----------------------------------------------------------------------',
    '## PDF output in NOPS format (1 file per exam)',
    '## -> for exams that can be printed, scanned, and automatically evaluated in R',
    '',
    '',
    '## --- Step 1 ---',
    '## exam generation',
    '',
    '## define an exam (= list of exercises)',
    'myexam <- list(',
    '  "tstat2.Rnw",',
    '  "ttest.Rnw",',
    '  "relfreq.Rnw",',
    '  "anova.Rnw",',
    '  c("boxplots.Rnw", "scatterplot.Rnw"),',
    '  "cholesky.Rnw"',
    ')',
    '',
    '## create multiple exams on the disk with different numbers of points',
    '## per exercise (see ?exams2nops for more examples)',
    'set.seed(403)',
    'ex1 <- exams2nops(myexam, n = 2, dir = ".", date = "2015-07-29",',
    '  points = c(1, 1, 1, 2, 2, 3), showpoints = TRUE)',
    'dir()',
    '',
    '## assume the PDF exams were already printed (and possibly backed up',
    '## in a different directory) so that they are not needed anymore',
    'file.remove(dir(pattern = "pdf$"))',
    '',
    '',
    '## --- Step 2 ---',
    '## scan results',
    '',
    '## assume two participants filled out the printed exam sheets',
    '## and the corresponding scans are in two PNG files,',
    'img <- dir(system.file("nops", package = "exams"), pattern = "nops_scan",',
    '  full.names = TRUE)',
    '',
    '## copy the PNG files to the working directory',
    'file.copy(img, to = ".")',
    '',
    '## read the scanned images (all locally available .png files) and collect',
    '## results in a ZIP archive (see ?nops_scan for more details)',
    'nops_scan()',
    'dir()',
    '',
    '## the ZIP archive contains copies of the PNG images so that these are',
    '## can be deleted here (possibly after backup in a different directory)',
    'file.remove(dir(pattern = "png$"))',
    '',
    '',
    '## -- Step 3 ---',
    '## evaluate results',
    '',
    '## three files are required: (a) an RDS file with the exam meta-information',
    '## (see Step 1), (b) a ZIP file with the scanned sheets (see Step 2), (c) a',
    '## CSV file with the student infomation (registration number, name, and some',
    '## for of ID/username)',
    '',
    '## here we create the CSV file on the fly but in practice this will typically',
    '## be processed from some registration service or learning management system etc',
    'write.table(data.frame(',
    '  registration = c("1501090", "9901071"),',
    '  name = c("Jane Doe", "Ambi Dexter"),',
    '  id = c("jane_doe", "ambi_dexter")',
    '), file = "Exam-2015-07-29.csv", sep = ";", quote = FALSE, row.names = FALSE)',
    'dir()',
    '## now the exam can be evaluated creating an output data frame (also stored',
    '## as CSV file) and individual HTML reports (stored in a ZIP file),',
    '',
    '## as there is only exactly on CSV/RDS/ZIP file in the current directory,',
    '## these are found automatically - furthermore an evaluation scheme without',
    '## partial points and differing points per exercise are used',
    'ev1 <- nops_eval(eval = exams_eval(partial = FALSE, negative = FALSE))',
    'dir()',
    '',
    '## inspect evaluated data',
    'ev1',
    '',
    '## inspect corresponding HTML reports',
    'unzip("nops_eval.zip")',
    'browseURL(file.path("jane_doe",    "Exam-2015-07-29.html"))',
    'browseURL(file.path("ambi_dexter", "Exam-2015-07-29.html"))',
    '',
    '',
    '## --- Options ---',
    '## below three typically needed options are discussed:',
    '## (a) using a different evaluation strategy (here with partial credits),',
    '## (b) using a different language (here de/German),',
    '## (c) an error of the participant when filling in the registration number.',
    '',
    '## as for (a): partial credits should only be used for multiple-choice questions',
    '## where at least one alternative is correct and at least one is false',
    '## [note that in this example this is not the case for the first question',
    '## (single-choice) and the third question for Jane Doe (no alternative correct)]',
    '',
    '## as for (c): for Ambi Dexter such an error was included in the PNG example',
    '## image, the actual number is "9911071" but the crosses indicate "9901071"',
    '',
    '## clean up previous evaluation',
    'file.remove(c("nops_eval.csv", "nops_eval.zip"))',
    '',
    '## write correct registration information',
    'write.table(data.frame(',
    '  registration = c("1501090", "9911071"),',
    '  name = c("Jane Doe", "Ambi Dexter"),',
    '  id = c("jane_doe", "ambi_dexter")',
    '), file = "Exam-2015-07-29.csv", sep = ";", quote = FALSE, row.names = FALSE)',
    '',
    '## call nops_eval() with modified options, where the error in the registration',
    '## number of Ambi Dexter will trigger an interactive prompt',
    'ev2 <- nops_eval(eval = exams_eval(partial = TRUE, rule = "false2"),',
    '  language = "de")',
    '',
    '## inspect evaluated data',
    'ev2',
    'cbind(ev1$points, ev2$points)',
    '',
    '## inspect corresponding HTML reports',
    'unzip("nops_eval.zip")',
    'browseURL(file.path("jane_doe",    "Exam-2015-07-29.html"))',
    'browseURL(file.path("ambi_dexter", "Exam-2015-07-29.html"))',
    '',
    '',
    '## ----------------------------------------------------------------------------------'
  )
  if(markup == "markdown") demo_nops <- gsub(".Rnw", ".Rmd", demo_nops, fixed = TRUE)

  writeLines(demo_all, file.path(dir, "demo-all.R"))
  if("exams2html" %in% writer)    writeLines(demo_html,    file.path(dir, "demo-html.R"))
  if("exams2pdf" %in% writer)     writeLines(demo_pdf,     file.path(dir, "demo-pdf.R"))
  if("exams2moodle" %in% writer)  writeLines(demo_moodle,  file.path(dir, "demo-moodle.R"))
  if(any(c("exams2qti12", "exams2qti21") %in% writer)) writeLines(demo_qti, file.path(dir, "demo-qti.R"))
  if("exams2arsnova" %in% writer) writeLines(demo_arsnova, file.path(dir, "demo-arsnova.R"))
  if("exams2nops" %in% writer)    writeLines(demo_nops,    file.path(dir, "demo-nops.R"))

  ## return all scripts invisibly
  invisible(list(
    all = demo_all,
    html = demo_html,
    pdf = demo_pdf,
    moodle = demo_moodle,
    qti = demo_qti,
    arsnova = demo_arsnova,
    nops = demo_nops
  ))
}
