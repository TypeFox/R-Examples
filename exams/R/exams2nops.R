exams2nops <- function(file, n = 1L, dir = NULL, name = NULL,
  language = "en", title = "Exam", course = "",
  institution = "R University", logo = "Rlogo.png", date = Sys.Date(), 
  replacement = FALSE, intro = NULL, blank = NULL, duplex = TRUE, pages = NULL,
  usepackage = NULL, header = NULL, encoding = "", startid = 1L, points = NULL, showpoints = FALSE, ...)
{
  ## pages could include formulary and distribution tables
  if(!is.null(pages)) pages <- sapply(pages, file_path_as_absolute)

  ## header: date, id, usepackage
  Date2ID <- function(date) function(i) paste(format(date, "%y%m%d"),
    formatC(i + startid - 1L, width = 5, flag = 0, format = "f", digits = 0), sep = "")
  if(!inherits(date, "Date")) date <- as.Date(date)
  d2id <- Date2ID(date)
  if(!is.null(usepackage)) {
    usepackage <- as.list(usepackage)
    names(usepackage) <- rep.int("usepackage", length(usepackage))
  }
  ## header: localization (titles, logos, etc.)
  if(missing(logo)) logo <- system.file(file.path("nops", "Rlogo.png"), package = "exams")
  if(course != "") course <- paste0("(", course, ")")
  loc <- list(
    nopsinstitution = institution,
    nopstitle = title,
    nopscourse = course,
    "newcommand{\\mylogo}" = logo
  )
  ## header: internationalization
  if(!file.exists(language)) language <- system.file(file.path("nops", paste0(language, ".dcf")), package = "exams")
  if(language == "") language <- system.file(file.path("nops", "en.dcf"), package = "exams")
  lang <- nops_language(language, markup = "latex")[c(
    "PersonalData", "FamilyName", "GivenName", "Signature", "RegistrationNumber", 
    "Checked", "NoChanges", "DocumentType", "DocumentID", "Scrambling", 
    "Replacement", "MarkCarefully", "NotMarked", "Or",
    "MarkExampleA", "MarkExampleB", "MarkExampleC", "MarkExampleD", "MarkExampleE",
    "Warning", "Answers", "FillAnswers", "Point", "Points")]
  ## header: collect everything
  header <- c(list(Date = date, ID = d2id), usepackage, loc, lang, header)

  ## determine number of alternative choice fors each exercise
  ufile <- unique(unlist(file))
  x <- exams_metainfo(xexams(ufile, driver = list(sweave = list(quiet = TRUE, encoding = encoding),
    read = NULL, transform = NULL, write = NULL), ...))[[1L]]
  names(x) <- ufile
  utype <- sapply(ufile, function(n) x[[n]]$type)
  wrong_type <- ufile[utype == "cloze"]
  if(length(wrong_type) > 0L) {
    stop(paste("the following exercises are cloze exercises:",
      paste(wrong_type, collapse = ", ")))
  }
  x <- sapply(ufile, function(n) x[[n]]$length)
  x[!(utype %in% c("schoice", "mchoice"))] <- 0L
  if(any(x == 1L | x > 5L)) {
    stop(paste("the following exercises have length < 2 or > 5:",
      paste(names(x)[x == 1L | x > 5], collapse = ", ")))
  }
  if(sum(x < 1L) > 3L) {
    stop(paste("currently only up to three exercises that are not schoice/mchoice are supported:",
      paste(names(x)[x < 1L], collapse = ", ")))
  }
  if(is.list(file)) {
    nchoice <- lapply(file, function(n) x[n])
    nchoice1 <- as.vector(sapply(nchoice, min))
    nchoice <- as.vector(sapply(nchoice, max))
    if(any(nchoice != nchoice1)) {
      stop(paste("the following groups of exercise do not have the same length:",
        paste(sapply(file, paste, collapse = "/")[nchoice != nchoice1], collapse = ", ")))
    }
  } else {
    nchoice <- as.vector(x[file])
  }

  ## generate appropriate template on the fly
  template <- file.path(tempdir(), "nops.tex")
  make_nops_template(length(file), replacement = replacement, intro = intro,
    blank = blank, duplex = duplex, pages = pages,
    file = template, nchoice = nchoice, encoding = encoding)

  ## if points should be shown generate a custom transformer
  transform <- if(showpoints) {
    to_latex <- make_exercise_transform_pandoc(to = "latex", base64 = FALSE)
    function(x) {
      x <- to_latex(x)
      x$question <- c(
        sprintf("\\emph{(%s %s)}", x$metainfo$points, ifelse(x$metainfo$points != 1, "\\myPoints", "\\myPoint")),
	x$question
      )
      return(x)
    }
  } else {
    NULL
  }

  if(is.null(dir)) {  
    rval <- exams2pdf(file, n = n, name = name, template = template,
      header = header, transform = transform, encoding = encoding,
      points = points, ...)
    names(rval) <- d2id(1:length(rval))
  } else {
    rval <- exams2pdf(file, n = n, dir = dir, name = name, template = template,
      header = header, transform = transform, encoding = encoding,
      points = points, ...)
    names(rval) <- d2id(1:length(rval))
    if(is.null(name)) name <- "metainfo"
    name <- paste(name, ".rds", sep = "")
    saveRDS(rval, file = file.path(dir, name))
  }
  
  invisible(rval)
}

make_nops_template <- function(n, replacement = FALSE, intro = NULL, blank = NULL,
  duplex = TRUE, pages = NULL, file = NULL, nchoice = 5, encoding = "")
{
page1 <- make_nops_page(n, nchoice = nchoice)
page2 <- if(replacement) {
  make_nops_page(n, nchoice = nchoice, replacement = TRUE)
} else {
  NULL
}
page3 <- if(any(nchoice < 1L)) {
  make_nops_string_page(which(nchoice < 1L), nchoice = 5) ## FIXME: possibly other nchoice values?
} else {
  NULL
}

## encoding
enc <- gsub("-", "", tolower(encoding), fixed = TRUE)
if(enc %in% c("iso8859", "iso88591")) enc <- "latin1"
if(enc == "iso885915") enc <- "latin9"

empty <- if(!duplex) {
"
\\newpage
"
} else {
"
\\newpage
\\thispagestyle{empty}
\\phantom{.}
"
}

if(is.null(blank)) blank <- ceiling(n/2)
blank <- rep("\\newpage\n\\phantom{.}", blank)

rval <- c("
\\documentclass[10pt,a4paper]{article} 
\\usepackage{graphicx,color}
\\usepackage{amsmath,amssymb,latexsym}
\\usepackage{verbatim,url,fancyvrb,ae}
\\usepackage{multicol,a4wide,pdfpages}
\\IfFileExists{sfmath.sty}{
  \\RequirePackage{sfmath}
}{}

\\DefineVerbatimEnvironment{Sinput}{Verbatim}{fontshape=sl}
\\DefineVerbatimEnvironment{Soutput}{Verbatim}{}
\\DefineVerbatimEnvironment{Scode}{Verbatim}{fontshape=sl}
\\newenvironment{Schunk}{}{}

\\usepackage[T1]{fontenc}",
if(enc != "") sprintf('\\usepackage[%s]{inputenc}', enc) else NULL,
"
\\renewcommand{\\rmdefault}{phv}
\\renewcommand{\\sfdefault}{phv}

\\setlength{\\parskip}{0.7ex plus0.1ex minus0.1ex}
\\setlength{\\parindent}{0em}
\\setlength{\\textheight}{29.6cm} 
\\setlength{\\oddsidemargin}{-2.54cm} 
\\setlength{\\evensidemargin}{-2.54cm} 
\\setlength{\\topmargin}{-2.54cm} 
\\setlength{\\headheight}{0cm} 
\\setlength{\\headsep}{0cm} 
\\setlength{\\footskip}{0cm} 
\\setlength{\\unitlength}{1mm} 
\\usepackage{chngpage}

%% for exams2pdf
\\newenvironment{question}{\\item}{}
\\newenvironment{solution}{\\comment}{\\endcomment}
\\newenvironment{answerlist}{\\renewcommand{\\labelenumi}{(\\alph{enumi})}\\begin{enumerate}}{\\end{enumerate}}

%% additional header commands
\\makeatletter
\\newcommand{\\ID}[1]{\\def\\@ID{#1}}
\\newcommand{\\Date}[1]{\\def\\@Date{#1}}
%
\\newcommand{\\nopsinstitution}[1]{\\def\\@nopsinstitution{#1}}
\\newcommand{\\nopstitle}[1]{\\def\\@nopstitle{#1}}
\\newcommand{\\nopscourse}[1]{\\def\\@nopscourse{#1}}
%
\\newcommand{\\PersonalData}[1]{\\def\\@PersonalData{#1}}
\\newcommand{\\FamilyName}[1]{\\def\\@FamilyName{#1}}
\\newcommand{\\GivenName}[1]{\\def\\@GivenName{#1}}
\\newcommand{\\Signature}[1]{\\def\\@Signature{#1}}
\\newcommand{\\RegistrationNumber}[1]{\\def\\@RegistrationNumber{#1}}
\\newcommand{\\Checked}[1]{\\def\\@Checked{#1}}
\\newcommand{\\NoChanges}[1]{\\def\\@NoChanges{#1}}
\\newcommand{\\DocumentType}[1]{\\def\\@DocumentType{#1}}
\\newcommand{\\DocumentID}[1]{\\def\\@DocumentID{#1}}
\\newcommand{\\Scrambling}[1]{\\def\\@Scrambling{#1}}
\\newcommand{\\Replacement}[1]{\\def\\@Replacement{#1}}
\\newcommand{\\MarkCarefully}[1]{\\def\\@MarkCarefully{#1}}
\\newcommand{\\NotMarked}[1]{\\def\\@NotMarked{#1}}
\\newcommand{\\Or}[1]{\\def\\@Or{#1}}
\\newcommand{\\MarkExampleA}[1]{\\def\\@MarkExampleA{#1}}
\\newcommand{\\MarkExampleB}[1]{\\def\\@MarkExampleB{#1}}
\\newcommand{\\MarkExampleC}[1]{\\def\\@MarkExampleC{#1}}
\\newcommand{\\MarkExampleD}[1]{\\def\\@MarkExampleD{#1}}
\\newcommand{\\MarkExampleE}[1]{\\def\\@MarkExampleE{#1}}
\\newcommand{\\Warning}[1]{\\def\\@Warning{#1}}
\\newcommand{\\Answers}[1]{\\def\\@Answers{#1}}
\\newcommand{\\FillAnswers}[1]{\\def\\@FillAnswers{#1}}
\\newcommand{\\Point}[1]{\\def\\@Point{#1}}
\\newcommand{\\Points}[1]{\\def\\@Points{#1}}

\\ID{YYMMDD00001}
\\Date{YYYY-MM-DD}
%
\\nopsinstitution{R University}
\\nopstitle{Exam}
\\nopscourse{}
%
\\PersonalData{Personal Data}
\\FamilyName{Family Name}
\\GivenName{Given Name}
\\Signature{Signature}
\\RegistrationNumber{Registration Number}
\\Checked{checked}
\\NoChanges{In this section no modifications of the data must be made!}
\\DocumentType{Type}
\\DocumentID{Exam ID}
\\Scrambling{Scrambling}
\\Replacement{Replacement}
\\MarkCarefully{Please mark the boxes carefully with a cross}
\\NotMarked{Not marked}
\\Or{or}
\\MarkExampleA{72}
\\MarkExampleB{80}
\\MarkExampleC{102}
\\MarkExampleD{109}
\\MarkExampleE{115}
\\Warning{This document is scanned automatically. Please keep clean and do not bend or fold. For filling in the document please use a \\textbf{blue or black pen}. \\\\ \\textbf{Only clearly marked and positionally accurate crosses will be processed!}}
\\Answers{Answers}
\\FillAnswers{In the following please fill in your answers.}
\\Point{point}
\\Points{points}

%% \\exinput{header}

\\newcommand{\\myID}{\\@ID}
\\newcommand{\\myDate}{\\@Date}
%
\\newcommand{\\myinstitution}{\\@nopsinstitution}
\\newcommand{\\mytitle}{\\@nopstitle}
\\newcommand{\\mycourse}{\\@nopscourse}
%
\\newcommand{\\myPersonalData}{\\@PersonalData}
\\newcommand{\\myFamilyName}{\\@FamilyName}
\\newcommand{\\myGivenName}{\\@GivenName}
\\newcommand{\\mySignature}{\\@Signature}
\\newcommand{\\myRegistrationNumber}{\\@RegistrationNumber}
\\newcommand{\\myChecked}{\\@Checked}
\\newcommand{\\myNoChanges}{\\@NoChanges}
\\newcommand{\\myDocumentType}{\\@DocumentType}
\\newcommand{\\myDocumentID}{\\@DocumentID}
\\newcommand{\\myScrambling}{\\@Scrambling}
\\newcommand{\\myReplacement}{\\@Replacement}
\\newcommand{\\myMarkCarefully}{\\@MarkCarefully}
\\newcommand{\\myNotMarked}{\\@NotMarked}
\\newcommand{\\myOr}{\\@Or}
\\newcommand{\\myMarkExampleA}{\\@MarkExampleA}
\\newcommand{\\myMarkExampleB}{\\@MarkExampleB}
\\newcommand{\\myMarkExampleC}{\\@MarkExampleC}
\\newcommand{\\myMarkExampleD}{\\@MarkExampleD}
\\newcommand{\\myMarkExampleE}{\\@MarkExampleE}
\\newcommand{\\myWarning}{\\@Warning}
\\newcommand{\\myAnswers}{\\@Answers}
\\newcommand{\\myFillAnswers}{\\@FillAnswers}
\\newcommand{\\myPoint}{\\@Point}
\\newcommand{\\myPoints}{\\@Points}

\\makeatother

\\markboth{\\textsf{{\\mytitle}: {\\myID}}}{\\textsf{{\\mytitle}: {\\myID}}}
\\pagestyle{myheadings}
\\begin{document} 
",
page1,
empty,
if(replacement) {
  c("\n\\newpage\n", page2, empty)
},
if(length(page3)) {
  c("\n\\newpage\n", page3, empty)
},
"
\\setcounter{page}{0}

\\setlength{\\textheight}{24cm} 
\\setlength{\\oddsidemargin}{0cm} 
\\setlength{\\evensidemargin}{0cm} 
\\setlength{\\topmargin}{0cm} 
\\setlength{\\headheight}{0cm} 
\\setlength{\\headsep}{1cm} 
\\setlength{\\footskip}{1cm} 

\\newpage
",

intro,

"
\\begin{enumerate}

%% \\exinput{exercises}

\\end{enumerate}
",
if(!is.null(pages)) paste("\\newpage\n\\includepdf[pages=1-]{", pages, "}", sep = "", collapse = "\n"),
"",
blank,
"
\\end{document}
")

if(!is.null(file)) writeLines(rval, file)
invisible(rval)
}

make_nops_page <- function(n, replacement = FALSE, nchoice = 5) {

## number of alternative choices
nchoice <- rep(nchoice, length.out = n)

## helper function for abcde labels
abcde <- function(i, above = FALSE, nchoice = 5) {
  ix <- (i - 1) %/% 15
  iy <- (i - 1) %% 15 + 1
  ix <- 19 + 64 * ix - as.numeric(ix >= 2) * 4
  iy <- 129 - 7 * iy - 3 * ((iy - 1) %/% 5) + above * 10
  nchoice <- max(nchoice)
  if(nchoice == 5) {
    sprintf(paste("\\put(%i,%i){\\makebox(0,0)[b]{\\textsf{", letters[1:5],"}}}", sep = "", collapse = "\n"),
      ix + 1 * 8, iy, ix + 2 * 8, iy, ix + 3 * 8, iy, ix + 4 * 8, iy, ix + 5 * 8, iy)  
  } else if(nchoice == 4) {
    sprintf(paste("\\put(%i,%i){\\makebox(0,0)[b]{\\textsf{", letters[1:4],"}}}", sep = "", collapse = "\n"),
      ix + 1 * 8, iy, ix + 2 * 8, iy, ix + 3 * 8, iy, ix + 4 * 8, iy)
  } else if(nchoice == 3) {
    sprintf(paste("\\put(%i,%i){\\makebox(0,0)[b]{\\textsf{", letters[1:3],"}}}", sep = "", collapse = "\n"),
      ix + 1 * 8, iy, ix + 2 * 8, iy, ix + 3 * 8, iy)
  } else if(nchoice == 2) {
    sprintf(paste("\\put(%i,%i){\\makebox(0,0)[b]{\\textsf{", letters[1:2],"}}}", sep = "", collapse = "\n"),
      ix + 1 * 8, iy, ix + 2 * 8, iy)
  } else {
    stop("'nchoice' must be one of 5, 4, 3, 2")
  }
}

qbox <- function(i, nchoice = 5) {
  ix <- (i - 1) %/% 15
  iy <- (i - 1) %% 15 + 1
  ix <- 19 + 64 * ix - as.numeric(ix >= 2) * 4
  iy <- 129 - 7 * iy - 3 * ((iy - 1) %/% 5)
  
  if(nchoice > 0) {
    sprintf("\\put(%i,%i){\\makebox(0,0){\\textsf{%i}}}\n\\multiput(%i,%i)(8,0){%i}{\\framebox(4,4){}}",
      ix + 2, iy + 6, i, ix + 6, iy + 4, nchoice)
  } else {
    sprintf("\\put(%i,%i){\\makebox(0,0){\\textsf{%i}}}", ix + 2, iy + 6, i)
  }
}

c("
\\thispagestyle{empty}
\\begin{picture}(210,290) 
\\thicklines 

% position marks for scanning
\\put(17.5,13){\\line(1,0){5}} \\put(20,10.5){\\line(0,1){5}} 
\\put(187.5,13){\\line(1,0){5}} \\put(190,10.5){\\line(0,1){5}} 
\\put(157.5,270){\\line(1,0){5}} \\put(160,267.5){\\line(0,1){5}} 
\\put(27.5,270){\\line(1,0){5}} \\put(30,267.5){\\line(0,1){5}} 

% personal data box
\\put(72.5,244){\\makebox(0,0){\\textsf{\\myPersonalData}}} 
\\put(20,198){\\framebox(105,43){}} \\thinlines 
\\multiput(20,217)(0,12){2}{\\line(1,0){90}} \\thicklines 
\\put(21,236){\\makebox(0,5)[l]{\\textsf{\\myFamilyName:}}} 
\\put(21,224){\\makebox(0,5)[l]{\\textsf{\\myGivenName:}}} 
\\put(21,212){\\makebox(0,5)[l]{\\textsf{\\mySignature:}}} 
\\put(123,200){\\makebox(0,0)[rb]{\\footnotesize{\\textsf{\\myChecked}}}} 

% registration number box
\\put(159,244){\\makebox(0,0){\\textsf{\\myRegistrationNumber}}} 
\\put(131,233){\\framebox(56,8){}} \\thinlines 
\\multiput(139,233)(8,0){6}{\\line(0,1){1.5}} \\thicklines 
\\multiput(133,163)(8,0){7}{\\begin{picture}(0,0) 
\\multiput(0,0)(0,7){10}{\\framebox(4,4){}}\\end{picture}}",
if(replacement) "\\setcounter{nr3}{0}" else "\\newcounter{nr3}",
"
\\multiput(129,228)(0,-7){10}{\\begin{picture}(0,0) 
\\multiput(0,0)(60,0){2}{\\makebox(0,0){\\textsf{\\arabic{nr3}}}}
\\end{picture} \\stepcounter{nr3}} 
% general instructions and logo
\\IfFileExists{\\mylogo}{\\put(175,251){\\includegraphics[height=2.51cm,keepaspectratio]{\\mylogo}}}{}
\\put(40,270){\\makebox(0,0)[bl]{\\textsf{\\textbf{\\LARGE{\\myinstitution}}}}}
\\put(20,147){\\parbox{170mm}{\\textsf{\\myWarning}}} 

% mark examples
\\put(20,158){\\makebox(0,0)[l]{\\textsf{\\myMarkCarefully:}}}
\\put(\\myMarkExampleB,158){\\makebox(0,0)[l]{\\textsf{\\myNotMarked:}}}
\\put(\\myMarkExampleD,158){\\makebox(0,0)[l]{\\textsf{\\myOr}}}
\\put(\\myMarkExampleA,157){\\framebox(4,4){}} 
\\put(\\myMarkExampleA,157){\\line(1,1){4}} \\put(\\myMarkExampleA,161){\\line(1,-1){4}} 
\\put(\\myMarkExampleA.2,157){\\line(1,1){3.8}} \\put(\\myMarkExampleA.2,161){\\line(1,-1){3.8}} 
\\put(\\myMarkExampleA,157.2){\\line(1,1){3.8}} \\put(\\myMarkExampleA,160.8){\\line(1,-1){3.8}} 
\\put(\\myMarkExampleC,157){\\framebox(4,4){}} 
\\put(\\myMarkExampleE,158){\\colorbox{black}{\\framebox(2,2){}}} 


% title and date
\\put(40,262){\\parbox[t]{120mm}{\\large{\\textsf{\\textbf{{\\mytitle} {\\myDate}}}}}}

% boxes for answers (inlcuding labels and separators)
",
## first column
## title
sprintf("\\put(43,138){\\makebox(0,0){\\textsf{{\\myAnswers} 1 - %s}}}", min(15, n)),
## labels
abcde(1, above = TRUE, nchoice = nchoice[1:min(15, n)]),
abcde(min(15, n), above = FALSE, nchoice = nchoice[1:min(15, n)]),
"",
## if second column
if(n > 15) {
  c(
  ## separator
  "\\put(72,18){\\line(0,1){121}}",
  ## title
  sprintf("\\put(107,138){\\makebox(0,0){\\textsf{{\\myAnswers} 16 - %i}}}", min(30, n)),
  ## labels
  abcde(16, above = TRUE, nchoice = nchoice[16:min(30, n)]),
  abcde(min(30, n), above = FALSE, nchoice = nchoice[16:min(30, n)]))
},
## if third column
if(n > 30) {
  c(
  ## separator
  "\\put(134,18){\\line(0,1){121}}",
  ## title
  sprintf("\\put(167,138){\\makebox(0,0){\\textsf{{\\myAnswers} 31 - %i}}}", n),
  ## labels
  abcde(31, above = TRUE, nchoice = nchoice[31:n]),
  abcde(n, above = FALSE, nchoice = nchoice[31:n]))
},
## box for each question
sapply(1:n, function(i) qbox(i, nchoice = nchoice[i])),
"
% block with id, scrambling, type, replacement box
\\linethickness{0.5mm} \\put(20,164){\\framebox(105,28){}} \\thicklines  
\\put(32,177){\\makebox(0,0)[t]{\\textsf{\\myDocumentType}}} 
\\put(25,166){\\framebox(14,7){}} 
\\put(67,177){\\makebox(0,0)[t]{\\textsf{\\myDocumentID \\mycourse}}}
\\put(46,166){\\framebox(42,7){}} \\put(25,183.5){\\parbox{70mm}{%
\\textsf{\\myNoChanges}}}
\\thinlines \\put(113,180){\\line(0,1){1.5}} \\thicklines 
\\put(113,191){\\makebox(0,0)[t]{\\textsf{\\textbf{\\myScrambling}}}} 
\\put(106,180){\\framebox(14,7){}}
\\put(67,169.5){\\makebox(0,0){\\Large{\\textsf{\\myID}}}}
% scrambling is currently always zero
\\put(109.5,183.5){\\makebox(0,0){\\Large{\\textsf{0}}}}
\\put(116.5,183.5){\\makebox(0,0){\\Large{\\textsf{0}}}}",

## type: the number of questions rounded up in steps of 5 
## (needed for uibk scanning services)
sprintf("\\put(32,169.5){\\makebox(0,0){\\Large{\\textsf{%s}}}}",
  formatC(5 * ((n - 1) %/% 5 + 1), width = 3, flag = "0")),

## replacement?
if(replacement) {
"
% replacement sheet
\\put(116,170){\\framebox(4,4){}}
\\put(114,172){\\makebox(0,0)[r]{\\textsf{\\myReplacement:}}}

% cross in replacement box
\\put(116,170){\\line(1,1){4}} \\put(116.1,174.15){\\line(1,-1){4}} 
\\put(116.2,169.9){\\line(1,1){3.8}} \\put(116.2,174){\\line(1,-1){3.8}} 
\\put(116,170.2){\\line(1,1){3.8}} \\put(116,173.8){\\line(1,-1){3.8}} 
"
},

"
\\end{picture}
")
}

make_nops_string_page <- function(labels, nchoice = 5) {

## number of alternative choices
n <- length(labels)
nchoice <- rep(nchoice, length.out = n)

c("
\\thispagestyle{empty}
\\begin{picture}(210,290) 
\\thicklines 

% position marks for scanning
\\put(17.5,13){\\line(1,0){5}} \\put(20,10.5){\\line(0,1){5}} 
\\put(187.5,13){\\line(1,0){5}} \\put(190,10.5){\\line(0,1){5}} 
\\put(157.5,270){\\line(1,0){5}} \\put(160,267.5){\\line(0,1){5}} 
\\put(27.5,270){\\line(1,0){5}} \\put(30,267.5){\\line(0,1){5}} 

% general instructions and logo
\\IfFileExists{\\mylogo}{\\put(175,251){\\includegraphics{\\mylogo}} }{}
\\put(40,270){\\makebox(0,0)[bl]{\\textsf{\\textbf{\\LARGE{\\myinstitution}}}}}
\\put(20,210){\\parbox{170mm}{\\textsf{\\myFillAnswers}}}",

switch(n,
"1" = sprintf("
\\put(23,204){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20, 22){\\framebox(170,183){}}", labels[1L]),
"2" = sprintf("
\\put(23,204){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20,116){\\framebox(170,89){}}
\\put(23,110){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20, 22){\\framebox(170,89){}}", labels[1L], labels[2L]),
"3" = sprintf("
\\put(23,204){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20,148){\\framebox(170,57){}}
\\put(23,141){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20, 85){\\framebox(170,57){}}
\\put(23, 78){\\makebox(0,0)[t]{\\textsf{%s}}} 
\\put(20, 22){\\framebox(170,57){}}", labels[1L], labels[2L], labels[3L])
),

"
% title and date
\\put(40,262){\\parbox[t]{120mm}{\\large{\\textsf{\\textbf{{\\mytitle} {\\myDate}}}}}}

% box for text answers
",

sprintf("\\put(106,237){\\makebox(0,0){\\textsf{%s}}}", labels[1L]),
sprintf("\\multiput(110,235)(8,0){%s}{\\framebox(4,4){}}", nchoice[1L]),

if(n > 1L) { c(
sprintf("\\put(106,230){\\makebox(0,0){\\textsf{%s}}}", labels[2L]),
sprintf("\\multiput(110,228)(8,0){%s}{\\framebox(4,4){}}", nchoice[2L])
)} else NULL,

if(n > 2L) { c(
sprintf("\\put(106,223){\\makebox(0,0){\\textsf{%s}}}", labels[3L]),
sprintf("\\multiput(110,221)(8,0){%s}{\\framebox(4,4){}}", nchoice[3L])
)} else NULL,
 
"
% block with id, scrambling, type
\\linethickness{0.5mm} \\put(20,217){\\framebox(140,28){}} \\thicklines  
\\put(32,230){\\makebox(0,0)[t]{\\textsf{\\myDocumentType}}} 
\\put(25,219){\\framebox(14,7){}} 
\\put(67,230){\\makebox(0,0)[t]{\\textsf{\\myDocumentID}}} 
\\put(46,219){\\framebox(42,7){}}
\\put(25,236.5){\\parbox{70mm}{%
\\textsf{\\myNoChanges}}}
\\put(32,222.5){\\makebox(0,0){\\Large{\\textsf{999}}}}
\\put(67,222.5){\\makebox(0,0){\\Large{\\textsf{\\myID}}}}

\\end{picture} 
")
}

nops_language <- function(file, markup = c("latex", "html"))
{
  ## read file
  lang <- drop(read.dcf(file))
  
  ## necessary fields for a correct lanuage specification
  langfields <- c("PersonalData", "FamilyName", "GivenName", "Signature", "RegistrationNumber", 
    "Checked", "NoChanges", "DocumentType", "DocumentID", "Scrambling", 
    "Replacement", "MarkCarefully", "NotMarked", "Or",
    "MarkExampleA", "MarkExampleB", "MarkExampleC", "MarkExampleD", "MarkExampleE",
    "Warning", "Answers", "FillAnswers", "Point", "Points",
    "ExamResults", "Evaluation", "Mark", "Question", "GivenAnswer", "CorrectAnswer",
    "ExamSheet")
  if(!all(langfields %in% names(lang))) stop("invalid language specification")
  
  ## desired output markup
  markup <- match.arg(tolower(markup), c("latex", "html"))
  if(markup == "html") lang <- structure(tth::tth(lang), .Names = names(lang))
  
  return(as.list(lang))
}
