## generate exams in .xml format
exams2lops <- function(file, n = 1L, nsamp = NULL, dir = ".",
  name = NULL, quiet = TRUE, edir = NULL, tdir = NULL, sdir = NULL, verbose = FALSE,
  solution = TRUE, doctype = NULL, head = NULL, resolution = 100,
  width = 4, height = 4, encoding = "", converter = "tex2image", base64 = FALSE,
  auto_scramble = TRUE, ...)
{
  ## set up .xml transformer and writer function
  htmltransform <- make_exercise_transform_html(converter = converter, base64 = base64, ...)
  lopswrite <- make_exams_write_lops(name, auto_scramble = auto_scramble, ...)

  ## create final .xml exam
  xexams(file, n = n, nsamp = nsamp, driver = list(
      sweave = list(quiet = quiet, pdf = FALSE, png = TRUE, resolution = resolution, width = width, height = height, encoding=encoding),
      read = NULL,
      transform = htmltransform,
      write = lopswrite),
    dir = dir, edir = edir, tdir = tdir, sdir = sdir, verbose = verbose)
}


## writes the final .xml site
make_exams_write_lops <- function(name = NULL, auto_scramble = TRUE, ...)
{
  function(x, dir, info)
  {
    args <- list(...)

    if(is.null(name))
      name <- "lopsexam"
    name <- paste(name, as.character(Sys.Date()), info$id, sep = "-")

    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
    dir.create(test_dir <- file.path(tdir, name))

    ## start to write the exam.xml
    exam.xml <- c(
      '<?xml version=\'1.0\' encoding=\'utf-8\'?>',
      '<learning_resources>',
      paste('<exam area="sd" restype="exam" shortname="', if(is.null(args$shortname)) name else args$shortname, '">', sep = ""),
      '<metadata>',
      paste('<author>', if(is.null(args$author)) 'Mr. Exam' else args$author, '</author>', sep = ""),
      paste('<title>', if(is.null(args$title)) name else args$title, '</title>', sep = ""),
      paste('<language>', if(is.null(args$language)) 'en' else args$language, '</language>', sep = ""),
      '<category area="sd" restype="csp" shortname="sd_root_concept"/>',
      '</metadata>',
      '<description/>',
      paste('<grading area="', if(is.null(args$gradingarea)) 'sd' else args$gradingarea,
        '" shortname="', if(is.null(args$gradingshortname)) 'teilgenommen' else args$gradingshortname, '"/>', sep = ""),
      paste('<usage_details allowed_time="0" final_exam="false" practice="false" auto_scramble="', tolower(as.character(auto_scramble)), '"/>', sep = ""),
      '<directory>'
    )

    ## cycle trough th questions
    for(ex in x) {
      type <- ex$metainfo$type
      class(ex) <- c(if(type == "schoice") "mchoice" else type, "list")
      exam.xml <- c(exam.xml, write_lops_question(ex, test_dir, ...))
    }

    ## finish the exam xml
    exam.xml <- c(
      exam.xml,
      '</directory>',
      '</exam>',
      '</learning_resources>'
    )

    ## create exam directory and write the xml file into it
    dir.create(exam_dir <- file.path(test_dir, paste("sd_exam", name, sep = "_")))
    writeLines(exam.xml, file.path(exam_dir, paste(basename(exam_dir), "xml", sep = ".")))

    ## compress
    owd <- getwd()
    setwd(test_dir)
    zip(zipfile = zipname <- paste(name, "zip", sep = "."), files = list.files(test_dir))
    setwd(owd)

    ## copy the final .zip file
    file.copy(file.path(test_dir, zipname), file.path(dir, zipname), overwrite = TRUE)

    invisible(NULL)
  }
}


## generic WU question writer function
write_lops_question <- function(x, dir, ...)
{
  UseMethod("write_lops_question")
}


## write WU mchoice question
## only question type at the moment

## FIXME: add identical schoice method?

write_lops_question.mchoice <- write_lops_question.schoice <- function(x, dir, ...)
{
  args <- list(...)
  time <- gsub("-", "", paste(as.character(Sys.Date()), make_id(10), sep = ""))
  name <- if(is.null(x$metainfo$name)) paste("mchoice.lops", time, sep = "") else paste(x$metainfo$name, time)
  name <- paste("sd_excs", gsub(" ", "_", name), sep = "_")
  presentation <- if(is.null(args$presentation)) "fix" else args$presentation

  xml <- c(
    '<?xml version=\'1.0\' encoding=\'UTF-8\'?>',
    '<learning_resources>',
    paste('<exercise area="sd" restype="excs" shortname="', name, '">', sep = ""),
    '<metadata>',
    paste('<title>', if(is.null(x$metainfo$title)) name else x$metainfo$title, '</title>', sep = ""),
    '<time>1</time>',
    paste('<author>', if(is.null(args$author)) 'Mr. Exam' else args$author, '</author>', sep = ""),
    paste('<language>', if(is.null(args$language)) 'en' else args$language, '</language>', sep = ""),
    '<category restype="csp" shortname="sd_root_concept" area="sd"/>',
    '</metadata>',
    '<question_data>',
    paste('<multiplechoice presentation="', presentation, '" points="100">', sep = '')
  )

  ## insert the possible answers
  for(i in seq_along(x$questionlist)) {
    xml <- c(xml,
      paste('<answer value="', if(x$metainfo$solution[i]) 'true' else 'false', '">', sep = ""),
      '<answer_text>',
      x$questionlist[i],
      '</answer_text>',
      '</answer>'
    )
  }

  ## insert the question and finish
  xml <- c(xml,
    '<problem_text>',
    x$question,
    '</problem_text>',
    '<feedback/>',
    '<assessment schema="wi_score_multiplechoice"/>',
    '</multiplechoice>',
    '</question_data>',
    '</exercise>',
    '</learning_resources>'
  )

  ## path replacements
  xml <- gsub(paste(attr(x$supplements, "dir"), .Platform$file.sep, sep = ""), "", xml)

  ## <img and src= tag replacements, otherwise there will
  ## be an error message when uploading
  xml <- gsub('<img', '<image', xml)
  xml <- gsub('src="', 'name="', xml)
  
  ## create the question directory and store all files
  dir.create(exdir <- file.path(path.expand(dir), name))
  writeLines(xml, file.path(exdir, paste(name, "xml", sep = ".")))
  if(!is.null(x$supplements) && length(x$supplements))
    for(i in x$supplements)
      if(!grepl("solution", basename(i)))
        file.copy(i, file.path(exdir, basename(i)))

  ## prepare return string
  points <- if(is.null(x$metainfo$points)) 1 else x$metainfo$points
  rval <- paste('<excs area="sd" restype="excs" shortname="', name,
    '" points="', points, '"/>', sep = "")
  
  rval
}
