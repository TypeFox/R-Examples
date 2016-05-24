## generate exams in Moodle XML format
## http://docs.moodle.org/en/Moodle_XML_format
exams2moodle <- function(file, n = 1L, nsamp = NULL, dir = ".",
  name = NULL, quiet = TRUE, edir = NULL, tdir = NULL, sdir = NULL, verbose = FALSE,
  resolution = 100, width = 4, height = 4, encoding = "", 
  iname = TRUE, stitle = NULL, testid = FALSE, zip = FALSE,
  num = NULL, mchoice = NULL, schoice = mchoice, string = NULL, cloze = NULL,
  points = NULL, rule = NULL, pluginfile = TRUE,
  converter = NULL, ...)
{
  ## default converter is "ttm" if all exercises are Rnw, otherwise "pandoc"
  if(is.null(converter)) {
    converter <- if(any(tolower(tools::file_ext(unlist(file))) == "rmd")) "pandoc" else "ttm"
  }
  ## set up .html transformer
  htmltransform <- make_exercise_transform_html(converter = converter, ..., base64 = !pluginfile)

  ## generate the exam
  if(encoding == "") encoding <- "UTF-8"
  exm <- xexams(file, n = n, nsamp = nsamp,
   driver = list(
       sweave = list(quiet = quiet, pdf = FALSE, png = TRUE,
         resolution = resolution, width = width, height = height,
         encoding = encoding),
       read = NULL, transform = htmltransform, write = NULL),
     dir = dir, edir = edir, tdir = tdir, sdir = sdir, verbose = verbose)

  ## get the possible moodle question body functions and options
  moodlequestion = list(num = num, mchoice = mchoice, schoice = schoice, cloze = cloze, string = string)

  for(i in c("num", "mchoice", "schoice", "cloze", "string")) {
    if(is.null(moodlequestion[[i]])) moodlequestion[[i]] <- list()
    if(is.list(moodlequestion[[i]])) {
      if(is.null(moodlequestion[[i]]$eval))
        moodlequestion[[i]]$eval <- list("partial" = TRUE, "rule" = rule)
      if(is.list(moodlequestion[[i]]$eval)) {
        if(!moodlequestion[[i]]$eval$partial) stop("Moodle can only process partial credits!")
        if(i == "cloze" & is.null(moodlequestion[[i]]$eval$rule))
          moodlequestion[[i]]$eval$rule <- "none"
      }
      moodlequestion[[i]] <- do.call("make_question_moodle", moodlequestion[[i]])
    }
    if(!is.function(moodlequestion[[i]])) stop(sprintf("wrong specification of %s", sQuote(i)))
  }

  ## create a temporary directory
  dir <- path.expand(dir)
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  } else {
    tdir <- path.expand(tdir)
  }
  if(!file.exists(tdir)) dir.create(tdir)

  ## obtain the number of exams and questions
  nx <- length(exm)
  nq <- length(exm[[1L]])

  ## create a name
  if(is.null(name)) name <- "moodlequiz"

  ## function for internal ids
  make_test_ids <- function(n, type = c("test", "section", "item"))
  {
    switch(type,
      "test" = if(testid) paste(name, make_id(9), sep = "_") else name,
      paste(type, formatC(1:n, flag = "0", width = nchar(n)), sep = "_")
    )
  }

  ## generate the test id
  test_id <- make_test_ids(type = "test")

  ## create the directory where the test is stored
  dir.create(test_dir <- file.path(tdir, name))

  exsecs <- rep(NA, length = nq)
  if(!is.null(stitle)) {
    if((ks <- length(stitle)) > nq) stop("more section titles than exercises specified!")
    exsecs[1:ks] <- stitle
  }

  ## points setting
  if(!is.null(points))
    points <- rep(points, length.out = nq)

  ## encoding
  enc <- gsub("-", "", tolower(encoding), fixed = TRUE)
  if(enc %in% c("iso8859", "iso88591")) enc <- "latin1"
  if(enc == "iso885915") enc <- "latin9"
  charset <- encoding
  if(enc == "utf8")
    charset <- "UTF-8"
  if(enc == "latin1")
    charset <- "ISO-8859-1"
  if(enc == "latin9")
    charset <- "ISO-8859-15"

  ## start the quiz .xml
  xml <- c(paste('<?xml version="1.0" encoding="', charset, '"?>', sep = ''), '<quiz>\n')

  ## cycle through all questions and samples
  for(j in 1:nq) {
    ## search for \exsection{}
    exsec <- if(is.null(exm[[1]][[j]]$metainfo$section)) {
      paste("Exercise", formatC(j, flag = "0", width = nchar(nq)))
    } else exm[[1]][[j]]$metainfo$section

    ## if specified, overule the section
    if(!is.na(exsecs[j])) exsec <- exsecs[j]

    ## first, create the category tag for the question
    xml <- c(xml,
      '\n<question type="category">',
      '<category>',
      paste('<text>$course$/', if(iname) paste(test_id, '/', sep = '') else NULL, exsec, '</text>', sep = ''),
      '</category>',
      '</question>\n')

    ## create ids for all samples
    sample_ids <- paste(exsec, make_test_ids(nx, type = "sample"), sep = "_")

    ## create the questions
    for(i in 1:nx) {
      ## overule points
      if(!is.null(points)) exm[[i]][[j]]$metainfo$points <- points[[j]]

      ## get the question type
      type <- exm[[i]][[j]]$metainfo$type

      ## attach sample id to metainfo
      exm[[i]][[j]]$metainfo$id <- paste(sample_ids[i], type, sep = "_")

      ## add sample and question number to name
      sqname <- ""
      if(nx>1L) sqname <- paste("R", formatC(i, flag = "0", width = nchar(nx)), " ", sep = "")
      exm[[i]][[j]]$metainfo$name <-
          paste(sqname,
                "Q", formatC(j, flag = "0", width = nchar(nq)),
                " : ",
                if(!is.null(exm[[i]][[j]]$metainfo$title)) {
                    exm[[i]][[j]]$metainfo$title
                } else exm[[i]][[j]]$metainfo$file,
                sep="")
      ## create the .xml
      question_xml <- moodlequestion[[type]](exm[[i]][[j]])

      ## include supplements using base64 encoding, with either moodle's
      ## pluginfile mechanism or data URIs
      if(length(exm[[i]][[j]]$supplements) > 0) {
        for(si in seq_along(exm[[i]][[j]]$supplements)) {
          if(any(grepl(f <- basename(exm[[i]][[j]]$supplements[si]), question_xml))) {
            if(isTRUE(pluginfile)) {
              newfn   <- paste0("@@PLUGINFILE@@/", f)
              href    <- paste0("\"", f,"\"")
              newhref <- paste0("\"", newfn,"\"")
              filetag <- paste0("<file name=\"", f, "\" encoding=\"base64\">",
                                base64enc::base64encode(exm[[i]][[j]]$supplements[si]),
                                "</file>")

              # Prepend @@PLUGINFILE@@ to link target
              question_xml <- gsub(href, newhref, question_xml, fixed = TRUE)

              # Insert base64 encoded file at the end of <questiontext>
              idx <- which(grepl(newhref, question_xml, fixed = TRUE))
              textend <- which(grepl("</text>", question_xml, fixed = TRUE))
              textend <- head(textend[textend > idx], 1)

              question_xml <- append(question_xml, filetag, after = textend)
            } else {
              question_xml <- gsub(paste(f, '"', sep = ''),
                paste(fileURI(exm[[i]][[j]]$supplements[si]), '"', sep = ''),
                question_xml, fixed = TRUE)
            }
          }
        }
      }

      ## add question to quiz .xml
      xml <- c(xml, question_xml)
    }
  }

  ## finish the quiz
  xml <- c(xml, '</quiz>')

  ## write to dir
  writeLines(xml, file.path(test_dir, paste(name, "xml", sep = ".")))

  ## compress
  if(zip) {
    owd <- getwd()
    setwd(test_dir)
    zip(zipfile = zipname <- paste(name, "zip", sep = "."), files = list.files(test_dir))
    setwd(owd)
  } else zipname <- list.files(test_dir)

  ## copy the final .zip file
  file.copy(file.path(test_dir, zipname), dir, recursive = TRUE)

  ## assign test id as an attribute
  attr(exm, "test_id") <- test_id

  invisible(exm)
}


## Moodle question constructor function (originally for Moodle 2.3)
make_question_moodle <-
make_question_moodle23 <- function(name = NULL, solution = TRUE, shuffle = FALSE, penalty = 0,
  answernumbering = "abc", usecase = FALSE, cloze_mchoice_display = "MULTICHOICE",
  truefalse = c("True", "False"), enumerate = TRUE, abstention = NULL,
  eval = list(partial = TRUE, negative = FALSE, rule = "false2"))
{
  function(x) {
    ## how many points?
    points <- if(is.null(x$metainfo$points)) 1 else x$metainfo$points

    ## choice policy
    eval <- if(!all(names(exams_eval()) %in% names(eval))) {
      if(x$metainfo$type == "cloze" & is.null(eval$rule))
        eval$rule <- "none"
      do.call("exams_eval", eval)
    } else eval

    ## match question type
    type <- switch(x$metainfo$type,
      "num" = "numerical",
      "mchoice" = "multichoice",
      "schoice" = "multichoice",
      "cloze" = "cloze",
      "string" = "shortanswer"
    )

    ## question name
    if(is.null(name)) name <- x$metainfo$name

    ## extra abstention option
    if(is.null(abstention)) abstention <- x$metainfo$abstention
    if(is.null(abstention) || identical(abstention, FALSE)) abstention <- ""
    if(isTRUE(abstention)) abstention <- "Abstention"

    ## start the question xml
    xml <- c(
      paste('\n<question type="', type, '">', sep = ''),
      '<name>',
      paste('<text>', name, '</text>'),
      '</name>',
      '<questiontext format="html">',
      '<text><![CDATA[<p>', if(type != "cloze") x$question else '##QuestionText##', '</p>]]></text>',
      '</questiontext>'
    )

    ## insert the solution
    if((length(x$solution) | (nsol <- length(x$solutionlist))) && solution) {
      xml <- c(xml,
        '<generalfeedback format="html">',
        '<text><![CDATA[<p>', x$solution,
        if(!type %in% c("mchoice", "schoice") && nsol) {
          g <- rep(seq_along(x$metainfo$solution), sapply(x$metainfo$solution, length))
          soll <- sapply(split(x$solutionlist, g), paste, collapse = " / ")
          c(if(enumerate) '<ol type = "a">' else '</br>',
            paste(if(enumerate) "<li>" else NULL, soll, if(enumerate) "</li>" else NULL),
            if(enumerate) '</ol>' else NULL)
        } else NULL,
        '</p>]]></text>',
        '</generalfeedback>'
      )
    }

    ## penalty and points
    if(type == "cloze")
      points <- rep(points, length.out = length(x$metainfo$solution))
    xml <- c(xml,
      paste('<penalty>', penalty, '</penalty>', sep = ''),
      paste('<defaultgrade>', sum(points), '</defaultgrade>', sep = '')
    )

    ## multiple choice processing
    if(type == "multichoice") {
      xml <- c(xml,
        paste('<shuffleanswers>', if(shuffle) 'true' else 'false', '</shuffleanswers>', sep = ''),
        paste('<single>', if(x$metainfo$type == "schoice") 'true' else 'false', '</single>', sep = ''),
        paste('<answernumbering>', answernumbering, '</answernumbering>', sep = '')
      )

      frac <- as.integer(x$metainfo$solution)
      pv <- eval$pointvec(paste(frac, sep = "", collapse = ""))
      pv[pv == -Inf] <- 0 ## FIXME: exams_eval() return -Inf when rule = "none"?

      frac[x$metainfo$solution] <- pv["pos"]
      frac[!x$metainfo$solution] <- pv["neg"]
      frac <- moodlePercent(frac)

      for(i in seq_along(x$questionlist)) {
        xml <- c(
          xml,
          paste('<answer fraction="', frac[i], '" format="html">', sep = ''),
          '<text><![CDATA[<p>', x$questionlist[i], '</p>]]></text>',
          if(!is.null(x$solutionlist)) {
            c('<feedback format="html">',
            '<text><![CDATA[<p>', x$solutionlist[i], '</p>]]></text>',
            '</feedback>')
          } else NULL,
          '</answer>'
        )
      }

      ## add abstention option (if any)
      if(abstention != "") {
        xml <- c(xml,
          '<answer fraction="0" format="html">',
          '<text><![CDATA[<p>',
          abstention,
          '</p>]]></text>',
          '</answer>'
        )
      }
    }

    ## numeric question processing
    if(type == "numerical") {
      xml <- c(xml,
        '<answer fraction="100" format="moodle_auto_format">',
        paste('<text>', x$metainfo$solution, '</text>', sep = ''),
        paste('<tolerance>', max(x$metainfo$tolerance), '</tolerance>', sep = ''),
        '</answer>'
      )
    }

    ## string questions
    if(type == "shortanswer") {
      xml <- c(xml,
        paste('<usecase>', usecase * 1, '</usecase>', sep = ''),
        '<answer fraction="100" format="moodle_auto_format">',
        '<text>', x$metainfo$solution, '</text>',
        '</answer>'
      )
    }

    ## cloze type questions
    if(type == "cloze") {
      ## how many questions
      solution <- if(!is.list(x$metainfo$solution)) {
        list(x$metainfo$solution)
      } else x$metainfo$solution
      n <- length(solution)

      xml[grep('<defaultgrade>', xml, fixed = TRUE)] <- paste('<defaultgrade>', sum(points),
        '</defaultgrade>', sep = '')

      questionlist <- if(!is.list(x$questionlist)) {
        if(x$metainfo$type == "cloze") as.list(x$questionlist) else list(x$questionlist)
      } else x$questionlist
      if(length(questionlist) < 1) questionlist <- NULL

      ## split id for the questionlist
      sid <- unlist(sapply(1:n, function(i) rep(i, length(solution[[i]]))))

      ## tolerance of numerical questions
      tol <- if(!is.list(x$metainfo$tolerance)) {
        if(x$metainfo$type == "cloze") as.list(x$metainfo$tolerance) else list(x$metainfo$tolerance)
      } else x$metainfo$tolerance
      tol <- rep(tol, length.out = n)

      ## optionally fix the num answer field width
      ## by supplying an additional wrong answer
      numwidth <- if(is.null(x$metainfo$numwidth)) FALSE else TRUE
      if(numwidth) {
        nums <- NULL
        for(i in 1:n) {
          ql <- if(is.null(questionlist)) "" else questionlist[sid == i]
          k <- length(ql)
          if(x$metainfo$clozetype[i] == "num") {
            for(j in 1:k) {
              nums <- rbind(nums,
                c(solution[[i]][j] - max(tol[[i]]),
                solution[[i]][j] + max(tol[[i]])))
            }
          }
        }
        if(!is.null(nums)) {
          if(is.logical(x$metainfo$numwidth)) {
            fnums <- format(as.numeric(nums))
          } else {
            fnums <- if(!is.character(x$metainfo$numwidth)) {
              paste(rep("1", length = as.integer(x$metainfo$numwidth)), sep = "", collapse = "")
            } else x$metainfo$numwidth
          }
          num_w <- max(unlist(sapply(fnums, nchar)))
          do <- TRUE
          while(do) {
            fnums <- make_id(num_w)
            tolcheck <- NULL
            for(i in 1:nrow(nums)) {
              tolcheck <- c(tolcheck, fnums >= nums[i, 1] & fnums <= nums[i, 2])
            }
            do <- (fnums %in% nums) & any(tolcheck)
          }
        }
      }

      ## cycle through all questions
      qtext <- NULL; inum <- 1
      for(i in 1:n) {
        ql <- if(is.null(questionlist)) "" else questionlist[sid == i]
        k <- length(ql)
        tmp <- NULL
        if(grepl("choice", x$metainfo$clozetype[i])) {
          tmp <- paste('{', points[i], ':', cloze_mchoice_display, ':', sep = '')

          frac <- solution[[i]]
          if(length(frac) < 2)
            frac <- c(frac, !frac)
          frac2 <- frac
          pv <- eval$pointvec(frac)
          frac[frac2] <- pv["pos"]
          frac[!frac2] <- pv["neg"]
          p <- moodlePercent(frac)

          if(k < 2) {
            tmp <- paste(ql, tmp)
            p <- paste('%', p, '%', sep = '')
            p[2] <- paste('~', p[2], sep = '')
            ql <- paste(p, truefalse[rev(frac2 + 1)], sep = '', collapse = '')
          } else {
            ql2 <- NULL
            for(j in 1:k)
              ql2 <- paste(ql2, if(j > 1) '~' else NULL, paste('%',  p[j], '%', sep = ''), ql[j], sep = '')
            ql <- ql2
          }
          tmp <- paste(tmp, ql, sep = '')
          tmp <- paste(tmp, '}', sep = '')
        }
        if(x$metainfo$clozetype[i] == "num") {
          for(j in 1:k) {
            tmp <- c(tmp, paste(ql[j], ' {', points[i], ':NUMERICAL:=', solution[[i]][j],
              ':', max(tol[[i]]), if(numwidth) paste('~%0%', fnums, ":0", sep = '') else NULL,
              '}', sep = ''))
          }
        }
        if(x$metainfo$clozetype[i] == "string") {
          for(j in 1:k) {
            tmp <- c(tmp, paste(ql[j], ' {', points[i], ':SHORTANSWER:%100%', solution[[i]][j],
              if(!usecase) paste('~%100%', tolower(solution[[i]][j]), sep = '') else NULL,
              '}', sep = ''))
          }
        }
        if(x$metainfo$clozetype[i] == "verbatim") {
          for(j in 1:k) {
            tmp <- c(tmp, paste0(ql[j], ' {', points[i], solution[[i]][j], '}'))
          }
        }

        ## FIXME, there is a NULL when using boxhist2?
        tmp <- gsub('NULL', '', tmp)

        ## insert in ##ANSWERi## tag
        if(any(grepl(ai <- paste("##ANSWER", i, "##", sep = ""), x$question, fixed = TRUE))) {
          x$question <- gsub(ai, paste(tmp, collapse = ", "), x$question, fixed = TRUE)
        } else qtext <- c(qtext, tmp)
      }
      if(!is.null(qtext) & enumerate)
        qtext <- c('<ol type = "a">', paste('<li>', qtext, '</li>'), '</ol>')
      qtext <- c(x$question, qtext)
      xml <- gsub('##QuestionText##', paste(qtext, collapse = "\n"), xml)
    }

    ## end the question
    xml <- c(xml, '</question>\n')

    ## path replacements
    xml <- gsub(paste(attr(x$supplements, "dir"), .Platform$file.sep, sep = ""), "", xml)

    xml
  }
}


## "Numbers" Moodle currently accepts as fraction value
## for mchoice items
moodleFractions <- c(100,90,83.33333,80,75,70,
                     66.66667,60,50,40,
                     33.33333,30,25,20,16.66667,
                     14.28571, 12.5,11.11111, 10,5)

## Convert a number in [0, 1] to one of the percentages
## above if the difference is less then 1
moodlePercent <- function(p)
{
  p <- 100 * p
  z <- abs(outer(abs(p), moodleFractions, "-"))
  mp <- moodleFractions[max.col(-z)] * sign(p)
  if(any(abs(mp - p) > 1))
    stop("Percentage not in list of moodle fractions")
  as.character(mp)
}

