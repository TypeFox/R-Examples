## create IMS QTI 1.2 .xml files
## specifications and examples available at:
## http://www.imsglobal.org/question/qtiv1p2/imsqti_asi_bindv1p2.html
## http://www.imsglobal.org/question/qtiv1p2/imsqti_asi_bestv1p2.html#1466669
exams2qti12 <- function(file, n = 1L, nsamp = NULL, dir = ".",
  name = NULL, quiet = TRUE, edir = NULL, tdir = NULL, sdir = NULL, verbose = FALSE,
  resolution = 100, width = 4, height = 4, encoding  = "",
  num = NULL, mchoice = NULL, schoice = mchoice, string = NULL, cloze = NULL,
  template = "qti12",
  duration = NULL, stitle = "Exercise", ititle = "Question",
  adescription = "Please solve the following exercises.",
  sdescription = "Please answer the following question.", 
  maxattempts = 1, cutvalue = 0, solutionswitch = TRUE, zip = TRUE,
  points = NULL, eval = list(partial = TRUE, negative = FALSE),
  converter = NULL, ...)
{
  ## default converter is "ttm" if all exercises are Rnw, otherwise "pandoc"
  if(is.null(converter)) {
    converter <- if(any(tolower(tools::file_ext(unlist(file))) == "rmd")) "pandoc" else "ttm"
  }
  ## set up .html transformer
  htmltransform <- make_exercise_transform_html(converter = converter, ...)

  ## generate the exam
  is.xexam <- FALSE
  if(is.list(file)) {
    if(any(grepl("exam1", names(file))))
      is.xexam <- TRUE
  }
  if(!is.xexam) {
    exm <- xexams(file, n = n, nsamp = nsamp,
      driver = list(
        sweave = list(quiet = quiet, pdf = FALSE, png = TRUE,
          resolution = resolution, width = width, height = height,
          encoding = encoding),
        read = NULL, transform = htmltransform, write = NULL),
      dir = dir, edir = edir, tdir = tdir, sdir = sdir, verbose = verbose)
  } else {
    exm <- file
    rm(file)
  }

  ## start .xml assessement creation
  ## get the possible item body functions and options  
  itembody <- list(num = num, mchoice = mchoice, schoice = schoice, cloze = cloze, string = string)

  for(i in c("num", "mchoice", "schoice", "cloze", "string")) {
    if(is.null(itembody[[i]])) itembody[[i]] <- list()
    if(is.list(itembody[[i]])) {
      if(is.null(itembody[[i]]$eval))
        itembody[[i]]$eval <- eval
      if(i == "cloze" & is.null(itembody[[i]]$eval$rule))
        itembody[[i]]$eval$rule <- "none"
      itembody[[i]] <- do.call("make_itembody_qti12", itembody[[i]])
    }
    if(!is.function(itembody[[i]])) stop(sprintf("wrong specification of %s", sQuote(i)))
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

  ## the package directory
  pkg_dir <- find.package("exams")

  ## get the .xml template
  template <- path.expand(template)
  template <- ifelse(
    tolower(substr(template, nchar(template) - 3L, nchar(template))) != ".xml",
    paste(template, ".xml", sep = ""), template)
  template <- ifelse(file.exists(template),
    template, file.path(pkg_dir, "xml", basename(template)))
  if(!all(file.exists(template))) {
    stop(paste("The following files cannot be found: ",
      paste(basename(template)[!file.exists(template)], collapse = ", "), ".", sep = ""))
  }
  xml <- readLines(template[1L])

  ## check template for section and item inclusion
  ## extract the template for the assessement, sections and items
  if(length(section_start <- grep("<section ident", xml, fixed = TRUE)) != 1L ||
    length(section_end <- grep("</section>", xml, fixed = TRUE)) != 1L) {
    stop(paste("The XML template", template,
      "must contain exactly one opening and closing <section> tag!"))
  }
  section <- xml[section_start:section_end]
  if(length(item_start <- grep("<item ident", section, fixed = TRUE)) != 1 ||
    length(item_end <- grep("</item>", section, fixed = TRUE)) != 1) {
    stop(paste("The XML template", template,
      "must contain exactly one opening and closing <item> tag!"))
  }
  xml <- c(xml[1:(section_start - 1)], "##TestSections##", xml[(section_end + 1):length(xml)])
  item <- section[item_start:item_end]
  section <- section[1:(item_start - 1)]

  ## obtain the number of exams and questions
  nx <- length(exm)
  nq <- if(!is.xexam) length(exm[[1L]]) else length(exm)

  ## create a name
  if(is.null(name)) name <- file_path_sans_ext(basename(template))

  ## function for internal ids
  make_test_ids <- function(n, type = c("test", "section", "item"))
  {
    switch(type,
      "test" = paste(name, make_id(9), sep = "_"),
      paste(type, formatC(1:n, flag = "0", width = nchar(n)), sep = "_")
    )
  }

  ## generate the test id
  test_id <- make_test_ids(type = "test")

  ## create section ids
  sec_ids <- paste(test_id, make_test_ids(nq, type = "section"), sep = "_")

  ## create section/item titles and section description
  if(is.null(stitle)) stitle <- ""
  stitle <- rep(stitle, length.out = nq)
  if(!is.null(ititle)) ititle <- rep(ititle, length.out = nq)
  if(is.null(adescription)) adescription <- ""
  if(is.null(sdescription)) sdescription <- ""
  sdescription <- rep(sdescription, length.out = nq)

  ## points setting
  if(!is.null(points))
    points <- rep(points, length.out = nq)

  ## create the directory where the test is stored
  dir.create(test_dir <- file.path(tdir, name))

  ## cycle through all exams and questions
  ## similar questions are combined in a section,
  ## questions are then sampled from the sections
  items <- sec_xml <- NULL
  for(j in 1:nq) {
    ## first, create the section header
    sec_xml <- c(sec_xml, gsub("##SectionId##", sec_ids[j], section, fixed = TRUE))

    ## insert a section title -> exm[[1]][[j]]$metainfo$name?
    sec_xml <- gsub("##SectionTitle##", stitle[j], sec_xml, fixed = TRUE)

    ## insert a section description -> exm[[1]][[j]]$metainfo$description?
    sec_xml <- gsub("##SectionDescription##", sdescription[j], sec_xml, fixed = TRUE)

    ## special handler
    if(is.xexam) nx <- length(exm[[j]])

    ## create item ids
    item_ids <- paste(sec_ids[j], make_test_ids(nx, type = "item"), sep = "_")

    ## now, insert the questions
    for(i in 1:nx) {
      ## special handling of indices
      if(is.xexam) {
        if(i < 2)
          jk <- j
        j <- i
        i <- jk
      }

      ## overule points
      if(!is.null(points)) exm[[i]][[j]]$metainfo$points <- points[[j]]

      ## get and insert the item body
      type <- exm[[i]][[j]]$metainfo$type

      ## create an id
      iname <- paste(item_ids[if(is.xexam) j else i], type, sep = "_")

      ## attach item id to metainfo
      exm[[i]][[j]]$metainfo$id <- iname

      ibody <- gsub("##ItemBody##",
        paste(thebody <- itembody[[type]](exm[[i]][[j]]), collapse = "\n"),
        item, fixed = TRUE)

      ## insert possible solution
      enumerate <- attr(thebody, "enumerate")
      if(is.null(enumerate)) enumerate <- FALSE
      xsolution <- exm[[i]][[j]]$solution
      if(!is.null(exm[[i]][[j]]$solutionlist)) {
        if(!all(is.na(exm[[i]][[j]]$solutionlist))) {
          xsolution <- c(xsolution, if(length(xsolution)) "<br />" else NULL)
          if(enumerate) xsolution <- c(xsolution, '<ol type = "a">')
          if(exm[[i]][[j]]$metainfo$type == "cloze") {
            g <- rep(seq_along(exm[[i]][[j]]$metainfo$solution), sapply(exm[[i]][[j]]$metainfo$solution, length))
            ql <- sapply(split(exm[[i]][[j]]$questionlist, g), paste, collapse = " / ")
            sl <- sapply(split(exm[[i]][[j]]$solutionlist, g), paste, collapse = " / ")
          } else {
            ql <- exm[[i]][[j]]$questionlist
            sl <- exm[[i]][[j]]$solutionlist
          }
          nsol <- length(ql)
          xsolution <- c(xsolution, paste(if(enumerate) rep('<li>', nsol) else NULL,
            ql, if(length(exm[[i]][[j]]$solutionlist)) "<br />" else NULL,
            sl, if(enumerate) rep('</li>', nsol) else NULL))
          if(enumerate) xsolution <- c(xsolution, '</ol>')
        }
      }

      ibody <- gsub("##ItemSolution##", paste(xsolution, collapse = "\n"), ibody, fixed = TRUE)

      ## insert an item id
      ibody <- gsub("##ItemId##", iname, ibody)

      ## insert an item title
      ibody <- gsub("##ItemTitle##",
        if(is.null(ititle)) exm[[i]][[j]]$metainfo$name else ititle[j],
        ibody, fixed = TRUE)

      ## copy supplements
      if(length(exm[[i]][[j]]$supplements)) {
        if(!file.exists(media_dir <- file.path(test_dir, "media")))
          dir.create(media_dir)
        sj <- 1
        while(file.exists(file.path(media_dir, sup_dir <- paste("supplements", sj, sep = "")))) {
          sj <- sj + 1
        }
        dir.create(ms_dir <- file.path(media_dir, sup_dir))
        for(si in seq_along(exm[[i]][[j]]$supplements)) {
          file.copy(exm[[i]][[j]]$supplements[si],
            file.path(ms_dir, f <- basename(exm[[i]][[j]]$supplements[si])))
          if(any(grepl(dirname(exm[[i]][[j]]$supplements[si]), ibody))) {
            ibody <- gsub(dirname(exm[[i]][[j]]$supplements[si]),
              file.path('media', sup_dir), ibody, fixed = TRUE)
          } else {
            if(any(grepl(f, ibody))) {
              ibody <- gsub(paste(f, '"', sep = ''),
                paste('media/', sup_dir, '/', f, '"', sep = ''), ibody, fixed = TRUE)
            }
          }
        }
      }

      ## include body in section
      sec_xml <- c(sec_xml, ibody, "")
    }

    ## close the section
    sec_xml <- c(sec_xml, "", "</section>")
  }

  ## process duration to P0Y0M0DT0H1M35S format
  if(!is.null(duration)) {
    dursecs <- round(duration * 60)
    dur <- dursecs %/% 86400 ## days
    dursecs <- dursecs - dur * 86400
    duration <- paste("P0Y0M", dur, "DT", sep = "")
    dur <- dursecs %/% 3600 ## hours
    dursecs <- dursecs - dur * 3600
    duration <- paste(duration, dur, "H", sep = "")
    dur <- dursecs %/% 60 ## minutes
    dursecs <- dursecs - dur * 60
    duration <- paste("<duration>", duration, dur, "M", dursecs, "S", "</duration>", sep = "")
  } else {
    duration <- ""
  }

  ## process cutvalue/maximal number of attempts
  make_integer_tag <- function(x, type, default = 1) {
    if(is.null(x)) x <- Inf
    x <- round(as.numeric(x))
    if(x < default) {
      warning(paste("invalid ", type, " specification, ", type, "=", default, " used", sep = ""))
      x <- default
    }
    if(is.finite(x)) sprintf("%s=\"%i\"", type, x) else ""
  }
  maxattempts <- make_integer_tag(maxattempts, type = "maxattempts", default = 1)
  cutvalue <- make_integer_tag(cutvalue, type = "cutvalue", default = 0)

  ## finalize the test xml file, insert ids/titles, sections, and further control details
  feedbackswitch <- FALSE ## currently hard-coded
  hintswitch <- FALSE
  xml <- gsub("##TestIdent##", test_id, xml)
  xml <- gsub("##TestTitle##", name, xml)
  xml <- gsub("##TestDuration##", duration, xml)
  xml <- gsub("##TestSections##", paste(sec_xml, collapse = "\n"), xml)
  xml <- gsub("##MaxAttempts##", maxattempts, xml)
  xml <- gsub("##CutValue##", cutvalue, xml)
  xml <- gsub("##FeedbackSwitch##", if(feedbackswitch) "Yes" else "No", xml)
  xml <- gsub("##HintSwitch##",     if(hintswitch)     "Yes" else "No", xml)
  xml <- gsub("##SolutionSwitch##", if(solutionswitch) "Yes" else "No", xml)
  xml <- gsub("##AssessmentDescription##", adescription, xml)

  ## write to dir
  writeLines(xml, file.path(test_dir, "qti.xml"))

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


## QTI 1.2 item body constructor function
## includes item <presentation> and <resprocessing> tags
make_itembody_qti12 <- function(rtiming = FALSE, shuffle = FALSE, rshuffle = shuffle,
  minnumber = NULL, maxnumber = NULL, defaultval = NULL, minvalue = NULL,
  maxvalue = NULL, cutvalue = NULL, enumerate = TRUE, digits = NULL, tolerance = is.null(digits),
  maxchars = 12, eval = list(partial = TRUE, negative = FALSE), fix_num = TRUE)
{
  function(x) {
    ## how many points?
    points <- if(is.null(x$metainfo$points)) 1 else x$metainfo$points

    ## how many questions
    solution <- if(!is.list(x$metainfo$solution)) {
      list(x$metainfo$solution)
    } else x$metainfo$solution
    n <- length(solution)

    questionlist <- if(!is.list(x$questionlist)) {
      if(x$metainfo$type == "cloze") {
        g <- rep(seq_along(x$metainfo$solution), sapply(x$metainfo$solution, length))
        split(x$questionlist, g)
      } else list(x$questionlist)
    } else x$questionlist
    if(length(questionlist) < 1) questionlist <- NULL

    tol <- if(!is.list(x$metainfo$tolerance)) {
      if(x$metainfo$type == "cloze") as.list(x$metainfo$tolerance) else list(x$metainfo$tolerance)
    } else x$metainfo$tolerance
    tol <- rep(tol, length.out = n)

    q_points <- rep(points, length.out = n)
    if(x$metainfo$type == "cloze")
      points <- sum(q_points)

    ## set question type(s)
    type <- x$metainfo$type
    type <- if(type == "cloze") x$metainfo$clozetype else rep(type, length.out = n)

    ## evaluation policy
    if(is.null(eval) || length(eval) < 1L) eval <- exams_eval()
    if(!is.list(eval)) stop("'eval' needs to specify a list of partial/negative/rule")
    eval <- eval[match(c("partial", "negative", "rule"), names(eval), nomatch = 0)]
    if(x$metainfo$type %in% c("num", "string")) eval$partial <- FALSE
    if(x$metainfo$type == "cloze" & is.null(eval$rule)) eval$rule <- "none"
    eval <- do.call("exams_eval", eval) ## always re-call exams_eval

    ## character fields
    maxchars <- if(is.null(x$metainfo$maxchars)) {
        if(length(maxchars) < 2) {
           c(maxchars, NA, NA)
        } else maxchars[1:3]
    } else x$metainfo$maxchars
    if(!is.list(maxchars))
      maxchars <- list(maxchars)
    maxchars <- rep(maxchars, length.out = n)
    for(j in seq_along(maxchars)) {
      if(length(maxchars[[j]]) < 2)
        maxchars[[j]] <- c(maxchars[[j]], NA, NA)
    }

    ## start item presentation
    ## and insert question
    xml <- c(
      '<presentation>',
      '<flow>',
      if(!is.null(x$question)) {
        c(
          '<material>',
          '<matbreak/>',
          '<mattext texttype="text/html" charset="utf-8"><![CDATA[',
          x$question,
          ']]></mattext>',
          '<matbreak/>',
          '</material>'
        )
      } else NULL
    )

    ## insert responses
    ids <- el <- pv <- list()
    for(i in 1:n) {
      ## get item id
      iid <- x$metainfo$id

      ## generate ids
      ids[[i]] <- list("response" = paste(iid, "RESPONSE", make_id(7), sep = "_"),
        "questions" = paste(iid, make_id(10, length(x$metainfo$solution)), sep = "_"))

      ## evaluate points for each question
      pv[[i]] <- eval$pointvec(solution[[i]])
      pv[[i]]["pos"] <- pv[[i]]["pos"] * q_points[i]
      if(length(grep("choice", type[i])))
        pv[[i]]["neg"] <- pv[[i]]["neg"] * q_points[i]

      ## insert choice type responses
      if(length(grep("choice", type[i]))) {
        xml <- c(xml,
          paste('<response_lid ident="', ids[[i]]$response, '" rcardinality="',
            if(type[i] == "mchoice") "Multiple" else "Single", '" rtiming=',
            if(rtiming) '"Yes"' else '"No"', '>', sep = ''),
          paste('<render_choice shuffle=', if(shuffle) '"Yes">' else '"No">', sep = '')
        )
        for(j in seq_along(solution[[i]])) {
          xml <- c(xml,
            '<flow_label class="List">',
            paste('<response_label ident="', ids[[i]]$questions[j], '" rshuffle="',
              if(rshuffle) 'Yes' else 'No', '">', sep = ''),
            '<material>',
            '<mattext texttype="text/html" charset="utf-8"><![CDATA[',
             paste(if(enumerate) {
               paste(letters[if(x$metainfo$type == "cloze") i else j], ".",
                 if(x$metainfo$type == "cloze" && length(solution[[i]]) > 1) paste(j, ".", sep = "") else NULL,
                 sep = "")
             } else NULL, questionlist[[i]][j]),
            ']]></mattext>',
            '</material>',
            '</response_label>',
            '</flow_label>'
          )
        }

        ## finish response tag
        xml <- c(xml,
          '</render_choice>',
          '</response_lid>'
        )
      }
      if(type[i] == "string" || type[i] == "num") {
        for(j in seq_along(solution[[i]])) {
          soltext <- if(type[i] == "num") {
             if(!is.null(digits)) format(round(solution[[i]][j], digits), nsmall = digits) else solution[[i]][j]
          } else {
            if(!is.character(solution[[i]][j])) format(solution[[i]][j]) else solution[[i]][j]
          }
          xml <- c(xml,
            if(!is.null(questionlist[[i]][j])) {
              c('<material>',
                paste('<mattext><![CDATA[', paste(if(enumerate) {
                  paste(letters[i], ".", sep = '')
                } else NULL, questionlist[[i]][j]), ']]></mattext>', sep = ""),
                '</material>',
                '<material>', '<matbreak/>', '</material>'
              )
            } else NULL,
            paste(if(type[i] == "string") '<response_str ident="' else {
              if(!tolerance | fix_num) '<response_str ident="' else '<response_num ident="'
              }, ids[[i]]$response, '" rcardinality="Single">', sep = ''),
            paste('<render_fib',
              if(!is.na(maxchars[[i]][1])) {
                paste(' maxchars="', max(c(nchar(soltext), maxchars[[i]][1])), '"', sep = '')
              } else NULL,
              if(!is.na(maxchars[[i]][2])) {
                paste(' rows="', maxchars[[i]][2], '"', sep = '')
              } else NULL,
              if(!is.na(maxchars[[i]][3])) {
                paste(' columns="', maxchars[[i]][3], '"', sep = '')
              } else NULL, '>', sep = ''),
            '<flow_label class="Block">',
            paste('<response_label ident="', ids[[i]]$response, '" rshuffle="No"/>', sep = ''),
            '</flow_label>',
            '</render_fib>',
            if(type[i] == "string") '</response_str>' else {
              if(!tolerance | fix_num) '</response_str>' else '</response_num>'
            },
            '<material>', '<matbreak/>', '</material>'
          )
        }
      }
    }

    ## finish presentation
    xml <- c(xml, '</flow>', '</presentation>')

    if(is.null(minvalue)) {  ## FIXME: add additional switch, so negative points don't carry over?!
      if(eval$negative) {
        minvalue <- sum(sapply(pv, function(x) { x["neg"] }))
      } else minvalue <- 0
    }

    ## start resprocessing
    xml <- c(xml,
      '<resprocessing>',
      '<outcomes>',
      paste('<decvar varname="SCORE" vartype="Decimal" defaultval="',
        if(is.null(defaultval)) 0 else defaultval, '" minvalue="',
        if(is.null(minvalue)) 0 else minvalue, '" maxvalue="',
        if(is.null(maxvalue)) points else maxvalue, '" cutvalue="',
        if(is.null(cutvalue)) points else cutvalue, '"/>', sep = ''),
      '</outcomes>')

    correct_answers <- wrong_answers <- correct_num <- wrong_num <- vector(mode = "list", length = n)
    for(i in 1:n) {
      if(length(grep("choice", type[i]))) {
        for(j in seq_along(solution[[i]])) {
          if(solution[[i]][j]) {
            correct_answers[[i]] <- c(correct_answers[[i]],
              paste('<varequal respident="', ids[[i]]$response,
                '" case="Yes">', ids[[i]]$questions[j], '</varequal>', sep = '')
            )
          } else {
            wrong_answers[[i]] <- c(wrong_answers[[i]],
              paste('<varequal respident="', ids[[i]]$response,
                '" case="Yes">', ids[[i]]$questions[j], '</varequal>', sep = '')
            )
          }
        }
      }
      if(type[i] == "string" || type[i] == "num") {
        for(j in seq_along(solution[[i]])) {
          if(type[i] == "string") {
            soltext <- if(!is.character(solution[[i]][j])) {
              format(round(solution[[i]][j], digits), nsmall = digits)
            } else solution[[i]][j]
            correct_answers[[i]] <- c(correct_answers[[i]], paste('<varequal respident="', ids[[i]]$response,
              '" case="No"><![CDATA[', soltext, ']]></varequal>', sep = "")
            )
          } else {
            correct_answers[[i]] <- c(correct_answers[[i]],
              if(!tolerance) {
                paste('<varequal respident="', ids[[i]]$response,
                  '" case="No"><![CDATA[', if(!is.null(digits)) {
                    format(round(solution[[i]][j], digits), nsmall = digits)
                  } else solution[[i]][j],
                  ']]></varequal>', sep = "")
              } else {
                if(fix_num) {
                  correct_num[[i]] <- c(correct_num[[i]],
                    paste('<varequal respident="', ids[[i]]$response,
                      '" case="No"><![CDATA[', if(!is.null(digits)) {
                      format(round(solution[[i]][j], digits), nsmall = digits)
                      } else solution[[i]][j],
                      ']]></varequal>', sep = "")
                  )
                }
                wrong_num[[i]] <- paste(
                  '<and>\n',
                  paste('<vargte respident="', ids[[i]]$response, '"><![CDATA[',
                    solution[[i]][j] - max(tol[[i]]),
                    ']]></vargte>\n', sep = ""),
                  paste('<varlte respident="', ids[[i]]$response, '"><![CDATA[',
                    solution[[i]][j] + max(tol[[i]]),
                    ']]></varlte>\n', sep = ""),
                  '</and>', collapse = '\n', sep = ''
                )
              }
            )
          }
        }
      }
      if(!is.null(correct_answers[[i]])) {
        attr(correct_answers[[i]], "points") <- pv[[i]]
        attr(correct_answers[[i]], "type") <- type[i]
      }
      if(!is.null(wrong_answers[[i]]))
        attr(wrong_answers[[i]], "points") <- pv[[i]]
    }

    ## delete NULL list elements
    correct_answers <- delete.NULLs(correct_answers)
    wrong_answers <- delete.NULLs(wrong_answers) 
    correct_num <- unlist(delete.NULLs(correct_num))
    wrong_num <- delete.NULLs(wrong_num)
    if(length(wrong_num)) {
      wrong_num <- sapply(wrong_num, function(x) {
        paste('<not>', x, '</not>', collapse = '\n')
      })
      wrong_num <- unlist(wrong_num)
    }



    ## partial points
    if(eval$partial) {
      if(length(correct_answers)) {
        for(i in seq_along(correct_answers)) {
          for(j in correct_answers[[i]]) {
            xml <- c(xml,
              '<respcondition continue="Yes" title="Mastery">',
              '<conditionvar>',
              j,
              '</conditionvar>',
              paste('<setvar varname="SCORE" action="Add">',
                attr(correct_answers[[i]], "points")["pos"], '</setvar>', sep = ''),
              '</respcondition>'
            )
          }
        }
      }
      if(length(wrong_answers)) {
        for(i in seq_along(wrong_answers)) {
          for(j in wrong_answers[[i]]) {
            xml <- c(xml,
              '<respcondition continue="Yes" title="Fail">',
              '<conditionvar>',
              j,
              '</conditionvar>',
              paste('<setvar varname="SCORE" action="Add">',
                attr(wrong_answers[[i]], "points")["neg"], '</setvar>', sep = ''),
              '<displayfeedback feedbacktype="Solution" linkrefid="Solution"/>',
              '</respcondition>'
            )
          }
        }
      }
    }

    ## partial cloze incorrect num string answers
    if(eval$partial & x$metainfo$type == "cloze") {
      if(length(correct_answers)) {
        for(i in seq_along(correct_answers)) {
          ctype <- attr(correct_answers[[i]], "type")
          if(ctype == "string" || ctype == "num") {
            xml <- c(xml,
              '<respcondition title="Fail" continue="Yes">',
              '<conditionvar>',
              '<not>',
              correct_answers[[i]],
              '</not>',
              '</conditionvar>',
              paste('<setvar varname="SCORE" action="Add">',
                attr(correct_answers[[i]], "points")["neg"], '</setvar>', sep = ''),
              '<displayfeedback feedbacktype="Solution" linkrefid="Solution"/>',
              '</respcondition>'
            )
          }
        }
      }
    }

    ## scoring/solution display for the correct answers
    xml <- c(xml,
      '<respcondition title="Mastery" continue="Yes">',
      '<conditionvar>',
      if(!is.null(correct_answers) & (length(correct_answers) > 1 | grepl("choice", x$metainfo$type))) '<and>' else NULL
    )

    xml <- c(xml,
      unlist(correct_answers),
      if(!is.null(correct_answers) & (length(correct_answers) > 1 | grepl("choice", x$metainfo$type))) '</and>' else NULL,
      if(!is.null(wrong_answers)) {
        c('<not>', '<or>', unlist(wrong_answers), '</or>', '</not>')
      } else {
        NULL
      },
      '</conditionvar>',
      if(!eval$partial) {
        paste('<setvar varname="SCORE" action="Set">', points, '</setvar>', sep = '')
      } else NULL,
      '<displayfeedback feedbacktype="Response" linkrefid="Mastery"/>',
      '</respcondition>'
    )

    ## force display of correct answers of num exercises
    if(length(correct_num)) {
      for(j in correct_num) {
        xml <- c(xml,
          '<respcondition continue="Yes" title="Mastery">',
          '<conditionvar>',
          j,
          '</conditionvar>',
          paste('<setvar varname="SCORE" action="Add">', 0.001, '</setvar>', sep = ''),
          paste('<setvar varname="SCORE" action="Add">', -0.001, '</setvar>', sep = ''),
          '</respcondition>'
        )
      }
    }

    ## force display of all other correct answers
    if(length(correct_answers)) {
      for(j in seq_along(correct_answers)) {
        if(attr(correct_answers[[j]], "type") != "num") {
          xml <- c(xml,
            '<respcondition continue="Yes" title="Mastery">',
            '<conditionvar>',
            correct_answers[[j]],
            '</conditionvar>',
            paste('<setvar varname="SCORE" action="Add">', 0.001, '</setvar>', sep = ''),
            paste('<setvar varname="SCORE" action="Add">', -0.001, '</setvar>', sep = ''),
            '</respcondition>'
          )
        }
      }
    }


    ## handling incorrect answers
    correct_answers <- unlist(correct_answers)
    wrong_answers <- c(unlist(wrong_answers), unlist(wrong_num))

    xml <- c(xml,
      '<respcondition title="Fail" continue="Yes">',
      '<conditionvar>',
      if(!is.null(wrong_answers)) NULL else '<not>',
      if(is.null(wrong_answers)) {
        c(if(length(correct_answers) > 1) '<and>' else NULL,
          correct_answers,
          if(length(correct_answers) > 1) '</and>' else NULL)
      } else {
        c('<or>', wrong_answers, '</or>')
      },
      if(!is.null(wrong_answers)) NULL else '</not>',
      '</conditionvar>',
      if(!eval$partial) {
        paste('<setvar varname="SCORE" action="Set">', minvalue, '</setvar>', sep = '')
      } else NULL,
      '<displayfeedback feedbacktype="Solution" linkrefid="Solution"/>',
      '</respcondition>'
    )


    ## handle all other cases
    xml <- c(xml,
      '<respcondition title="Fail" continue="Yes">',
      '<conditionvar>',
      '<other/>',
      '</conditionvar>',
      paste('<setvar varname="SCORE" action="Set">', if(!eval$partial) minvalue else 0, '</setvar>', sep = ''),
      '<displayfeedback feedbacktype="Solution" linkrefid="Solution"/>',
      '</respcondition>'
    )

    ## handle unanswered cases
#    xml <- c(xml,
#      '<respcondition title="Fail" continue="Yes">',
#      '<conditionvar>',
#      '<unanswered/>',
#      '</conditionvar>',
#      '<setvar varname="SCORE" action="Set">0</setvar>',
#      '</respcondition>'
#    )

    ## end of response processing
    xml <- c(xml, '</resprocessing>')

    attr(xml, "enumerate") <- enumerate

    xml
  }
}


## function to create identfier ids
make_id <- function(size, n = 1L) {
  if(is.null(n)) n <- 1L
  rval <- matrix(sample(0:9, size * n, replace = TRUE), ncol = n, nrow = size)
  rval[1L, ] <- pmax(1L, rval[1L, ])
  colSums(rval * 10^((size - 1L):0L))
}


## delete NULL list elements
delete.NULLs <- function(x.list) {
  rval <- x.list[unlist(lapply(x.list, length) != 0)]
  rval <- if(length(rval)) rval else NULL
  rval
}


## get OLAT test results
read_olat_results <- function(file, xexam = NULL)
{
  ## checking
  stopifnot(file.exists(file <- path.expand(file)))

  ## read xexam (if any)
  if(!is.null(xexam)) {
    if(is.character(xexam)) xexam <- readRDS(xexam)
  }
 
  ## read data
  x <- readLines(file, warn = FALSE)
  x <- read.table(file, header = TRUE, sep = "\t",
    colClasses = "character", skip = 1, fill = TRUE,
    nrows = min(which(x == "")) - 3, quote = "\"")

  ## number of columns of person info
  nc <- min(grep("X1_", names(x), fixed = TRUE)) - 1

  ## only test results
  y <- x[, -(1:nc)]

  ## columns pertaining to items
  iid <- na.omit(as.numeric(unlist(sapply(
    strsplit(substr(names(y), 2, nchar(names(y))),
    "_", fixed = TRUE), head, 1))))

  ## item-wise data.frame
  y <- lapply(split(x = seq_along(iid), f = iid),
    function(ind) y[, ind, drop = FALSE])

  ## logical item x person matrix
  ipmat <- t(sapply(y, function(d) nchar(as.character(d[, ncol(d) - 1L])) > 10L))

  ## which person solved which items

  ## assume xexams object
  ## number of sections and items
#  ns <- length(xexam)
#  ni <- unique(sapply(xexam, length))
#  stopifnot(length(ni) == 1L)
  ix1 <- lapply(1:ncol(ipmat), function(i) which(as.vector(ipmat[,i])))
  ni <- max(unlist(lapply(ix1, length)))
  ns <- length(unique(iid)) / ni

  stopifnot(ns %% 1 == 0)

  ix2 <- lapply(ix1, function(i) {
    ix <- rep(NA, ni)
    ix[1L + (i - 1L) %/% ns] <- 1L + (i - 1L) %% ns
    ix
  })

  ## compute results
  process_item_result <- function(j)
  {
    rval <- lapply(1:length(ix2[[j]]), function(i) {
      id <- ix2[[j]][i]
      if(!is.na(id)) {
        ir <- y[[id + (i - 1) * ns]][j, ]
        k <- ncol(ir)
        points <- as.numeric(ir[, k - 2])
        points <- if(is.na(points)) 0 else points
        start <- ir[, k - 1]
        dur <- ir[, k]
        ssol <- ssol0 <- ir[, 1:(k - 3)]
        if(length(ssol) > 1) {
          ssol <- ssol0 <- try(gsub(".", "0", paste(ssol[-length(ssol)], collapse = ""), fixed = TRUE), silent = TRUE)
        } else {
          ssol <- try(gsub(",", ".", ssol, fixed = TRUE), silent = TRUE)
        }
        if(inherits(ssol, "try-error")) ssol <- ssol0 <- NA
        solx <- scheck <- NA
        if(!is.null(xexam) & !is.na(id)) {
          solx  <- xexam[[id]][[i]]$metainfo$solution
          tolx  <- xexam[[id]][[i]]$metainfo$tolerance
	        typex <- xexam[[id]][[i]]$metainfo$type[1]
	        ptsx  <- xexam[[id]][[i]]$metainfo$points
	        if(is.null(ptsx)) ptsx <- 1
          if(typex %in% c("mchoice", "schoice")) {
            solx <- exams::mchoice2string(solx)
            scheck <- (ssol == solx) * ptsx
          } else if(typex == "num") {
	          ssol <- as.numeric(ssol)
            scheck <- ((ssol >= solx - tolx[1L]) & (ssol <= solx + tolx[2L])) * ptsx
          } else if(typex == "cloze") {
            stop("OLAT cloze reader not yet implemented")
	        }
        }
        if(is.na(scheck)) scheck <- 0
	toPOSIXct <- function(x) ifelse(is.na(x) | x == "", NA, as.POSIXct(strptime(start, format = "%Y-%m-%dT%H:%M:%S")))
        res <- data.frame(id + (i - 1) * ns, as.numeric(points), scheck, ssol0, solx,
          toPOSIXct(start), as.numeric(dur), stringsAsFactors = FALSE)
      } else res <- data.frame(t(rep(NA, 7L)))
      res[res == ""] <- NA
      names(res) <- paste(c("id", "points", "check", "answer", "solution", "start", "duration"), i, sep = ".")
      return(res)
    })
    return(data.frame(rval, stringsAsFactors = FALSE))
  }

  res <- lapply(1:length(ix2), process_item_result)
  rval <- res[[1]]
  for(j in 2:length(res)) rval <- rbind(rval, res[[j]])
  res <- cbind(x[, 2:nc], rval)
  names(res) <- gsub("Institutionsnummer", "MatrNr", names(res))
  names(res) <- gsub("..s.", "", names(res), fixed = TRUE)

  true_false <- apply(res, 2, function(x) {
    if(is.character(x)) {
      x == "true" || x == "false"
    } else FALSE
  })
  if(any(true_false)) {
    for(i in which(true_false)) {
      wh <- FALSE
      if(!any(grepl("false", res[[i]]))) {
        res[[i]] <- rep(TRUE, length(res[[i]]))
        wh <- TRUE
      }
      if(!any(grepl("true", res[[i]]))) {
        wh <- TRUE
        res[[i]] <- rep(FALSE, length(res[[i]]))
      }
      if(!wh) res[[i]] <- res[[i]] == "true"
    }
  }

  res
}


## README: commented for now
## ## other functions, not in use yet
## ## functions to input test and item controls text
## controllist <- function(...) structure(list(...), class = "controllist")
## 
## as.character.controllist <- function(x, ...)
## {
##   paste(
##     names(x),
##     " = \"",
##     sapply(x, paste, collapse = " "),
##     "\"", sep = "", collapse = " "
##   )
## }
