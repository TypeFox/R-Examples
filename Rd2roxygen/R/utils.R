## extract tags
tag = function(x) attr(x, "Rd_tag")

## replace tags
untag = function(x) {
  if (is.null(x)) return(NULL)
  attr(x, "Rd_tag") = "TEXT"
  x
}

## construct strings from rd
reconstruct = function(rd) {
  if (is.null(rd)) return()

  if (is.list(rd)) {
    if (length(tag(rd)) && tag(rd) %in% c('\\item', '\\tabular', '\\eqn', '\\deqn', '\\link')) {
      if (tag(rd) == '\\link')
        return(paste('\\link', sprintf('[%s]', attr(rd, 'Rd_option')), '{', rd, '}', sep = ""))
      if (length(rd) == 2) {
        return(paste(tag(rd), '{', rd[[1]], '}{',
                     paste(sapply(rd[[2]], reconstruct), collapse = ""),
                     '}', sep = "", collapse = ""))
      } else if (length(rd) == 0) return(tag(rd))
    }
    special = tag(rd) == toupper(tag(rd))
    singles = tag(rd) %in% c('\\tab', '\\cr')
    prefix = ifelse(special, "",
                     paste(tag(rd), ifelse(singles, "", "{"), sep = ""))
    suffix = ifelse(special, "", ifelse(singles, "", "}"))
    paste(prefix, paste(sapply(rd, reconstruct), collapse = ""), suffix,
          sep = "")
  } else {
    if (tag(rd) == 'TEXT') gsub('%', '\\%', rd, fixed = TRUE) else rd
  }
}

## wrap strings with comment prefix
comment_line = function(x, exdent = 0) {
  if (missing(x)) return(comment_prefix())

  strwrap(x, width = 80, exdent = exdent, prefix = comment_prefix())
}

## add comments
comment_tag = function(tag, value) {
  if (is.null(value) || value == "" || length(value) == 0) return()

  comment_line(paste(tag, value), exdent = 0)
}

## access the comment prefix
comment_prefix = function() {
  getOption("roxygen.comment", "#' ")
}

Rbin = function() shQuote(file.path(R.home('bin'), 'R'))

tidy_examples = function(rd, idx0, idx1, ..., path) {
  tmp = rd[idx0:idx1]
  if (length(tmp) > 1 && tmp[2] == '# !formatR') {
    rd = rd[-(idx0 + 1)]  # remove this token
    return(rd)
  }
  tmp[1] = sub('^\\\\examples\\{', '', tmp[1])
  nn = length(tmp)
  tmp[nn] = sub('\\}$', '', tmp[nn])
  txt = gsub('\\%', '%', tmp, fixed = TRUE) # will escape % later
  txt = sub('^\\\\+dont(run|test|show)', 'tag_name_dont\\1 <- function() ', txt)
  txt = tidy_code(txt, ...)
  if (!inherits(txt, 'try-error')) {
    txt = gsub('(^|[^\\])%', '\\1\\\\%', txt)
    txt = gsub('tag_name_dont(run|test|show) <- function\\(\\) \\{', '\\\\dont\\1{', txt)
    txt[txt == ''] = '\n'
    txt = unlist(strsplit(txt, '\n'))
    # remove the four spaces introduced by disguising \\dontrun as a function
    if (length(idx2 <- grep('\\\\dont(run|test|show)\\{', txt))) {
      for (i in idx2) {
        j = i + 1
        while (txt[j] != '}') {
          txt[j] = sub('^    ', '', txt[j])
          j = j + 1
        }
      }
    }
    txt = paste(txt, rep(c('\n', ''), c(length(txt) - 1, 1)), sep = '', collapse = '')
    txt[1] = paste('\\examples{', txt[1], sep = '')
    nn0 = length(txt)
    txt[nn0] = paste(txt[nn0], '}', sep = '')
    rd[idx0] = paste(txt, collapse = '\n')
    if (idx1 > idx0) rd = rd[-((idx0 + 1):idx1)]
  } else {
    message('(!) failed to reformat examples code in ', path)
    message(paste('   ', tmp, collapse = '\n'))
  }
  rd
}

in_dir = function(dir, expr) {
  owd = setwd(dir); on.exit(setwd(owd))
  expr
}
