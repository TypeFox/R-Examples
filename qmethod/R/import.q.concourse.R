import.q.concourse <- function(q.concourse.dir, languages=NULL) {

  # Input validation ===============
  # q.concourse.dir cannot be validated, there appears to be no file.exists equivalent for dirs
  # languages (if applicable) are validated further down


  # Find all meaningful item short handles across all languages ================
  q.item.handles <- c() #  set up empty vector

  if (is.null(languages)) {  # monolingual version
    q.item.handles <- list.files(
      path = q.concourse.dir,
      no.. = TRUE,  # not dotfiles
      pattern = "\\.tex$"
    )
    q.item.handles <- file_path_sans_ext(q.item.handles)
  } else {  # multilingual version
    for (lang in languages) {
      lang <- list.files(
        path = paste(q.concourse.dir, lang, "/", sep=""),
        no.. = TRUE, #  no dotfiles
        pattern = "\\.tex$"  # only tex
      )
      lang <- file_path_sans_ext(lang) #  kill extensions
      q.item.handles <- append(q.item.handles, lang) # append vector
    }
  }
  q.item.handles <- unique(q.item.handles) #  kill duplicates

  # Setting up empty matrix ====================================================
  q.concourse <- matrix(#  items and languages are 2 dimensions
    data = NA, #  no such thing yet, so all NAs (that should be the baseline!)
    byrow = TRUE, #  items are filled per row
    nrow = length(q.item.handles),
    ncol = length(languages),
    dimnames = list( #  dims should be called accordingly ...
      q.item.handles, #  rows are named by handle
      languages #  columns by language
    )
  )

  # Reading in full items ======================================================
  if (is.null(languages)) {  # monolingual version
    for (handle in q.item.handles) {  # loop over item handles
      path <- paste(q.concourse.dir, handle, ".tex", sep = "")  # establish path
      q.concourse[handle] <- readChar(path, file.info(path)$size) # assign the full text
    }
  } else {  # multilingual version
    for (lang in languages) {  # loop over the languages
      for (handle in q.item.handles) {  # loop over item handles
        path <- paste(q.concourse.dir, lang, "/", handle, ".tex", sep = "")  # establish path
        if (file.exists(path)) {  # if translation exists, read in
          q.concourse[handle,lang] <- readChar(path, file.info(path)$size) # assign the full text
        } else {  # error out if there is a missing translation
          stop(paste("There is no", lang, "version for", handle, "."))
        }
      }
    }
  }
  if (is.null(languages)) {  # monolingual version
    q.concourse <- as.matrix(q.concourse)  # other functions expect matrix
    message(
      paste(
        "Gathered",
        length(q.concourse),
        "full items."
      )
    )
  } else {  # multilingual version
    message(
      paste(
        "Gathered",
        length(q.concourse)/length(languages),
        "full items, each in",
        length(languages),
        "languages."
      )
    )
  }
  return(q.concourse)
}