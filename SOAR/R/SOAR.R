.enc  <- local({

  bad <- c(" ", "<", ">", ":", "\"", "/", "\\", "|", "?", "*")
  rpl <- paste("@", 0:9, sep = "")

  regex <- paste("(", paste(LETTERS, collapse = "|"), ")", sep = "")

  function (x) {
    x <- gsub("@", "\t", x, fixed = TRUE)
    for (i in seq(along = bad))
      x <- gsub(bad[i], rpl[i], x, fixed = TRUE)
    x <- gsub(regex, "@\\1", x)
    x <- gsub("\t", "@@", x, fixed = TRUE)
    paste(x, "@.RData", sep = "")
  }
})

.dec  <- local({

  bad <- c(" ", "<", ">", ":", "\"", "/", "\\", "|", "?", "*")
  rpl <- paste("@", 0:9, sep = "")

  regex <- paste("@(", paste(LETTERS, collapse = "|"), ")", sep = "")

  function (x) {
    x <- gsub("@@", "\t", x, fixed = TRUE)
    x <- sub("@\\.RData$", "", x)
    x <- sub("\\.RData$", "", x)
    x <- gsub(regex, "\\1", x)
    for (i in seq(along = bad))
      x <- gsub(rpl[i], bad[i], x, fixed = TRUE)
    gsub("\t", "@", x, fixed = TRUE)
  }
})

.mostFiles <- function(path)
  setdiff(dir(path, all.files = TRUE), c(".", ".."))

.pathAttributes <- function () {
  s <- search()
  paths <- lapply(1:length(s),
                  function(i) attr(as.environment(i), "path"))
  paths[[length(s)]] <- system.file()
  m <- sapply(paths, is.null)
  paths[m] <- s[m]
  unlist(paths)
}

.attach <- function(directory, pos = 2,
                    warn = !file.exists(directory),
                    readonly) {
  env <- attach(NULL, pos, basename(directory))
  attr(env, "path") <- directory
  attr(env, "readonly") <- readonly
  if (file.exists(directory)) {
    fils <- .mostFiles(directory)
    objs <- .dec(fils)
    fils <- file.path(directory, fils)
    for(i in seq(along = objs))
      eval(substitute(delayedAssign(OBJECT, {
        load(file = FILE)
        get(OBJECT)
      }), list(OBJECT = objs[i], FILE = fils[i])),
           envir = env)
  } else if (warn)
    warning(paste(directory,
                  "does not currently exist. ",
                  "A call to 'Store' will create it."),
            call. = FALSE)
}

.makeClone <- function(Name, Which) {
  dsn <- deparse(substitute(Name))
  dsw <- deparse(substitute(Which))
  f <- function(...) {}
  body(f) <- substitute({
    Call <- match.call()
    Call[[1]] <- quote(SOAR::NAME)
    if(is.null(Call[["lib"]]))
      Call[["lib"]] <- Sys.getenv(LIB, unset = WHICH)
    if(is.null(Call[["lib.loc"]]))
      Call[["lib.loc"]] <- Sys.getenv(LIB_LOC, unset = path.expand("~"))
    eval.parent(Call)
  }, list(NAME = as.name(dsn),
          WHICH = paste(".R", dsw, sep = "_"),
          LIB = paste("R_CENTRAL", toupper(dsw), sep = "_"),
          LIB_LOC = "R_CENTRAL_LIB_LOC"))
  f
}

NAME <- function() NULL

Attach <- function(lib = Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache"),
                   lib.loc = Sys.getenv("R_LOCAL_LIB_LOC", unset = "."),
                   pos = 2, uniquely = TRUE, readonly = FALSE, ...) {
  if(class((x <- substitute(lib))) == "name")
    lib <- deparse(x) else lib <- lib
  if(!(file.exists(lib.loc) && file.info(lib.loc)$isdir))
    stop(lib.loc, " is not an existing directory.", call. = FALSE, domain = NA)
  path <- file.path(lib.loc, lib)
  if(file.exists(path)) {
    if(!file.info(path)$isdir)
      stop(path, " exists but is not a directory!", call. = FALSE, domain = NA)
    fils <- .mostFiles(path)
    if(any(i <- grep("@\\.RData$", fils, invert = TRUE))) {
      warning(paste("Converting", path, "from old filename format to new"),
              call. = FALSE)
      newF <- .enc(.dec(fils[i]))
      for(j in i) file.rename(file.path(path, fils[j]),
                              file.path(path, newF[j]))    }
  }
  if(uniquely)
    if(any(m <- (.pathAttributes() == path)))
      for(j in rev(sort(which(m)))) detach(pos = j)
  .attach(directory = path, pos = pos, readonly = readonly, ...)
}

AttachData <- .makeClone(Attach, Data)
AttachUtils <- .makeClone(Attach, Utils)

Store <- function(..., list = character(0),
                  lib = Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache"),
                  lib.loc = Sys.getenv("R_LOCAL_LIB_LOC", unset = "."),
                  remove = TRUE) {
  if(class((.x <- substitute(lib))) == "name")
    lib <- deparse(.x) else lib <- lib
  if(!(file.exists(lib.loc) && file.info(lib.loc)$isdir))
    stop(lib.loc, " is not an existing directory.", call. = FALSE, domain = NA)
  path <- file.path(lib.loc, lib)
  if(file.exists(path) && !file.info(path)$isdir)
      stop(path, " exists but is not a directory!", call. = FALSE, domain = NA)
  if(m <- match(path, .pathAttributes(), nomatch = FALSE)) {
    e <- as.environment(m[1])
    if(!is.null(n <- attr(e, "readonly")) && n)
      stop(path, " is attached as read only!", call. = FALSE, domain = NA)
  }
  nam <- list
  if(!is.null(m <- match.call(expand.dots = FALSE)$...))
    nam <- c(nam,
             sapply(m, function(x)
                    switch(class(x),
                           name = deparse(x),
                           call = {
                             o <- eval(x, envir = parent.frame(n = 4))
                             if(!is.character(o))
                               stop("non-character name!",
                                    call. = FALSE, domain = NA)
                             o
                           },
                           character = x,
                           stop("garbled call to 'Store'",
                                call. = FALSE, domain = NA))))
  if(length(nam) == 0) return()
  nam <- drop(sort(unique(nam)))
  if(file.exists(path)) {
    if(!file.info(path)$isdir)
      stop(path, " exists but is not a directory!", call. = FALSE, domain = NA)
    fils <- .mostFiles(path)
    if(any(i <- grep("@\\.RData$", fils, invert = TRUE))) {
      warning(paste("Converting", path, "from old filename format to new"),
              call. = FALSE, domain = NA)
      newF <- .enc(.dec(fils[i]))
      for(j in i) file.rename(file.path(path, fils[j]),
                              file.path(path, newF[j]))
    }
  } else dir.create(path)
  for(n in nam) {
    no <- !exists(n, inherits = FALSE, envir = parent.frame())
    comm <- if (no) substitute({
      assign(N, get(N))
      save(list = N, file = FILE)
      rm(list = N)
    }, list(N = n, FILE = file.path(path, .enc(n)))) else
    substitute({
      save(list = N, file = FILE)
    }, list(N = n, FILE = file.path(path, .enc(n))))
    eval.parent(comm)
  }
  pos <- if(any(m <- (.pathAttributes() == path))) {
    m <- which(m)[1]
    detach(pos = m)
    m
  } else 2
  .attach(path, pos = pos, warn = FALSE, readonly = FALSE)
  if(remove) {
    o <- intersect(eval.parent(quote(objects(all.names = TRUE))), nam)
    if(length(o) > 0)
      eval.parent(substitute(remove(list = O), list(O = o)))
  }
}

StoreData <- .makeClone(Store, Data)
StoreUtils <- .makeClone(Store, Utils)

Objects <- function(lib = Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache"),
                    lib.loc = Sys.getenv("R_LOCAL_LIB_LOC", unset = "."),
                    all.names = FALSE, pattern = ".*", readonly = FALSE) {
  if(class((x <- substitute(lib))) == "name")
    lib <- deparse(x) else lib <- lib
  if(!(file.exists(lib.loc) && file.info(lib.loc)$isdir))
    stop(lib.loc, " is not an existing directory.", call. = FALSE, domain = NA)
  path <- file.path(lib.loc, lib)
  if(file.exists(path) && !file.info(path)$isdir)
      stop(path, " exists but is not a directory!", call. = FALSE, domain = NA)
  if(!is.na(pos <- match(path, .pathAttributes())))
    objects(pos = pos, all.names = all.names, pattern = pattern) else {
      .attach(path, pos = 2, warn = FALSE, readonly = readonly)
      objects(2, all.names = all.names, pattern = pattern)
    }
}

ObjectsData <- .makeClone(Objects, Data)
ObjectsUtils <- .makeClone(Objects, Utils)
Ls <- Objects
LsUtils <- .makeClone(Ls, Utils)
LsData <- .makeClone(Ls, Data)


Remove <- function(..., list = character(0),
                   lib = Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache"),
                   lib.loc = Sys.getenv("R_LOCAL_LIB_LOC", unset = ".")) {
  if(class((.x <- substitute(lib))) == "name")
    lib <- deparse(.x) else lib <- lib
  if(!(file.exists(lib.loc) && file.info(lib.loc)$isdir))
    stop(lib.loc, " is not an existing directory.", call. = FALSE, domain = NA)
  path <- file.path(lib.loc, lib)
  if(file.exists(path) && !file.info(path)$isdir)
      stop(path, " exists but is not a directory!", call. = FALSE, domain = NA)
  if(m <- match(path, .pathAttributes(), nomatch = FALSE)) {
    e <- as.environment(m[1])
    if(!is.null(n <- attr(e, "readonly")) && n)
      stop(path, " is attached as read only!", call. = FALSE, domain = NA)
  }
  nam <- list
  if(!is.null(m <- match.call(expand.dots = FALSE)$...))
    nam <- c(nam,
             sapply(m, function(x)
                    switch(class(x),
                           name = deparse(x),
                           call = {
                             o <- eval(x, envir = parent.frame(n = 4))
                             if(!is.character(o))
                               stop("non-character name!",
                                    call. = FALSE, domain = NA)
                             o
                           },
                           character = x,
                           stop("garbled call to 'Remove'"))))
  if(length(nam) == 0) return()
  nam <- .enc(sort(unique(nam)))
  fils <- .mostFiles(path)
  if(any(i <- grep("@\\.RData$", fils, invert = TRUE))) {
    warning(paste("Converting", path, "from old filename format to new"),
            call. = FALSE, domain = NA)
    newF <- .enc(.dec(fils[i]))
    for(j in i) file.rename(file.path(path, fils[j]),
                            file.path(path, newF[j]))
    fils <- .mostFiles(path)
  }
  mnam <- intersect(nam, fils)
  if(length(d <- setdiff(nam, mnam)) > 0)
    warning(gettextf("Object %s not found and hence not removed\n",
                     format(.dec(d))), call. = FALSE, domain = NA)
  file.remove(file.path(path, mnam))
  pos <- if(any(m <- (.pathAttributes() == path))) {
    m <- which(m)[1]
    detach(pos = m)
    m
  } else 2
  .attach(path, pos = pos, readonly = FALSE)
}

RemoveData <- .makeClone(Remove, Data)
RemoveUtils <-.makeClone(Remove, Utils)

Search <- local({
  .zf <- function (s)
      paste(substring(paste(rep(0, (m <-
                                    max(n <-
                                        nchar(s <-
                                              as.character(s))))),
                            collapse = ""), 0, m - n), s, sep = "")

  function (abbrev = FALSE) {
    d <- search()
    wd <- normalizePath(getwd(), winslash = "/")
    f <- ""
    e <- globalenv()
    for (j in 2:length(d)) {
      e <- parent.env(e)
      p <- attr(e, "path")
      if (is.null(p))
          p <- ""
      else p <- dirname(p)
      f <- c(f, p)
    }
    f[length(f)] <- normalizePath(dirname(system.file()), winslash="/")
    rhome <- normalizePath(R.home(), winslash = "/")
    home <- normalizePath("~", winslash = "/")
    f <- gsub(wd, ".",
              gsub(rhome, "R_HOME",
                   gsub(home, "~",
                        f, fixed = TRUE), fixed = TRUE), fixed = TRUE)
    abbrev <- abbrev[1]
    if (is.logical(abbrev) && abbrev)
        abbrev <- 50
    else abbrev <- as.numeric(abbrev)
    if (!is.na(abbrev) && abbrev > 0)
        if (any(w <- ((n <- nchar(f)) > abbrev))) {
          f[w] <- paste("...", substring(f[w], n[w] - abbrev +
                                         1, n[w]))
        }
    d <- cbind(name = d, lib.loc = f)
    dimnames(d)[[1]] <- .zf(1:nrow(d))
    noquote(d)
  }
})
