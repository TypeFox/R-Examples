qmocMethods <- function(x) {
  metaObject <- qmetaObject(x)
  if (is.null(metaObject))
    stop("'x' should be a QObject instance or class")
  methods <- .Call("qt_qmocMethods", metaObject, PACKAGE="qtbase")
  methods[[1]] <- c("method", "signal", "slot", "constructor")[methods[[1]] + 1]
  methods <- c(list(sub("\\(.*", "", methods[[2]])), methods)
  names(methods) <- c("name", "type", "signature", "return", "nargs")
  as.data.frame(methods, stringsAsFactors=FALSE, row.names = methods$signature)
}

qsignals <- function(x) {
  m <- qmocMethods(x)
  s <- m[m$type == "signal", c("name", "signature")]
  rownames(s) <- NULL
  s
}
qslots <- function(x) {
  m <- qmocMethods(x)
  s <- m[m$type == "slot", c("name", "signature", "return")]
  rownames(s) <- NULL
  s
}

qnormalizedSignature <- function(x) {
  .Call("qt_qnormalizedSignature", as.character(x), PACKAGE="qtbase")
}

qresolveSignature <- function(x, sig, type, nargs) {
  ## The user can specify a signal by name or by signature
  methods <- qmocMethods(x)
  if (!missing(type))
    methods <- subset(methods, type == type)
  sigs <- as.character(methods$signature)
  if (!length(grep("\\(", sig))) {
    sigs <- sigs[grep(paste("^", sig, "\\(", sep = ""), sigs)]
    if (length(sigs) > 1) {
      argmatch <- (if (missing(nargs)) 0 else nargs) == methods[sigs, "nargs"]
      if (sum(argmatch) == 0) # just take first, may fail below
        argmatch[1] <- TRUE
      if (sum(argmatch) > 1)
        stop("ambiguous method selection: ", paste(sigs, collapse=", "))
      sigs <- sigs[argmatch]
    }
    if (!length(sigs))
      stop("method does not exist")
    sig <- sigs[1]
  } else {
    sig <- qnormalizedSignature(sig)
    if (!(sig %in% sigs))
      stop("method does not exist")
  }
  ## slots cannot have default arguments, usually signals are multiplied
  if (!missing(nargs) && methods[sig, "nargs"] > nargs)
    stop("method has more arguments than provided signature")
  sig
}

qproperties <- function(x) {
  metaObject <- qmetaObject(x)
  if (!is.null(metaObject)) {
    props <- .Call("qt_qproperties", metaObject, PACKAGE="qtbase")
    names(props) <- c("name", "type", "readable", "writable")
    name <- props$name
    props$name <- NULL
    as.data.frame(props, row.names=name)
  } else data.frame(name = character(), type = character(),
                    readable = logical(), writable = logical())
}

qmetadata <- function(x) {
  attr(x, "metadata")
}

## Every time metadata is set, we set the accessor methods to
## regenerate it. After regneration, the accessors merely retrieve the
## compiled form.
"qmetadata<-" <- function(x, value) {
  if (!isQObjectClass(x)) # not a QObject class, no compilation
    return(x)
  
  ### We use direct substitution here to avoid any conflicts with the
  ### instance environment. These methods are set on any user-defined,
  ### QObject-drived class, so we need to be careful. We could go up
  ### the environment chain, but this approach seems more robust.
  
  ### FIXME: it is possible for this compiled metadata to end up in a
  ### lazily-loaded package. How could we avoid this, short of having
  ### the user package reset the metadata for each class on load?
  createMetaObject <- eval(substitute(function() {
    compiledMetaObject <- qtbase:::compileMetaObject(class, metadata)
    getMetaObject <- eval(substitute(function() { x },
                                     list(x = compiledMetaObject)))
    qtbase::qsetMethod("metaObject", class, getMetaObject)
    qtbase::qsetMethod("staticMetaObject", class, getMetaObject)
    compiledMetaObject
  }, list(class = x, metadata = value)))
  qsetMethod("metaObject", x, createMetaObject)
  qsetMethod("staticMetaObject", x, createMetaObject)
  qsetMethod("qt_metacall", x, # just ensure we have this method defined
             function(call, id, args)
             .Call("qt_qmetacall", this, call, id, args, PACKAGE="qtbase"))
  x
}

isQObjectClass <- function(x) !is.null(x$staticMetaObject)

qmetaObject <- function(x, ...) UseMethod("qmetaObject")
qmetaObject.default <- function(x, ...) NULL
qmetaObject.QMetaObject <- function(x) x
qmetaObject.QObject <- function(x) x$metaObject()
qmetaObject.RQtClass <- function(x) {
  if (isQObjectClass(x))
    x$staticMetaObject()
  else NULL
}
qmetaObject.RQtUserClass <- function(x) {
### FIXME: hack because we do not support user static methods
  staticMetaObject <- attr(x, "instanceEnv")$staticMetaObject
  if (!is.null(staticMetaObject))
    staticMetaObject()
  else NextMethod()
}

## Does not handle 'const' types yet -- but do we want this?
qmetaMethod <- function(signature, access = c("public", "protected", "private"),
                        argNames)
{
  access <- match.arg(access)
  if (!grepl(")$", signature))
    signature <- paste(signature, "()", sep = "")
  signature <- qnormalizedSignature(signature)
  args <- strsplit(gsub(".*?\\((.*?)\\)", "\\1", signature), ",",
                   fixed=TRUE)[[1]]
  haveNames <- grepl("[^ ] [^ ]", args)
  if (!all(haveNames)) {
    if (any(haveNames))
      stop("If any arguments are named, all must be named")
    if (missing(argNames))
      argNames <- paste("x", seq(args), sep = "")
    argTypes <- args
  } else {
    argTokens <-
      matrix(as.character(unlist(strsplit(args, " ", fixed=TRUE))), 2)
    argTypes <- argTokens[1,]
    argNames <- argTokens[2,]
  }
  lhs <- strsplit(sub("\\(.*", "", signature), " ", fixed=TRUE)[[1]]
  name <- tail(lhs, 1)
  returnType <- paste(head(lhs, -1), collapse = " ")
  if (!length(returnType))
    returnType <- ""
  signature <- paste(name, "(", paste(argTypes, collapse=","), ")", sep="")
  list(signature = signature, args = argNames, type = returnType,
       access = access, name = name)
}

### Some stuff derived from QtRuby for creating a MetaData blob

## 'x' is a class
## 'metadata' is an environment of classinfos, signals and slots

### It lets an R class define:
### - Slots: useful for providing dbus services, and easy
### - Signals: very useful
### - Class info: useful for describing dbus services

### TODO:
### - Properties: Encapsulated fields are nice, but do we want to be bound
###   by QValue? Could just implement our own custom properties.

compileMetaObject <- function(x, metadata) {
  metalist <- as.list(metadata)
  
  infos <- metadata$classinfos
  signals <- metadata$signals
  slots <- metadata$slots
  props <- as.list(metadata$properties)
  ## we only export to metadata if a type is provided
  props <- Filter(function(p) !is.null(p$type), props)
  metalist$properties <- lapply(props, `[`, c("name", "type"))
  
  ## generate 'stringdata' table
  allNames <- unique(c(attr(x, "name"), unlist(metalist), ""))
  offsets <- cumsum(c(0, head(nchar(allNames), -1) + 1))
  names(offsets) <- allNames
  
  stringdata <- charToRaw(paste(allNames, collapse=";"))
  stringdata[stringdata == charToRaw(";")] <- as.raw(0)

  ## generate 'data' table

  ##
  ## From the enum MethodFlags in qt-copy/src/tools/moc/generator.cpp
  ##
  AccessPrivate <- 0x00
  AccessProtected <- 0x01
  AccessPublic <- 0x02
  MethodMethod <- 0x00
  MethodSignal <- 0x04
  MethodSlot <- 0x08
  MethodCompatibility <- 0x10
  MethodCloned <- 0x20
  MethodScriptable <- 0x40
  
  ##
  ## From the enum PropertyFlags in qmetaobject_p.h
  ##
  PropertyFlags <- list(Readable = 0x1,
                        Writable = 0x2,
                        Resettable = 0x4,
                        EnumOrFlag = 0x8,
                        Constant = 0x400,
                        Final = 0x800,
                        Designable = 0x1000,
                        ResolveDesignable = 0x2000,
                        Scriptable = 0x4000,
                        ResolveScriptable = 0x8000,
                        Stored = 0x10000,
                        ResolveStored = 0x20000,
                        Editable = 0x40000,
                        ResolvedEditable = 0x80000,
                        User = 0x100000,
                        ResolveUser = 0x200000,
                        Notify = 0x400000)

  access <- c(private = AccessPrivate, protected = AccessProtected,
              public = AccessPublic)

  nmethods <- length(signals) + length(slots)
  ninfos <- length(infos)
  nprops <- length(props)
  headerLen <- 10L
  data <-
    c(1, # revision
      offsets[attr(x, "name")], 	# classname
      length(infos), headerLen, # class info
      nmethods, headerLen + 2*ninfos, # methods
      nprops, headerLen + 2*ninfos + 5*nmethods, # properties
      0, 0) # enums

  ## the class info
  data <- c(data, offsets[rbind(names(infos), infos)])

  methodData <- function(methods, flag) {
    do.call("c", lapply(methods, function(method) {
      ## cannot use character extraction due to empty strings
      if (!length(method$args))
        method$args <- ""
      offm <- match(c(method$signature, method$args, method$type, ""),
                    names(offsets))
      c(offsets[offm], access[method$access] + flag + MethodScriptable)
    }))
  }
  
  ## the signals
  data <- c(data, methodData(signals, MethodSignal))
  
  ## the slots
  data <- c(data, methodData(slots, MethodSlot))

  ## the properties
  data <- c(data, do.call(c, lapply(props, function(prop) {
    flags <- with(PropertyFlags, {
      flags <- Designable + Scriptable + Readable
      if (!is.null(prop$write))
        flags <- flags + Writable
      if (is.character(prop$notify))
        flags <- flags + Notify
      ## if (!is.null(prop$reset))
      ##   flags <- flags + Resettable
      if (isTRUE(props$constant))
        flags <- flags + Constant
      if (isTRUE(props$final))
        flags <- flags + Final
      ## if (isTRUE(props$designable))
      ##   flags <- flags + Designable
      ## if (isTRUE(props$scriptable))
      ##   flags <- flags + Scriptable
      if (isTRUE(props$stored))
        flags <- flags + Stored
      if (isTRUE(props$user))
        flags <- flags + User
      flags
    })
    c(offsets[c(prop$name, prop$type)], flags)
  })))

  notifies <- sapply(props, function(p) {
    if (is.null(p$notify))
      NA
    else p$notify
  })
  signalNames <- sapply(signals, `[[`, "signature")
  notify_ids <- match(notifies, signalNames) - 1L 
  notify_ids[is.na(notify_ids)] <- 0
  data <- c(data, notify_ids)
  
  .Call("qt_qnewMetaObject", x, stringdata, data, PACKAGE="qtbase")
}
