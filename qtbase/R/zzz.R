
.noGenerics <- TRUE

globalVariables("this")

.onUnload <- function(libpath) {
 .Call("cleanupQtApp", PACKAGE="qtbase")
 library.dynam.unload("qtbase", libpath)
}

.onLoad <- function(libname, pkgname) 
{
  if (.Platform$OS.type=="windows") {
    dllpath <- Sys.getenv("QTBASE_QT_PATH")
    if (!nzchar(dllpath))
      dllpath <- .windows_qt_path()
    if (!file.exists(dllpath))
      .install_system_dependencies()
    library.dynam("qtbase", pkgname, libname, DLLpath = dllpath)
  } else library.dynam("qtbase", pkgname, libname)

  qlibrary(Qt, NULL)

  ### HACK: populate QGlobalSpace now, as other Smoke libs will have
  ### the same class. Really, Smoke should namespace this thing, but
  ### maybe we cannot assume that every class has a unique name?
  Qt$QGlobalSpace
  
  ## Prefer OpenGL1.x engine, rather than the OpenGL2 ES engine
  ## R is usually running on non-mobile platforms
  ## This must be called before QApplication is constructed!
  ## Currently disabled, since 1.x engine seems unmaintained.
  if (!is.null(Qt$QGL$setPreferredPaintEngine))
    Qt$QGL$setPreferredPaintEngine(Qt$QPaintEngine$OpenGL)

  locale <- saveLocale()
  
  .Call("addQtEventHandler", PACKAGE="qtbase")

  restoreLocale(locale)
  
  if (!is.null(Qt$QGL$setPreferredPaintEngine))
    Qt$QGL$setPreferredPaintEngine(Qt$QPaintEngine$OpenGL2)
  
### Temporarily disabled until we figure out why this crashes (on my machine)
  ## reg.finalizer(getNamespace("qtbase"), function(ns)
  ##               {
  ##                 if ("qtbase" %in% loadedNamespaces())
  ##                   .onUnload(file.path(libname, pkgname))
  ##               }, onexit=TRUE)
}

## Qt will muck up the current locale when it is initialized, so we
## work around that here. Ideally, we would translate the current
## locale into the default QLocale, but that does not seem possible in
## general, as a QLocale (except the one based on the internal
## QSystemLocale) only supports a single language for all categories.

saveLocale <- function() {
  categories <- c("LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_NUMERIC",
                  "LC_TIME", "LC_MESSAGES", "LC_PAPER", "LC_MEASUREMENT")
  nonEmptyString <- function(x) {
    is.character(x) && nzchar(x)
  }
  Filter(nonEmptyString, sapply(categories, Sys.getlocale, simplify = FALSE))
}

restoreLocale <- function(locale) {
  suppressWarnings(mapply(Sys.setlocale, names(locale), locale))
}
