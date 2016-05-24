cat("*** R VERSION\n")
cat(".R.version:\n")
print(R.version)

cat("\n\n*** PLATFORM\n")
cat(".Platform:\n")
print(.Platform)

cat("\n\n*** MACHINE\n")
cat(".Machine:\n")
print(.Machine)

cat("\n\n*** DIRECTORIES\n")
cat("R.home() directory: ", R.home(), "\n", sep="")
for (dir in c("bin", "doc", "etc", "modules", "share")) {
  cat(sprintf("R.home('%s') directory: %s\n", dir, R.home(dir)))
}
cat("Working directory: ", getwd(), "\n", sep="")
cat("Temporary directory: ", tempdir(), "\n", sep="")
cat("User home directory: ", normalizePath("~"), "\n", sep="")

cat("\n\n*** LIBRARY PATHS\n")
cat(".libPaths():\n")
print(.libPaths())

cat("\n\n*** ENVIRONMENT VARIABLE\n")
cat("Sys.getenv():\n")
print(as.list(Sys.getenv()))

cat("\n\n*** LOCALE\n")
cat("Sys.getlocale():\n")
print(Sys.getlocale())

cat("\n\n*** OPTIONS\n")
cat("options():\n")
print(options())

cat("\n\n*** LOADED NAMESPACES\n")
cat("loadedNamespaces():\n")
print(loadedNamespaces())

cat("\n\n*** SEARCH PATH\n")
cat("search():\n")
print(search())

cat("\n\n*** INTERACTIVITY\n")
cat("interactive(): ", interactive(), "\n", sep="")

cat("\n\n*** CAPABILITIES\n")
cat("capabilities():\n")  ## May try to open X11
print(capabilities())  ## May try to open X11

cat("\n\n*** COMPLETION SETTINGS\n")
cat("rc.status():\n")
print(rc.status())
cat("rc.options():\n")
print(rc.options())

cat("\n\n*** SESSION INFORMATION\n")
cat("utils::sessionInfo():\n")
print(utils::sessionInfo())
