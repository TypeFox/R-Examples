.onAttach <- function(...) {
   packageStartupMessage("Welcome to archivist.github (ver: ", utils::packageVersion("archivist.github"), ").")
}

.onLoad <- function(...) {
  
  archivist::aoptions("readmeDescription", "A Repository of Artifacts supported by [archivist](https://github.com/pbiecek/archivist) \n\n
  [`Repository`](https://github.com/pbiecek/archivist/wiki/archivist-package-Repository) stores specific values of an artifact,
  different for various artifact's classes and artifacts themselves. To learn more about artifacts visit [wiki](https://github.com/pbiecek/archivist/wiki).")
  archivist::aoptions("repoDescription", "A Repository of Artifacts supported by archivist https://github.com/pbiecek/archivist")
  archivist::aoptions("response",  FALSE)
  archivist::aoptions("commitMessage", NULL)
  archivist::aoptions("alink", FALSE)
}

## no S4 methodology here; speedup :
.noGenerics <- TRUE
