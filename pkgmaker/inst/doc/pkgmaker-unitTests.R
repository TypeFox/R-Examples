
## ----setup, include=FALSE------------------------------------------------
pkg <- 'pkgmaker'
require( pkg, character.only=TRUE )
prettyVersion <- packageDescription(pkg)$Version
prettyDate <- format(Sys.Date(), '%B %e, %Y')
authors <- packageDescription(pkg)$Author


