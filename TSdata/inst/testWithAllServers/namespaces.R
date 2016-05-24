z <- loadedNamespaces()

unloadNamespace("tseries") # loads zoo ?

#if (!all(z == loadedNamespaces())
#    stop("unloading a not loaded namespace is changing loadedNamespaces.")

require("RMySQL")
require("RPostgreSQL")
unloadNamespace("RMySQL")
unloadNamespace("RPostgreSQL")

require("TSmisc")
##require("TSjson")
unloadNamespace("TSmisc")
##unloadNamespace("TSjson")

require("TSmisc")
require("TSSQLite")
unloadNamespace("TSmisc")
unloadNamespace("TSSQLite")

require("TSmisc")
require("TSPostgreSQL")
unloadNamespace("TSPostgreSQL")
unloadNamespace("TSmisc")

require("TSmisc")
require("TSMySQL")
unloadNamespace("TSmisc")
unloadNamespace("TSMySQL")

require("TSMySQL")
require("TSSQLite")
unloadNamespace("TSMySQL")
unloadNamespace("TSSQLite")

require("TSMySQL")
require("TSPostgreSQL")

# next previously failed with  argument "where" is missing, with no default 
#  in Jan 2015 R-devel and MySQL 0.10.2
# works April 21, 2015 with R-devel  and MySQL 0.10.3
unloadNamespace("TSMySQL")
unloadNamespace("TSPostgreSQL")

#detach("package:TScompare", unload=TRUE)
#detach("package:TSmisc", unload=TRUE)
##detach("package:TSjson", unload=TRUE)
#detach("package:WriteXLS", unload=TRUE)

loadedNamespaces()
search()

