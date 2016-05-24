.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  if (interactive()) {
    packageStartupMessage('WARNING!\nThere are some quite big changes in this version of rmongodb.
mongo.bson.to.list, mongo.bson.from.list (which are workhorses of many other rmongofb high-level functions) are rewritten.
Please, \nTEST IT BEFORE PRODUCTION USAGE.\nAlso there are some other important changes, please see NEWS file and release notes at
https://github.com/mongosoup/rmongodb/releases/')
  }
}
