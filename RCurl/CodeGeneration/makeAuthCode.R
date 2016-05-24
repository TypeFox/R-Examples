library(RGCCTranslationUnit)
curl.h = "/usr/local/include/curl/curl.h"
cpp = getCppDefines(curl.h)
defs = RGCCTranslationUnit:::processDefines(cpp, tu = tu, filter = NULL)
i = grep("CURLAUTH", names(defs$macros))
auth = defs$macros[i] 
cat("#include <curl/curl.h>",
      "#include <Rdefines.h>",
      "SEXP R_getAuthValues()", "{", 
      sprintf("\tSEXP ans = NEW_NUMERIC(%d);", length(auth)),
      "\tint i = 0;",
      sprintf("\tREAL(ans)[i++] = %s;", names(auth)), 
      "\treturn(ans);",
      "}", sep = '\n', file = "auth.c")
system("make auth.so")
dyn.load("auth.so")
vals = structure(.Call("R_getAuthValues"), names = names(auth))

con = textConnection("foo", 'w', local = TRUE)
dput(vals, con)
def = textConnectionValue(con)
close(con)


# Add this to 'setClass("CURLAuth", contains = "BitwiseValue")',
#  bitClasses.R.
# We get a warning
#   class "CURLAuth" is inheriting an inconsistent superclass structure from class "BitwiseValue", inconsistent with "SymbolicConstant"
#  if we define the class CURLAuth after we source() xbits.R
code =
  c(# 'setClass("CURLAuth", contains = "BitwiseValue")',
    paste('CURLAUTHValues =', paste(def, collapse = " ")),
    mapply(function(shortName, name, val)
            sprintf("%s <- %s <- BitwiseValue(%s, '%s', 'CURLAuth')",
                      shortName, name, as.character(val), name),
             gsub("^CURL", "", names(vals)), names(vals), vals)
    )

cat(code, file = "../R/curlAuthConstants.R", sep = "\n")
  


