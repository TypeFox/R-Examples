files = list.files("~/Rpackages/XML/R-ex",  pattern = "[^~]$", full = TRUE)
#files = files[ - grep("xmlFlatList", files)]
# Seems to work.
# Taking dtd out of this exclusion list yielded a crash
#[ 1] "/Users/duncan/Rpackages/XML/R-ex/addNode.R"             
#[ 2] "/Users/duncan/Rpackages/XML/R-ex/dtdElement.R"          
#[ 3] "/Users/duncan/Rpackages/XML/R-ex/dtdElementValidEntry.R"
#[ 4] "/Users/duncan/Rpackages/XML/R-ex/dtdIsAttribute.R"      
#[ 5] "/Users/duncan/Rpackages/XML/R-ex/dtdValidElement.R"     
#[ 6] "/Users/duncan/Rpackages/XML/R-ex/parseDTD.R"            
#[ 7] "/Users/duncan/Rpackages/XML/R-ex/xmlContainsEntity.R"   
#[ 8] "/Users/duncan/Rpackages/XML/R-ex/xmlStopParser.R"

#
# Put parseDTD.R back in seems okay.
#
# addNode.R - fails.
# But run with gctorture, succeeds the first time and then on the second run
#  get a no applicable method for xmlChildren as tt$top returns NULL.

# Add all the ^dtd.* back
#    works if xmlFlatList is not there.
# grep("ContainsEntity|xmlStopParser.R|addNode|xmlFlatList", files)

# Need xmlStopParser to have the call to xmlStopParser commented out.
# That causes a problem.

# grep("addNode|xmlFlatList", files) works

sapply(files, function(x) {cat("***", x, "\n") ; gc() ; try(source(x))})

if(FALSE)
 sapply(files, function(x) {cat("***", x, "\n") ; gc() ; sapply(1:5, function(i) {cat(i, " ") ; try(source(x));}) ;cat("\n")})


# okay
# (i.e. run with gctorture())
# xmlFlatListTree.R

# not okay xmlNamespaceDefinitions.R


# XXX xmlStopParser.R
#gctorture()
#source("/Users/duncan/Rpackages/XML/R-ex/xmlStopParser.R")

