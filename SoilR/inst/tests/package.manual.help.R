#
# vim:set ff=unix expandtab ts=2 sw=2:
test.op=function(){
    library(SoilR)
    str <- as.character(?SoilR)
    #print(str)
    # note that is not the help text but (in case there is help)
    # the path to the file.
    # If the help file is missing then str will not contain the word "SoilR"
    # grep will not report 1 as the position in the character vector 
    # str but integer(0)  the test will fail
    checkEquals(grep("SoilR",str),1)
 }
