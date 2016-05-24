# Copyright (c) 2007 Simon Urbanek
#
# proj4 R package, License: GPL v2

# convert projection argument into PROJ.4 arguments string
# accepts either a named vector c(proj='merc',units='m'),
# a vector of parameters c('+proj=merc','+units=m')
# or a single string "+proj= +units=m"
#
# datum.default is added as "+datum=.." if not NA/NULL and +datum
# or +ellps doesn't exist in proj
# the same applies to ellps.default (in that order, i.e. datum has
# a higher precedence)
.proj2char <- function(proj, ellps.default=NA, datum.default=NA) {
    if (length(names(proj))) {
        proj <- paste('+',names(proj),'=',proj,sep='',collapse='\n')
        # remove spaces in all arguments
        proj <- gsub('\n',' ',gsub(' ','',proj))
    } else {
        if (length(proj) > 1)
            proj <- paste(as.character(proj), collapse=' ')
    }
    if (!is.character(proj)) proj <- as.character(proj)
    if (!is.null(datum.default) && !is.na(datum.default) && !length(grep("\\+datum=",proj)) && !length(grep("\\+ellps=",proj)))
        proj <- paste(proj," +datum=",datum.default,sep='')
    if (!is.null(ellps.default) && !is.na(ellps.default) && !length(grep("\\+datum=",proj)) && !length(grep("\\+ellps=",proj)))
        proj <- paste(proj," +ellps=",ellps.default,sep='')
    proj
}
