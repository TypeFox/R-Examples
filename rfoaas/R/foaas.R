##  rfoaas -- An R interface to the FOAAS service
##
##  Copyright (C) 2014 - 2015  Dirk Eddelbuettel <edd@debian.org>
##
##  This file is part of rfoaas
##
##  rfoaas is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 2 of the License, or
##  (at your option) any later version.
##
##  rfoaas is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with rfoaas.  If not, see <http://www.gnu.org/licenses/>.

.foaas <- function(..., filter=.filter(), language=.language()) {

    ## -- The following used to work when foaas.com was running with text/plain in default
    ##
    #req <- URLencode(paste(..., sep="/"))     	         	# collate arguments and encode
    #con <- url(paste0("http://foaas.herokuapp.com/", req)) 	# form url and create connection
    #res <- readLines(con, n=n, warn=FALSE)       		# read one line from connection
    #close(con)                                                 # clean connection
    #Encoding(res) <- "UTF-8"    				# server-side is UTF-8, needed on Windows 
    #res
    ##
    ## -- but now we have to explicitly request it via accept headers, so we need http::GET
    
    #srv <- "http://foaas.herokuapp.com"
    srv <- "http://foaas.com"
    req <- paste(srv, ..., sep="/")		     	        # collate normal arguments 

    ## deal with optional arguments by test and conditional appends
    supargs <- c(ifelse(language=="none", "", paste0("i18n=", language)),
                 switch(filter,
                        "none" = "",
                        "shoutcloud" = "shoutcloud"))
    if (any(supargs != "")) {           			# if we have arguments
        supargs <- paste(supargs, collapse="&")                 # collate them, but ...
        supargs <- gsub("&$", "", gsub("^&", "", supargs)) 	# ... nuke leading or trailing '&'
        req <- paste(req, supargs, sep="?")                     # and append
    }

    req <- URLencode(req)					# encode as a URL just in case
    res <- GET(req, accept("text/plain"))
    txt <- content(res, "text", encoding="utf-8")
    txt
}

.from <- function() {
    getOption("rfoaasFrom", Sys.info()["user"])
}

.filter <- function(filter=c("none", "shoutcloud")) {
    filter <- getOption("rfoaasFilter", "none")
    filter <- match.arg(filter)
}

.language <- function(language) {
    if (missing(language)) language <- "none"
    language                            # could do other tests here
}


## 'meta' query one -- returns a version string
version     <- function()                      { .foaas("version") }

## 'meta' query two -- returns JSON object descriting name, url, and fields on available queries
## As this returns JSON, use RJSONIO or jsonlite to deal with the result
operations  <- function()                      { .foaas("operations") }

off         <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("off", name, from, filter=filter, language=language) }
you         <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("you", name, from, filter=filter, language=language) }
this        <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("this", from, filter=filter, language=language) }
that        <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("that", from, filter=filter, language=language) }
everything  <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("everything", from, filter=filter, language=language) }
everyone    <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("everyone", from, filter=filter, language=language) }
donut       <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("donut", name, from, filter=filter, language=language) }
shakespeare <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("shakespeare", name, from, filter=filter, language=language) }
linus       <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("linus", name, from, filter=filter, language=language) }
king        <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("king", name, from, filter=filter, language=language) }
pink        <- function(name, filter=.filter(), language=.language())                              { .foaas("pink", name, filter=filter, language=language) }
life        <- function(name, filter=.filter(), language=.language())                              { .foaas("life", name, filter=filter, language=language) }
chainsaw    <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("chainsaw", name, from, filter=filter, language=language) }
outside     <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("outside", name, from, filter=filter, language=language) }
thanks      <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("thanks", from, filter=filter, language=language) }
flying      <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("flying", from, filter=filter, language=language) }
fascinating <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("fascinating", from, filter=filter, language=language) }
madison     <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("madison", name, from, filter=filter, language=language) }
cool        <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("cool", from, filter=filter, language=language) }
field       <- function(name, from=.from(), reference, filter=.filter(), language=.language())     { .foaas("field", name, from, reference, filter=filter, language=language) }
nugget      <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("nugget", name, from, filter=filter, language=language) }
yoda        <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("yoda", name, from, filter=filter, language=language) }
ballmer     <- function(name, company, from=.from(), filter=.filter(), language=.language())       { .foaas("ballmer", name, company, from, filter=filter, language=language) }
what        <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("what", from, filter=filter, language=language) }
because     <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("because", from, filter=filter, language=language) }
caniuse     <- function(tool, from=.from(), filter=.filter(), language=.language())                { .foaas("caniuse", tool, from, filter=filter, language=language) }
bye         <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("bye", from, filter=filter, language=language) }
diabetes    <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("diabetes", from, filter=filter, language=language) }
bus         <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("bus", from, filter=filter, language=language) }
xmas        <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("xmas", name, from, filter=filter, language=language) }
awesome     <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("awesome", from, filter=filter, language=language) }
tucker      <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("tucker", from, filter=filter, language=language) }
bucket      <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("bucket", from, filter=filter, language=language) }
bday        <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("bday", name, from, filter=filter, language=language) }
family_     <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("family", from, filter=filter, language=language) }
shutup      <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("shutup", name, from, filter=filter, language=language) }
zayn        <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("zayn", from, filter=filter, language=language) }
dalton      <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("dalton", name, from, filter=filter, language=language) }
dosomething <- function(do, something, from=.from(), filter=.filter(), language=.language())       { .foaas("dosomething", do, something, from, filter=filter, language=language) }
#off_with    <- function(behaviour, from=.from(), filter=.filter(), language=.language())           { .foaas("off_with", behaviour, from, filter=filter, language=language) }
retard      <- function(from=.from(), filter=.filter(), language=.language())                      { .foaas("retard", from, filter=filter, language=language) }
thumbs      <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas("thumbs", name, from, filter=filter, language=language) }

## catch-all 
thing       <- function(name, from=.from(), filter=.filter(), language=.language())                { .foaas(name, from, filter=filter, language=language) }






