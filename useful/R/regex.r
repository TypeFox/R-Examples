## Written by Jared P. Lander
## See LISCENSE for copyright information
## make getting regex results easier
## looks like Hadley's stringr pretty much does this already, but ok

## takes the result of a regex and extracts the desired tect
## @results () the return of either regexpr or gregexpr
## @text (character) the the text being substr'd
ProcessRegex <- function(results, text)
{
    # find the start and stop positions
    theFrame <- data.frame(Start=results, Stop=results + attr(results, "match.length") - 1, Text=text, stringsAsFactors=FALSE)

    # just keep the text where the pattern was found
    theFrame <- theFrame[theFrame$Start != -1, ]

    # extract the desired text and return it
    return(with(theFrame, substr(Text, Start, Stop)))
}


## @pattern (character) the regular expression pattern to search for
## @text (character) 
## @ignore.case (logical) If FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.
## @perl (logical) Should perl-compatible regexps be used? Has priority over extended.
## @fixed (logical) If TRUE, pattern is a string to be matched as is. Overrides all conflicting arguments.
## @useBytes (logical) If TRUE the matching is done byte-by-byte rather than character-by-character.
## See ?regexpr for details
regex <- function(pattern, text, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
{
    # run the regex
    theResult <- regexpr(pattern=pattern, text=text, ignore.case=ignore.case, perl=perl, fixed=fixed, useBytes=useBytes)

    # extract the text
    return(ProcessRegex(results=theResult, text=text))
}


## @pattern (character) the regular expression pattern to search for
## @text (character) 
## @ignore.case (logical) If FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.
## @perl (logical) Should perl-compatible regexps be used? Has priority over extended.
## @fixed (logical) If TRUE, pattern is a string to be matched as is. Overrides all conflicting arguments.
## @useBytes (logical) If TRUE the matching is done byte-by-byte rather than character-by-character.
## See ?regexpr for details
OneRegex <- function(text, pattern, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
{
    # get the gregex for a single text entry
    theResult <- gregexpr(pattern=pattern, text=text, ignore.case=ignore.case, perl=perl, fixed=fixed, useBytes=useBytes)[[1]]

    # extract the text
    return(ProcessRegex(results=theResult, text=text))
}


## same as regex but for global pattern recognition
## does one gregexpr for each element of text
## results are a list
## @pattern (character) the regular expression pattern to search for
## @text (character) 
## @ignore.case (logical) If FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.
## @perl (logical) Should perl-compatible regexps be used? Has priority over extended.
## @fixed (logical) If TRUE, pattern is a string to be matched as is. Overrides all conflicting arguments.
## @useBytes (logical) If TRUE the matching is done byte-by-byte rather than character-by-character.
## See ?regexpr for details
#' @importFrom plyr alply
gregex <- function(pattern, text, ignore.case=FALSE, perl=FALSE, fixed=FALSE, useBytes=FALSE)
{
    # get the results from the gregexpr's
    theText <- alply(text, 1, OneRegex, pattern=pattern, ignore.case=ignore.case, perl=perl, fixed=fixed, useBytes=useBytes)

    # give the list elements the names of the text
    names(theText) <- text
    
    # get rid of empty returns
    theText <- theText[which(llply(theText, length) != 0)]

    return(theText)
}
