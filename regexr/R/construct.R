#' Construct Human Readable Regular Expressions
#' 
#' This function is used to construct human readable regular expressions from
#' sub-expressions.  The user may provide additional meta information about each 
#' sub-expression. This meta information is an optional name and comment for the 
#' sub-expressions.  This allows one to write regular expressions in a fashion 
#' similar to writing code, that is the regular expression is written top to 
#' bottom, the syntax is broken up into manageable chunks, the sub-expressions 
#' can be indented to give structural insight such as nested groups.  Finally,
#' sub-expressions can be commented to provide linguistic grounding for more 
#' complex sub-expressions.
#' 
#' @param \ldots A series of comma separated character strings (sub-expressions) 
#' that may optionally be named, commented (see \code{?`\%:)\%`}, and indented.
#' @return Returns a character vector of the class \code{regexr}. The attributes 
#' of the returned object retain the original name and comment properties.
#' @keywords regex
#' @export
#' @examples
#' ## Minimal Example
#' minimal <- construct("a", "b", "c")
#' minimal
#' unglue(minimal)
#' comments(minimal)
#' subs(minimal)
#' test(minimal)
#' summary(minimal)
#' 
#' ## Example 1
#' m <- construct(
#'     space =   "\\s+"              %:)%  "I see",
#'     simp =    "(?<=(foo))",
#'     or =      "(;|:)\\s*"         %:)%  "comment on what this does",
#'     is_then = "[ia]s th[ae]n"
#' )
#' 
#' m
#' unglue(m)
#' summary(m)
#' subs(m)
#' comments(m)
#' subs(m)[4] <- "(FO{2})|(BAR)"
#' summary(m)
#' test(m)
#' \dontrun{
#' subs(m)[5:7] <- c("(", "([A-Z]|(\\d{5})", ")")
#' test(m)
#' }
#' 
#' library(qdapRegex)
#' ## explain(m)
#' 
#' ## Example 2 (Twitter Handle 2 ways)
#' ## Bigger Sub-expressions
#' twitter <- construct(
#'   no_at_wrd = "(?<![@@\\w])"            %:)%  "Ensure doesn't start with @@ or a word",
#'   at =        "(@@)"                    %:)%  "Capture starting with @@ symbol",
#'   handle =    "(([a-z0-9_]{1,15})\\b)"  %:)%  "Any 15 letters, numbers, or underscores"  
#' )
#' 
#' ## Smaller Sub-expressions
#' twitter <- construct(
#'   no_at_wrd = "(?<![@@\\w])"          %:)%  "Ensure doesn't start with @@ or a word",
#'   at =        "(@@)"                  %:)%  "Capture starting with @@ symbol",
#'   
#'   s_gr1 =     "("                     %:)%  "GROUP 1 START",            
#'       handle =    "([a-z0-9_]{1,15})"       %:)%  "Any 15 letters, numbers, or underscores", 
#'       boundary =  "\\b",
#'   e_gr1 =     ")"                      %:)%"GROUP 1 END"
#' )
#' 
#' twitter
#' unglue(twitter)
#' comments(twitter)
#' subs(twitter)
#' summary(twitter)
#' test(twitter)
#' ## explain(twitter)
#' 
#' x <- c("@@hadley I like #rstats for #ggplot2 work.",
#'     "Difference between #magrittr and #pipeR, both implement pipeline operators for #rstats:
#'         http://renkun.me/r/2014/07/26/difference-between-magrittr-and-pipeR.html @@timelyportfolio",
#'     "Slides from great talk: @@ramnath_vaidya: Interactive slides from Interactive Visualization
#'         presentation #user2014. http://ramnathv.github.io/user2014-rcharts/#1",
#'     "tyler.rinker@@gamil.com is my email",
#'     "A non valid Twitter is @@abcdefghijklmnopqrstuvwxyz"
#' )
#' 
#' library(qdapRegex)
#' rm_default(x, pattern = twitter, extract = TRUE)
#' 
#' ## Example 3 (Modular Sub-expression Chunks)
#' combined <- construct(
#'     twitter = twitter               %:)%"Twitter regex created previously",
#'     or =      "|"                   %:)%"Join handle regex & hash tag regex",
#'     hash =    grab("@@rm_hash")     %:)%"Twitter hash tag regex"
#' )
#' 
#' combined
#' unglue(combined)
#' comments(combined)
#' subs(combined)
#' summary(combined)
#' test(combined)
#' ## explain(combined)
#' 
#' ## Different Structure (no names): Example from Martin Fowler: 
#' ## *Note: Fowler argues for improved choices in regex representation
#' ## and names that make the regex functionality more evident, commenting
#' ## only where needed. See:
#' ## browseURL("http://martinfowler.com/bliki/ComposedRegex.html")
#' 
#' pattern <- construct(
#'     '@@"^score',
#'     '\\s+',
#'     '(\\d+)'          %:)% 'points',
#'     '\\s+',
#'     'for',
#'     '\\s+',
#'     '(\\d+)'          %:)% 'number of nights',
#'     '\\s+',
#'     'night'           ,
#'     's?'              %:)% 'optional plural',
#'     '\\s+',
#'     'at',
#'     '\\s+',
#'     '(.*)'            %:)% 'hotel name',
#'     '";'
#' )
#'   
#' summary(pattern)
construct <- function(...){
    out <- paste0(...)
    class(out) <- c("regexr", class(out))
    attributes(out)[["subs"]] <- list(...)
    attributes(out)[["comments"]] <- lapply(list(...), get_comment)
    out
}
