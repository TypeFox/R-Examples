# Project: pkgmaker
# 
# Author: renaud gaujoux
# Created: Oct 23, 2013
###############################################################################


#' Knitr Extensions
#' 
#' \code{knit_ex} is a utility function for running small knitr examples, 
#' e.g., to illustrate functionalities or issues. 
#' 
#' @param x text to knit as a character vector
#' @param ... arguments passed to \code{\link[knitr]{knit2html}} or \code{\link[knitr]{knit}}
#' @param quiet logical that indicates if knitting should be quiet (no progress bars etc..).
#' @param open logical, only used when \code{x} is in .Rmd format, that indicates 
#' if the generated document result should be open in a browse, instead of 
#' being printed on screen.
#' Not that a browser will not open in non-interactive sessions, and the result will
#' be returned invisibly.
#' 
#' @return 
#' \code{knit_ex} returns the generated code, although invisibly when \code{open=TRUE}.
#' 
#' @export
#' @examples
#' 
#' library(knitr)
#' knit_ex("1 + 1")
#' 
knit_ex <- function(x, ..., quiet = TRUE, open = FALSE){
    
    library(knitr)
    # substitute special markup for Rmd markup (necessary for knit_ex examples)
    x <- gsub("^^^", "```", x, fixed = TRUE)
    
    if( !(html_chunks <- any(grepl("```{", x, fixed = TRUE))) ){
        if( all(!grepl(">>=", x, fixed = TRUE)) ){
            x <- c("```{r}", x, "```")
            html_chunks <- TRUE   
        }
    }
    x <- paste0(x, collapse = "\n")
    if( any(html_chunks) ){
        res <- knit2html(text = x, ..., fragment.only = TRUE, quiet = quiet)
        if( open ){
            tmp <- tempfile("knit_ex", fileext = '.html')
            cat(res, file = tmp, sep = "\n") 
            if( interactive() ) browseURL(tmp)
            return(invisible(res))
        }
    }else{
        res <- knit(text = x, ..., quiet = quiet)
    }
    cat(res)
}

try_message <- function(signal = FALSE){
    function(expr){
        tryCatch(expr, error = function(e){
                if( signal ) message(e)
                else message('Error: ', conditionMessage(e))
                invisible()
            })
    }
}

#' \code{hook_try} is a knitr hook to enable showing error 
#' messages thrown by \code{\link{try}}.
#' The function is not meant to be called directly, but only registered 
#' using \code{\link{knit_hooks}} (see details on this dedicated man page).
#' 
#' \code{hook_try} simply defines a function \code{try} in \code{envir} that prints 
#' the error message if any, and is called instead of base \code{\link{try}}. 
#' 
#' @param before logical that indicates when the hook is being called: 
#' before or after the chunk is processed.
#' @param options list of current knitr chunk options 
#' @param envir environment where the chunk is evaluated
#' 
#' @rdname knit_ex
#' @export
#' @examples
#' 
#' library(knitr)
#' 
#' # standard error message is caught
#' knit_ex("stop('ah ah')")
#' 
#' # with try the error is output on stderr but not caughted by knitr
#' knit_ex("try( stop('ah ah') )")
#' 
#' # no message caught
#' knit_ex("
#' ^^^{r, include = FALSE}
#' knit_hooks$set(try = pkgmaker::hook_try)
#' ^^^
#' 
#' ^^^{r, try=TRUE}
#' try( stop('ah ah') )
#' ^^^")
#' 
hook_try <- local({
    .try_defined <- FALSE
    function(before, options, envir){
    
        # remove hacked version of try
        if( !before ){
            if( .try_defined && exists('try', envir = envir, inherits = FALSE) ){
                remove(list = 'try', envir = envir)
            }
            .try_defined <<- FALSE
            return(invisible())
        }
        
        if( !is.null(options$try) ){
            
            # signal
            do.signal <- isFALSE(options$try)
            if( isManualVignette() && isTRUE(options$try) ){
                do.signal <- TRUE
            }
            # define hacked version of try()
            .try <- try_message(do.signal)
            assign('try', .try, envir)
            .try_defined <<- TRUE
        }
    }
})


chunkOutputHook <- function(name, hook, type = c('output', 'source', 'chunk')){
    type <- match.arg(type)
    function(){
        .hook_bkp <- NULL
        function(before, options, envir){
            # do nothing if the option is not ON
            if( is.null(options[[name]]) ) return()
            
            # set/unset hook
            if( before ){
                # store current hook function
                if( is.null(.hook_bkp) ) .hook_bkp <<- knit_hooks$get(type)
                
                # define hook wrapper
                hook_wrapper <- function(x, options){
                    res <- .hook_bkp(x, options)
                    hook(res, options)
                }
                        
                args <- list()
                args[[type]] <- hook_wrapper
                do.call(knit_hooks$set, args)
            }else{
                args <- list()
                args[[type]] <- .hook_bkp
                do.call(knit_hooks$set, args)
                .hook_bkp <<- NULL
            }
        }
    }
}

#' Knitr Hook for Handling Backspace Character
#' 
#' \code{hook_backspace} is a chunk hook that enables the use of backspace
#' characters in the output (e.g., as used in progress bars), and still 
#' obtain a final output as in the console.
#' 
#' @rdname knit_ex
#' @export 
#' @examples 
#' 
#' # Correctly formatting backspaces in chunk outputs
#' tmp <- tempfile(fileext = '.Rmd')
#' cat(file = tmp, "
#' ^^^{r, include = FALSE}
#' library(knitr)
#' knit_hooks$set(backspace = pkgmaker::hook_backspace())
#' ^^^
#' Default knitr does not handle backspace and adds a special character:
#' ^^^{r}
#' cat('abc\bd')
#' ^^^
#' 
#' Using the hook backspace solves the issue:
#' ^^^{r, backspace=TRUE}
#' cat('abc\bd')
#' ^^^
#' ")
#' 
#' # knit
#' out <- knitr::knit2html(tmp, fragment.only = TRUE)
#' # look at output
#' \dontrun{
#'   browseURL(out)
#'   edit( file = out)
#' }
#' # cleanup
#' unlink(c(tmp, out))
#' 
#' 
hook_backspace <- chunkOutputHook('backspace', 
        function(x, options){
            if( !isTRUE(options$backspace) ) x
            str_bs(x)
        }
)

#' \code{str_bs} substitutes backspace characters (\\b) to produce
#' a character string as it would be displayed in the console.
#'
#' @author
#' Renaud Gaujoux
#'  
#' \code{str_bs} was adapted from a proposal from Yihui Xie.
#'  
#' @rdname str_out 
#' @export
#' @examples 
#' 
#' # Backspace substitution
#' str_bs("abc")
#' str_bs("abc\b")
#' str_bs("abc\b\b")
#' str_bs("abc\bd")
#' str_bs("abc\b\bde\b")
#' 
#' # more complex example
#' x <- "\bab\nc\bd\n\babc\b\bd"
#' cat(x, "\n")
#' y <- str_bs(x)
#' y
#' cat(y, "\n")
#' 
str_bs <- function(x){
    # remove leading backspaces
    x <- gsub("^\b+", "", x)
    # remove backspaces at beginning of line
    x <- gsub("\n\b+", '\n', x)
    while( length(grep('\b', x, fixed = TRUE)) ) 
        x <- gsub('[^\n\b][\b]', '', x)
    
    x
}

md_toggleCode <- function(){
    cat(
"
<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js\"></script>
<script type=\"text/javascript\">
// toggle visibility of R source blocks in R Markdown output
function toggle_vis(o) {
    var d = o.style.display;
    o.style.display = (d == 'block' || d == '') ? 'none':'block';
  }
function toggle_R(what) {
  var x = document.getElementsByClassName('r');
  if (x.length == 0) return;

  for (i = 0; i < x.length; i++) {
    var y = x[i];
    switch (y.tagName.toLowerCase()) {
      case 'pre':
        toggle_vis(y);
        if( what == 'setup' ){
            y.id = \"Rcode_\" + i;
        }
        break;
      case 'code':
        var z = y.parentNode;
        // pandoc uses the class 'sourceCode r' on both pre and code
        if (z.tagName.toLowerCase() == 'pre' && z.className != 'sourceCode r') {
          toggle_vis(z);
          if( what == 'setup' ){
              z.id = \"Rcode_\" + i;
              var newContent = $(\"<a href=\\\"\\\" onclick=\\\"toggle_vis(document.getElementById('\" + z.id + \"')); return false;\\\">Show/Hide R code</a>\");
              newContent.insertBefore(z);
          }
        }
        break;
    }
  }
}

$( document ).ready(function(){
    toggle_R('setup');
});
</script>"
    )
}


.js_include_jquery <- "<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js\"></script>\n"
.js_toggle_fun <- "<script type=\"text/javascript\">
function toggle_vis2(cl) {
	elt = document.getElementsByClassName(cl)
	for (i = 0; i < elt.length; i++) {
		var o = elt[i].parentNode;
        var d = o.style.display;
        o.style.display = (d == 'block' || d == '') ? 'none':'block';
	}
  }

window.onload = function(){
	toggle_vis2('Chunk_none')	
}
</script>
"

.js_def_toggle_code <- "<script type=\"text/javascript\">
// toggle visibility of R source blocks in R Markdown output
function toggle_R(what) {
  var x = document.getElementsByClassName('r');
  if (x.length == 0) return;
  function toggle_vis(o) {
    var d = o.style.display;
    o.style.display = (d == 'block' || d == '') ? 'none':'block';
  }

  for (i = 0; i < x.length; i++) {
    var y = x[i];
    switch (y.tagName.toLowerCase()) {
      case 'pre':
        toggle_vis(y);
        if( what == 'setup' ){
        y.id = \"Rcode_\" + i;
        }
        break;
      case 'code':
        var z = y.parentNode;
        // pandoc uses the class 'sourceCode r' on both pre and code
        if (z.tagName.toLowerCase() == 'pre' && z.className != 'sourceCode r') {
            toggle_vis(z);
            if( what == 'setup' ){
                z.id = \"Rcode_\" + i;
                var newContent = $(\"<a href=\\\"\\\" onclick=\\\"$('#\" + z.id + \"').toggle(); return false;\\\">Show/Hide R code</a>\");
                newContent.insertBefore(z);
            }
    	}
    	break;
    }
  }
}

$( document ).ready(function(){
    toggle_R('setup');
});
</script>"

#' \code{hook_toggle} is a chunk hook that adds clickable elements to toggle \emph{indvidual}
#' code chunks in HTML documents generated from .Rmd files.
#' 
#' @rdname knit_ex
#' @export
#' @examples
#' 
#' knit_ex("
#' 
#' Declare chunk hook:
#' ^^^{r, setup}
#' library(knitr)
#' knit_hooks$set(toggle = hook_toggle())
#' ^^^
#' 
#' The R code of this chunk can be toggled on/off, and starts visible:
#' ^^^{r, toggle=TRUE}
#' print(1:10)
#' ^^^
#' The R code of this chunk can be toggled on/off, and starts hidden:
#' ^^^{r, toggle=FALSE}
#' print(1:10)
#' ^^^
#' 
#' This is a plain chunk that cannot be toggled on/off:
#' ^^^{r}
#' print(1:10)
#' ^^^
#' 
#' Now all chunks can be toggled and start visible:
#' ^^^{r, toggle_all}
#' opts_chunk$set(toggle = TRUE)
#' ^^^
#' 
#' ^^^{r}
#' sample(5)
#' ^^^
#' 
#' To diable the toggle link, one can pass anything except TRUE/FALSE:
#' ^^^{r, toggle = NA}
#' sample(5)
#' ^^^
#' 
#' ", open = TRUE)
#' 
hook_toggle <- function(){
    .init <- TRUE
    .last_label <- NULL
    
    fn <- chunkOutputHook('toggle', type = 'source', function(x, options){
        opt <- options$toggle
        label <- options$label
        if( !isTRUE(opt) && !isFALSE(opt) ) return(x)
#        print(x)
#        x <- gsub("^\n", '', x)
        # add javascript define for toggle function
        if( .init ){
            x <- paste0(.js_toggle_fun, x)
            .init <<- FALSE
        }
        
        disp <- if( opt ) 'Chunk_block'
                else 'Chunk_none'
        id <- paste0("Rcode_", label)
        subst <- paste0("```{", id, " \\1 ", disp, "}\n")
        if( !identical(label, .last_label) ){
            .last_label <<- label
            subst <- paste0("<a href=\"\" onclick=\"toggle_vis2('", id, "'); return false;\">Show/Hide R code</a>\n", subst)
        }
        sub("```([^\n]*)\n", sprintf(subst, 'block'), x)
    })
    fn()
}
