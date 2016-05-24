## ----echo=FALSE, warning=FALSE-------------------------------------------
knitr::opts_knit$set(out.format = "latex")
so.d <- knitr::knit_theme$get("solarized-dark")
so.l <- knitr::knit_theme$get("solarized-light")
knitr::knit_theme$set(so.l)
knitr::opts_chunk$set(concordance=TRUE)
# knitr::opts_chunk$set(message = FALSE, warning = FALSE)
knitr::opts_chunk$set(out.width = '0.5\\linewidth', fig.align = "center", fig.show = 'asis')

## ----poppr_funk, eval = TRUE, echo = FALSE-------------------------------
print_command <- function(funk){
  fargs <- formals(funk)
  
  lapply(names(fargs), function(arg_name, fargs){
    arg <- fargs[[arg_name]]
    if (missing(arg)){
      fargs[[arg_name]] <<- as.symbol(arg_name)
      names(fargs)[names(fargs) == arg_name] <<- ""
    }
  }, fargs)
  fargs$call <- as.symbol(funk)
  fargs <- fargs[c(length(fargs), 1:(length(fargs) - 1))]
  return(as.call(fargs))
}


# This is an internal function from knitr
my.color.block = function(color1 = '', color2 = '') {
  function(x, options) {
    x = gsub('\n*$', '', x)
    x = my_escape_latex(x, newlines = TRUE, spaces = TRUE)
    # babel might have problems with "; see http://stackoverflow.com/q/18125539/559676
    x = gsub('"', '"{}', x)
    sprintf('\n\n{\\ttfamily\\noindent%s%s%s}', color1, x, color2)
  }
}

# This is also an internal function from knitr
my_escape_latex = function(x, newlines = FALSE, spaces = FALSE) {
  x = gsub('\\\\', '\\\\textbackslash', x)
  x = gsub('([#$%&_{}])', '\\\\\\1', x)
  x = gsub('\\\\textbackslash', '\\\\textbackslash{}', x)
  x = gsub('~', '\\\\textasciitilde{}', x)
  x = gsub('\\^', '\\\\textasciicircum{}', x)
  if (newlines) x = gsub('(?<!\n)\n(?!\n)', '\\\\\\\\', x, perl = TRUE)
  if (spaces) x = gsub('  ', '\\\\ \\\\ ', x)
  x
}

## ----sl------------------------------------------------------------------
rnorm(10) # Works!

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.d)

## ----sd, error = TRUE----------------------------------------------------
try(rnorm("A")) # Throws a warning and error.

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.l)

## ----old_genind, message = TRUE, warning = TRUE, error = TRUE------------
library("poppr")
data("old_partial_clone", package = "poppr") # From version 1.1.5
data("partial_clone", package = "poppr")     # From version 2.0.0

## ----show_genind, message = TRUE, warning = TRUE, error = TRUE-----------
names(attributes(old_partial_clone)) # Has pop.names, ind.names, etc.
names(partial_clone)                 # Has strata slot.
# This is ultimately what we want
partial_clone

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.d)

## ----show_old_genind, message = TRUE, warning = TRUE, error = TRUE-------
try(old_partial_clone)

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.l)

## ----old2new_genind, message = TRUE, warning = TRUE, error = TRUE--------
opc <- old2new_genind(old_partial_clone)
opc # It prints!

## ----old_Pinf, message = TRUE, warning = TRUE, error = TRUE--------------
data("old_Pinf", package = "poppr") # From version 1.1.5
data("Pinf", package = "poppr")     # From version 2.0.0
names(attributes(old_Pinf)) # No strata slot
names(Pinf)                 # Has strata slot
Pinf # What we want

## ----old2new_genclone, message = TRUE, warning = TRUE, error = TRUE------
opi <- old2new_genclone(old_Pinf)
opi # It prints!

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.d)

## ----hier_warning--------------------------------------------------------
head(gethierarchy(Pinf, ~Continent/Country)) # Throws warning

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.l)

## ----hier_strata---------------------------------------------------------
head(strata(Pinf, ~Continent/Country))       # No warning :)

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.d)

## ----setpop_warning------------------------------------------------------
setpop(Pinf) <- ~Country # Ambiguous warning

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.l)

## ----setPop_nowarning----------------------------------------------------
setPop(Pinf) <- ~Country # No warning

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.d)

## ----hier_no, eval = FALSE-----------------------------------------------
#  myData@hierarchy <- myData@hierarchy[-2] # wat

## ----Pinf_wrong_strata, error = TRUE, message = TRUE, warning = TRUE-----
newPinf        <- Pinf
the_strata     <- head(newPinf@strata)
newPinf@strata <- the_strata     # Only setting strata for six samples!
newPinf@strata                   # How is this allowed?
try(setPop(newPinf) <- ~Country) # Oh.

## ----echo = FALSE--------------------------------------------------------
knitr::knit_hooks$set(message = my.color.block('\\bfseries\\color{errorcolor}{', '}'))
message(geterrmessage())

## ----echo = FALSE--------------------------------------------------------
knitr::knit_theme$set(so.l)

## ----Pinf_right_strata, error = TRUE, message = TRUE, warning = TRUE-----
the_strata       <- head(strata(Pinf))
try(strata(Pinf) <- the_strata)

## ----echo = FALSE--------------------------------------------------------
knitr::knit_hooks$set(message = my.color.block('\\bfseries\\color{errorcolor}{', '}'))
message(geterrmessage())

## ----echo = FALSE--------------------------------------------------------
knitr::knit_hooks$set(message = my.color.block('\\itshape\\color{messagecolor}{', '}'))

