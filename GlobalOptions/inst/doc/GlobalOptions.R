## ---- echo = FALSE, message = FALSE---------------------------------------------------------------
library(markdown)
options(markdown.HTML.options = c(options('markdown.HTML.options')[[1]], "toc"))

library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = FALSE,
    message = FALSE,
    fig.align = "center")
options(markdown.HTML.stylesheet = "custom.css")

options(width = 100)

## -------------------------------------------------------------------------------------------------
library(GlobalOptions)
opt = setGlobalOptions(
    a = 1,
    b = "text"
)

## -------------------------------------------------------------------------------------------------
opt()
opt("a")
opt$a
op = opt()
op
opt(a = 2, b = "new text")
opt()
opt$b = ""
opt()
opt(op)
opt()

## -------------------------------------------------------------------------------------------------
opt(a = 2, b = "new text")
opt(RESET = TRUE)
opt()

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = list(.value = 1,
             .length = c(1, 3),
             .class = "numeric")
)

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = list(.value = 1,
             .read.only = TRUE),
    b = 2
)
opt(READ.ONLY = TRUE)
opt(READ.ONLY = FALSE)
opt(READ.ONLY = NULL)  # default, means return both

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    verbose = 
        list(.value = TRUE,
             .filter = function(x) {
                 if(is.null(x)) {
                     return(FALSE)
                 } else if(is.na(x)) {
                     return(FALSE)
                 } else {
                     return(x)
                 }
              })
)
opt(verbose = FALSE); opt("verbose")
opt(verbose = NA); opt("verbose")
opt(verbose = NULL); opt("verbose")

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    margin = 
        list(.value = c(1, 1, 1, 1),
             .length = c(1, 2, 4),
             .filter = function(x) {
                if(length(x) == 1) {
                    return(rep(x, 4))
                } else if(length(x) == 2) {
                    return(rep(x, 2))
                } else {
                    return(x)
                }
            })
)
opt(margin = 2); opt("margin")
opt(margin = c(2, 4)); opt("margin")

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    prefix = ""
)
opt(prefix = function() paste("[", Sys.time(), "] ", sep = " "))
opt("prefix")
Sys.sleep(2)
opt("prefix")

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    test_fun = list(.value = function(x1, x2) t.test(x1, x2)$p.value,
                    .class = "function")
)
opt(test_fun = function(x1, x2) cor.test(x1, x2)$p.value)
opt("test_fun")

## -------------------------------------------------------------------------------------------------
lt = setGlobalOptions(
    a = list(.value = 1),
    b = list(.value = function() 2 * get_opt_value('a')),
    get_opt_value_fun = TRUE
)
opt = lt$opt_fun
# the variable name here should be same as the one in above
get_opt_value = lt$get_opt_value

opt("b")
opt(a = 2)
opt("b")
opt(a = 2, b = 3)
opt()

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = 1
)

opt(LOCAL = TRUE)
opt(a = 2)
opt("a")
opt(LOCAL = FALSE)
opt("a")

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = 1
)

f1 = function() {
    opt(LOCAL = TRUE)
    opt(a = 2)
    return(opt("a"))
}
f1()
opt$a

f2 = function() {
    opt(LOCAL = TRUE)
    opt(a = 4)
    return(opt("a"))
}
f2()
opt$a

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = 1
)

f1 = function() {
    opt(LOCAL = TRUE)
    opt(a = 2)
    return(f2())
}

f2 = function() {
    opt("a")
}

f1()
opt("a")

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
	a = list(.value = 1,
	         .visible = FALSE),
	b = 2
)
opt()
opt("a")
opt(a = 2)
opt("a")
opt()

## -------------------------------------------------------------------------------------------------
opt = setGlobalOptions(
    a = list(.value = 1,
             .private = TRUE)
)
require(stats)
ns = getNamespace("stats")
environment(opt)$options$a$`__generated_namespace__` = ns

## -------------------------------------------------------------------------------------------------
args(opt)

## -------------------------------------------------------------------------------------------------
opt1 = setGlobalOptions(
    a = list(.value = 1)
)
opt2 = setGlobalOptions(
    a = list(.value = 1)
)
opt1(a = 2)
opt1("a")
opt2("a")

## -------------------------------------------------------------------------------------------------
sessionInfo()

