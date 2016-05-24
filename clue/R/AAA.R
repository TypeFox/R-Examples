## Things which must come first in the package code.

### * Internal utilities.

.false <- function(x) FALSE
.true <- function(x) TRUE

## A fast version of structure().
.structure <-
function(x, ...)
    `attributes<-`(x, c(attributes(x), list(...)))
