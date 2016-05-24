.as_int<-function(x)
{if(!is.integer(x)) {
            xo = round(x)
            if(any(x > .Machine$integer.max))
                stop("the input has entries too large to be integer")
            if(!identical(TRUE, (ax <- all.equal(xo, x))))
                warning("the argument is not integer")
            else x=xo
        }
x
}
