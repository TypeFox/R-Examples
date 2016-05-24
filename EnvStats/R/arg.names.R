arg.names <-
function (x) 
{
    switch(mode(x), `function` = names(formals(x)), call = names(x)[-1], 
        character = {
            names(formals(x))
        }, stop("mode must be 'function', 'call', or 'character'"))
}
