convertArgs2String = function()
##title<< Save function argument values to a character string
##description<< convertArgs2Strings saves function argument settings to a character string.
##details <<
## convertArgs2Strings saves function argument settings to a character string.
## This is used in the ncdf routines that save these argument values as
## attributes of the target ncdf files.
{
    args.all  <- ls(parent.frame())
    args.call <- list()
    for (g in 1:length(args.all))
        args.call[[g]]<- get( args.all[g],envir=parent.frame())
    names(args.call)<- args.all
    extract.names=function(x)
    {
        if (class(x)[1]=='function') {
            text=body(x)
        } else {
            text = x
        }
        text
    }
    string.args  <- paste(paste(names(args.call),sapply(args.call,extract.names)
                                   ,sep=': '),collapse='; ')
    ##value<< character string: arguments given to function converted into a string.         
    return(string.args)
}
