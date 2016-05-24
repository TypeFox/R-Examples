wdGetProperty <-
function (property,object=wdapp[['Selection']],wdapp=.R2wd)
{
    if (length(property)>1)
        for (prop in property) {
            if (is.character(prop)) object<-object[[prop]]
            if (is.numeric(prop)) object$item(prop)
        }
     object
}
