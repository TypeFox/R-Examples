wdSetProperty <-
function (property,value,object=wdapp[['Selection']],wdapp=.R2wd)
{
    if (length(property)>1)
      for (prop in property[1:(length(property)-1)]) {
            if (is.character(prop)) object<-object[[prop]]
            if (is.numeric(prop)) object$item(prop)
        }
        object[[property[length(property)]]]<-value
    invisible()
}
