wdNewDoc <-
function (name=NULL,wdapp=.R2wd)
{
    if (is.null(wdapp[['Documents']])) stop('broken handle to Word. Reset using wdGet')
    wdapp[['Documents']]$Add()
    if (!is.null(name)) wdSave(name)
    invisible()
}

