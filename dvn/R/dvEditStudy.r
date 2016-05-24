dvEditStudy <-
function(   objectid, xmlfile, dv=getOption('dvn'),
            user=getOption('dvn.user'), pwd=getOption('dvn.pwd'), browser=FALSE, ...){
    if(is.null(user) | is.null(pwd))
        stop('Must specify username (`user`) and password (`pwd`)')
    if(is.null(xmlfile) || !is.character(xmlfile))
        stop('`xmlfile` must be xml character string or path to xml file')
    if(file.exists(xmlfile))
        filetosend <- charToRaw(paste(readLines(xmlfile),collapse=""))
    else
        filetosend <- charToRaw(xmlfile)
    xml <- dvDepositQuery(  query=paste('edit/study/',objectid,sep=''), user=user, pwd=pwd, dv=dv, browser=browser,
                            httpverb='PUT', postfields=filetosend,
                            httpheader=c('Content-Type'='application/atom+xml'))
    if(is.null(xml))
        return(NULL)
    if(browser==FALSE)
        .dvParseAtom(xml)
}
