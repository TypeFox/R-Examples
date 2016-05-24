dvCreateStudy <-
function(   dvname, xmlfile, dv=getOption('dvn'),
            user=getOption('dvn.user'), pwd=getOption('dvn.pwd'), browser=FALSE, ...){
    if(is.null(user) | is.null(pwd))
        stop('Must specify username (`user`) and password (`pwd`)')
    if(is.null(xmlfile) || !is.character(xmlfile))
        stop('`xmlfile` must be xml character string or path to xml file')
    if(inherits(dvname,'dvServiceDoc')){
        tmp <- dvname$dataverses$dvn
        if(length(tmp)>1)
            stop('Multiple dataverses available for this user. Please supply dataverse name as character string.')
        dvname <- tmp[1]
    }
    if(file.exists(xmlfile))
        filetosend <- charToRaw(paste(readLines(xmlfile),collapse=""))
    else
        filetosend <- charToRaw(xmlfile)
    xml <- dvDepositQuery(  query=paste('collection/dataverse/',dvname,sep=''), user=user, pwd=pwd, dv=dv, browser=browser,
                            httpverb='POST', postfields=filetosend,
                            httpheader=c('Content-Type'='application/atom+xml'))
    if(is.null(xml))
         return(NULL)
    if(browser==FALSE)
         .dvParseAtom(xml)
}
