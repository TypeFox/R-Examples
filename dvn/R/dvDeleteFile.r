dvDeleteFile <-
function(   fileid, dv=getOption('dvn'),
            user=getOption('dvn.user'), pwd=getOption('dvn.pwd'), browser=FALSE, ...){
    if(is.null(user) | is.null(pwd))
        stop('Must specify username (`user`) and password (`pwd`)')
    xml <- dvDepositQuery(  query=paste('edit-media/file/',fileid,sep=''), user=user, pwd=pwd, dv=dv, browser=browser,
                            httpverb='DELETE')
    if(is.null(xml))
		invisible(NULL)
	if(browser==FALSE){
		if(xml=='')
            message('Operation appears to have succeeded.')
        return(xml)
    }
}
