dvAddFile <-
function(   objectid, filename=NULL, dataframe=NULL, dv=getOption('dvn'),
            user=getOption('dvn.user'), pwd=getOption('dvn.pwd'), browser=FALSE, ...){
    if(is.null(user) | is.null(pwd))
        stop('Must specify username (`user`) and password (`pwd`)')
    if(is.null(filename) & is.null(dataframe))
        stop(   "Must specify `filename` as .zip or a vector of filenames,",
                " or `dataframe` as an R dataframe")
    if(!is.null(dataframe)){
        dataframe <- sapply(dataframe, function(x) {
            tmp <- file.path(tempdir(),paste(x,'RData',sep='.'))
            save(list=x, file=tmp)
            invisible(tmp)
        })
    }
    filename <- c(filename, dataframe)
    if(length(filename)>1 || (length(filename)==1 && !tools::file_ext(filename)=='zip')){
        tmpzip <- tempfile(fileext='.zip')
        zip(tmpzip, filename)
    }
    xml <- dvDepositQuery(  query=paste('edit-media/study/',objectid,sep=''),
            user=user, pwd=pwd, dv=dv, browser=browser,
            httpverb='POST', #upload=TRUE,
            postfields=readBin(tmpzip,what='raw',file.info(tmpzip)$size),
            httpheader=c(   'Content-Disposition'=paste('attachment; filename',tmpzip,sep='='),
                            'Content-Type'='application/zip',
                            'Packaging'='http://purl.org/net/sword/package/SimpleZip'))
    if(is.null(xml))
		return(NULL)
	if(browser==FALSE)
		.dvParseAtom(xml)
}
