clean.cdf=function(where=c("memory","disk","both"),path=getwd())
{
    where=match.arg(where)
    if(where=='memory'){
        return(rm('.pi0cdfp',envir=globalenv()))
    }else if(where=='disk'){
        return(invisible(!file.remove(file.path(dirname(path),basename(path),'.pi0cdfp.RData'))))
    }else{
        if(exists('.pi0cdfp',envir=globalenv()))rm('.pi0cdfp',envir=globalenv())
        return(invisible(!file.remove(file.path(dirname(path),basename(path),'.pi0cdfp.RData'))))
    }
}

save.cdf=function(path=getwd())
{
    if(exists('.pi0cdfp',envir=globalenv())){
        save(list='.pi0cdfp',envir=globalenv(),
            file=file.path(dirname(path),basename(path),'.pi0cdfp.RData')
            )
        return(invisible(TRUE))
    }else {
        warning(".pi0cdfp does not exist in the Global Environment!")
        return(invisible(FALSE))
    }
}

load.cdf=function(path=getwd()){
    try(load(file.path(dirname(path),basename(path),'.pi0cdfp.RData'),envir=globalenv()))
}
