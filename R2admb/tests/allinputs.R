library("R2admb")
type <- "tpl"
if (FALSE) {
    ## grab relevant files from ADMB tree
    findstr <- paste0('find ~/admb/admb-trunk/examples -path "*/autodif/*" -prune -o -name "*.',type,'"')
    cpstr <- "-exec cp {} inputs \\;"
    printstr <- "-print"
    system(paste(findstr,cpstr))
}

check_all <- function(type,skipfiles=NULL,fun) {
    res <- list()
    files <- list.files("inputs",pattern=paste0(".",type))
    filelist <- gsub(paste0("\\.",type,"$"),"",files)
    filelist <- setdiff(filelist,skipfiles)
    for (i in seq_along(filelist)) {
        cat(filelist[[i]],"\n")
        res[[i]] <- fun(file.path("inputs",filelist[[i]]))
    }
}

check_all("tpl",fun=R2admb:::read_tpl)
check_all("par",fun=read_pars)
