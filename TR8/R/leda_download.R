##'  Allows the user to retrieve the data files from the LEDA
##'  Traitbase website, merge them in a single \code{R} dataset and store
##'  the result in a local file; this file could be then used whenever the
##'  \code{tr8()} function is used in order to speed up the process of
##'  retrieving traits data.
##'
##'   The function uses a GUI created via the \code{gWidgets} package, to
##'   let the user select a folder where the datasets has to be stored.
##' @title A utility to download a local copy of the LEDA data files.
##' @param directory    is the directory where the downloaded data will be
##' stored (in order  to be used in future R sessions); default is NULL.
##' @return   The function save a local copy of LEDA data in a file called
##' \code{leda_database.Rda}
##' @references Kleyer, M., Bekker, R.M., Knevel, I.C., Bakker, J.P,
##' Thompson, K., Sonnenschein, M., Poschlod, P., Van
##' Groenendael, J.M., Klimes, L., Klimesova, J., Klotz, S.,
##' Rusch, G.M., Hermy, M., Adriaens, D., Boedeltje, G.,
##' Bossuyt, B., Dannemann, A., Endels, P., Götzenberger, L.,
##' Hodgson, J.G., Jackel, A-K., Kühn, I., Kunzmann, D.,
##' Ozinga, W.A., Römermann, C., Stadler, M., Schlegelmilch,
##' J., Steendam, H.J., Tackenberg, O., Wilmann, B.,
##' Corneliss
##' n, J.H.C., Eriksson, O., Garnier, E., Peco, B.
##' (2008): The LEDA Traitbase: A database of life-history
##' traits of Northwest European flora. Journal of Ecology 96:1266-1274. \samp{  http://www.leda-traitbase.org/LEDAportal/data_files.jsp}
##' @author   Gionata Bocci <boccigionata@@gmail.com>
leda_download_to_local_directory<-function(directory){
    ## gmessage(title="","Downloading LEDA files is a time-consuming activity, \nthus you are suggested to download the dataset once\nand store them in a local directory.\n\nYou will now be asked to choose such a directory.")
    ##directory<-gfile(type="selectdir")
    ## load the list containing the data (url, names, etc..) of the txt files
    ## data(leda_lookup)
    ## convert it to a dframe
    env<-new.env(parent = parent.frame())
    data(leda_lookup,envir=env)
    leda_lookup<-get("leda_lookup",envir=env)
    DF<-ldply(leda_lookup)
    ## build the first dataframe
    first<-leda_download(url=DF[1,2],skip_row=DF[1,3],column=DF[1,4],out_name=DF[1,5])
    rearranged<-first
    ## download all the other txt files and merge each one to the
    ## previous one
    for(i in 2:nrow(DF)){
        temp<-leda_download(url=DF[i,2],skip_row=DF[i,3],column=DF[i,4],out_name=DF[i,5])
        rearranged<-merge(rearranged,temp,by.x=0,by.y=0,all.x=TRUE,all.y=TRUE)
        row.names(rearranged)<-rearranged$Row.names
        rearranged<-rearranged[,2:ncol(rearranged)]
    }
    ## save the compleate dataset in a file called "leda_database.Rda"
    ## in the directory chosen by the user
    remove(list=c("leda_lookup"), envir = env)    
    save(file=file.path(directory,"leda_database.Rda"),rearranged)
}
