PLANTS_download_to_local_directory<-function(directory){
    env<-new.env(parent=parent.frame())
    data(ref_PLANTS,envir=env)
    ref_PLANTS<-get("ref_PLANTS",envir=env)
    url<-"http://bricol.net/downloads/data/PLANTSdatabase/09-02-02PLANTSdata.csv"
    PLANTS<-read.table(url,fileEncoding="iso-8859-1")
    
    ## remove commas in species names since they conflict with TNRS

    PLANTS$Scientific.Name<-as.character(PLANTS$Scientific.Name)
    PLANTS$Scientific.Name<-gsub(","," ",PLANTS$Scientific.Name)
    PLANTS<-merge(ref_PLANTS,PLANTS,,by.x="Symbol",by.y="Symbol")
    save(file=file.path(directory,"PLANTS.Rda"),PLANTS)
    remove(list=c("ref_PLANTS"), envir = env)    

}
