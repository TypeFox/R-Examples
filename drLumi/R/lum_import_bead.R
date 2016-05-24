lum_import_bead <- function(x){    
    n <- nchar(basename(x))
    extension <- substring(basename(x), first = n-3, last=n)
    newname <- basename(x)
    newx <- x
    if(extension==".zip"){
        newname <- substring(basename(x), first = 1, last =n-4)
        newx <- dirname(x)
        if(file.exists(file.path(newx, newname))){
            stop(paste("'",newname,"'" ,"exists in .zip folder. No unzip."))
        } 
        dir.create(file.path(newx, newname))
        unzip(x, exdir= file.path(newx, newname), 
              overwrite = TRUE, junkpaths = TRUE)
        all_files <- list.files(file.path(newx, newname))
        all_path <- file.path(newx, newname,all_files)
    } else {
      all_files <- list.files(file.path(x))
      all_path <- file.path(newx, all_files)
    }   
    
    if(length(all_files)==0L) stop(paste("No files found in x"))
    
    extension <- substring(all_files, first =nchar(all_files)-2, 
                           last=nchar(all_files))    
    if(!all(extension%in%c("csv","CSV"))) stop("Files must be csv")    
    
    batch_name <- strsplit(newname,"_rcsv")[[1]]    
        
    bead_files_list <- lapply(all_path, function(y)  bead_fun(y, batch_name))
    bead_files <- ldply(bead_files_list)
    ans <- list(bead_files = bead_files, name_batch = batch_name)
    ans
}

