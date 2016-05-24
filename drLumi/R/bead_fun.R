bead_fun <- function(x, batch_name){
    bead_name <- strsplit(basename(x),".csv")[[1]]
   # bead_name <- strsplit(bead_name,batch_name)[[1]][2]
   # bead_name <- substring(bead_name,2,999)
    ncol <-  max(count.fields(x, sep=","))
    bead_file <- read.table(x, header = FALSE, sep = ",", 
                            col.names = paste0("V",seq_len(ncol)), fill = TRUE)
    names(bead_file) <- tolower(as.character(as.matrix(bead_file[2,])))
    bead_file <- bead_file[-c(1,2),]
    bead_file$well <- bead_name
    bead_file$batch <-  batch_name 
    vars <- c("batch","well","eventno")
    bead_file$batch_well_eventno <- apply(bead_file[,vars],1,
                            function(x) paste(x, collapse="*") )
    return(bead_file)
}


