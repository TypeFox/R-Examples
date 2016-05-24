batchimport <- function(raw_metadata, batch_vars, name_batch){
    select_info <- as.character(raw_metadata[,1] )%in%batch_vars
    batch_block <- raw_metadata[select_info,]
    batch_block[,1] <- as.character(batch_block[,1]) 
    batch_block[,2] <- as.character(batch_block[,2]) 
    batch_block[,3] <- as.character(batch_block[,3]) 
    batch_block[1,2] <- paste(as.character(batch_block[1,2]), 
                        batch_block[1,4],sep="*")
    batch_data <- batch_block[,c(1,2)]
    col_names <-  batch_data[,1] 
    batch_data <- as.data.frame(t(batch_data))
    batch_data$batch <- name_batch
    batch_data$crc <- unique(raw_metadata$crc)
    batch_data <- batch_data[-1,]
    names(batch_data) <- c( col_names, "batch","crc")
    names(batch_data) <- gsub(" ","_",(tolower(names(batch_data))))
    batch_data <- rename.vars(batch_data, from="protocoldevelopingcompany", 
                        to = "protocolcompany",info=FALSE)
    selec_names <- names(batch_data)[names(batch_data)%nin%c("batch")]
    vars <- c("batch",selec_names)
    batch_data <- batch_data[, vars]
    aux.chron <- chron(as.character(batch_data$date), 
                        format=c(dates="d/m/y"), 
                        out.format=c(dates="m/d/y"))
    batch_data$date <- as.character(aux.chron)
    batch_data
}


