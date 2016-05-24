lum_import_fluor <- function(x){         
    ncol <-  max(count.fields(x, sep=","),na.rm=TRUE)
    raw_file <- read.table(x, header = FALSE, sep = ",",
                col.names = paste0("V",seq_len(ncol)), fill = TRUE)
    srow <- as.numeric(rownames(raw_file[raw_file[,1]=="Results",] ))
    raw_metadata <- raw_file[1:(srow-1),]
    name_batch <- as.character(raw_metadata[raw_metadata[,1]=="Batch",2])
    name_batch <- fix_name(name_batch)
    srow <- as.numeric(rownames(raw_file[raw_file[,1]=="DataType:",] ))
    lastblock <- as.numeric(rownames(raw_file[raw_file[,1]=="-- CRC --",] ))
    if(length(lastblock)==0) lastblock <- nrow(raw_file)+1
    snames <- as.character(raw_file[srow,2])
    srow <- c(srow, lastblock)
    i <- 1:(length(srow)-1)
    data_type_blocks <-  lapply(i, 
                function(x) raw_file[(srow[x]+1):(srow[x+1]-1),])
    names(data_type_blocks) <- snames
    dtblock <- lapply(data_type_blocks, function(x) define_data(x)) 
    for(i in 1:length(snames)){ 
        attr(dtblock[[snames[i]]],"type_results_data") <- snames[i] 
    }
    for(i in 1:length(snames)){ 
        attr(dtblock[[snames[i]]],"batch_name") <- name_batch 
    }
    for(i in 1:length(snames)){
        auxdf <- dtblock[[snames[i]]]
        auxdf <- subset(auxdf, auxdf[,1]!="")
        dtblock[[snames[i]]]  <- auxdf
    }
    crc_value <- as.character(raw_file[nrow(raw_file),1])
    if(crc_value=="Value") crc_value <- "not_in_file"
    raw_metadata$crc <- crc_value
    well_vars <- c("Median","Net MFI","Count","Result",
                "Range","% Recovery","Comments")
    scurve_vars <- c("Analysis Coefficients","R^2")
    average_vars <- c("Avg Net MFI","Avg Result","Avg Range",
                "%CV Replicates","Standard Expected Concentration")
    batch_vars <- c("Program","Build","Date","SN","Batch","Version",
                "Operator","ComputerName","Country Code",
                "ProtocolName","ProtocolVersion","ProtocolDescription",
                "ProtocolDevelopingCompany","SampleVolume","DDGate",
                "SampleTimeout","BatchStartTime","BatchStopTime")
    region_vars <- c("Units")
    sample_vars <- c("Dilution Factor")
    well_vars <- well_vars[well_vars%in%names(dtblock)]
    scurve_vars <- scurve_vars[scurve_vars%in%names(dtblock)]
    average_vars <- average_vars[average_vars%in%names(dtblock)]
    region_vars <- region_vars[region_vars%in%names(dtblock)]
    sample_vars <- sample_vars[sample_vars%in%names(dtblock)]
    ans <- list(dtblock = dtblock, 
            raw_metadata = raw_metadata,
            well_vars = well_vars, 
            scurve_vars = scurve_vars, 
            average_vars = average_vars, 
            batch_vars = batch_vars,
            name_batch = name_batch,
            region_vars = region_vars,
            sample_vars = sample_vars)
    ans
}
