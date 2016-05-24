read_data <-
function(filepath, type=c("csv","delim","folder"), index=NULL, ...) {    
    if (type=="csv") { data <- read.csv(filepath,...) # CSV FILES
    } else if (type=="delim") { data <- read.delim(filepath,...) # TAB-DELIMITED FILES
	} else if (type=="folder") { 
        if (is.null(index)) stop("Must supply an index if using the folder option.")
        
        labels <- read.csv(index,header=FALSE)
        files <- list.files(path=filepath,full.names=TRUE)
        
        frame <- c()
        for (file in labels[,1]) {
            filename <- NULL
            for (file2 in files) {
                if (basename(file2) == file) {
                    filename <- file2
                    break
                }
            }
            
            if (is.null(filename)) stop("Could not corresponding file from index file in folder.")
            
            lines <- readLines(filename)
            text <- paste(lines,collapse="\n")
            frame <- append(frame,text)
        }
        
        if (nrow(labels) == length(files)) {
            data <- data.frame(Text.Data=frame,Labels=labels[,2])
        } else if (nrow(labels) < length(files)) {
            diff <- length(files)-nrow(labels)
            fill <- as.data.frame(rep(NA,diff))
            
            raw_labels <- as.data.frame(labels[,2])
            colnames(fill) <- colnames(raw_labels)
            print(raw_labels)
            print(fill)
            labels_fixed <- rbind(raw_labels,fill)
            print(labels_fixed)
            
            data <- data.frame(Text.Data=frame,Labels=labels_fixed)
        } else {
            stop("There are more labels than documents in the index file.")
        }
        
    }
	
    return(data)
}