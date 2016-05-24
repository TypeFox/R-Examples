split_data_into_GFF_attributes <- function(object, gff.file, chr, attribute){

chr <- as.character(chr)


 region      <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start       <- region[1]
 end         <- region[2]
 gff.table   <- read.table(gff.file,sep="\t",colClasses=c("NULL","NULL","character","numeric","numeric",rep("NULL",3),"character"),
                           skip = start - 1, nrows = end - start + 1)

 #print(head(gff.table))
 Attr        <- gff.table[,4]


 # iterate over attributes to verify length of POSList
 # Found first attribute

 for(xx in 1:length(Attr)){
	
        line <- Attr[xx]
	line <- strsplit(line ,";")[[1]]
        mm   <- grep(attribute, line)
        if(length(mm)!=0){ # found attribute
                 if(gff.table[xx,1]=="chromosome"){next}
        	 id  <- mm[1]
                 Str <- line[id]
		 Str <- strsplit(Str ,",")[[1]]
		 mm  <- grep(attribute, Str)
		 id  <- mm[1]
                 Str <- Str[id]
        print("Found attribute")
        break
        }
 }

 print(Str)
 SAVENAMES    <- rep(NA,length(Attr))
 SAVENAMES[1] <- Str

 # Count unique attributes to verify length of list
 for(xx in 1:length(Attr)){
	
        line <- Attr[xx]
	line <- strsplit(line ,";")[[1]]
        mm   <- grep(attribute, line)
        if(length(mm)!=0){ # found attribute
		 if(gff.table[xx,1]=="chromosome"){next}
        	 id  <- mm[1]         
	         Str2 <- line[id]
		 Str2 <- strsplit(Str2 ,",")[[1]]
		 mm  <- grep(attribute, Str2)
		 id  <- mm[1]
                 Str2 <- Str2[id]
                #if(Str2!=Str){SAVENAMES[xx] <- Str; Str <- Str2}
		 SAVENAMES[xx] <- Str2
        }
 }

 XSAVENAMES <- SAVENAMES
 SAVENAMES  <- SAVENAMES[!is.na(SAVENAMES)]
 SAVENAMES  <- unique(SAVENAMES)
 
 POSList <- vector("list",length(SAVENAMES))
 
 # Generate position of each unique attribute id
 for(xx in 1:length(SAVENAMES)){
    cat(SAVENAMES[xx],":",xx,"of",length(SAVENAMES),"\n")
    #ids   <- grep(SAVENAMES[xx],Attr)
    ids  <- SAVENAMES[xx] == XSAVENAMES
    ids  <- which(ids)
    pos  <- gff.table[ids,2:3]
    pos  <- unique(pos)
    #print(pos) 
    pos	 <- apply(pos,1,function(x){return(x[1]:x[2])})
    pos  <- unlist(pos)
    pos  <- unique(as.numeric(pos))
    pos  <- sort(pos)
    POSList[[xx]] <- pos
 }

print("Split the data into Attributes")
object               <- splitting.data(object, positions= POSList, type=2)
object@feature.names <- SAVENAMES 

 return(object)
}

