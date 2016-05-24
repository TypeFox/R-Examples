get.feature.names <-  function(object, gff.file, chr){

chr <- as.character(chr)


 region      <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start       <- region[1]
 end         <- region[2]
 gff.table   <- read.table(gff.file,sep="\t",colClasses=c("NULL","NULL","character","numeric","numeric",rep("NULL",3),"character"),
                           skip = start - 1, nrows = end - start + 1)

 region.pos <- strsplit(object@region.names," - ")
 # region.pos is a list of character vectors

 #print(head(gff.table))
 #
 feature.names <- character(length(object@region.names))

 for(xx in 1:length(feature.names)){

	pos        <- as.numeric(region.pos[[xx]])
	left       <- pos[1]==gff.table[,2]
        right 	   <- pos[2]==gff.table[,3]
        match.pos  <- which(left&right)
        if(length(match.pos)>0){
		#match.pos <- match.pos[1]
                attr <- paste(gff.table[match.pos,1],gff.table[match.pos,4],sep="-->")
	        attr <- paste(attr,collapse=";")
		# print(attr)
        }else{next}
        # print(match.pos)
	# set attribute to feature.names
	feature.names[xx] <- attr
 }

return(feature.names)
}
