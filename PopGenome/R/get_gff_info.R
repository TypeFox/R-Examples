get_gff_info <- function(object=FALSE, gff.file, chr, position, feature=FALSE, extract.gene.names = FALSE){

chr <- as.character(chr)

#if(nchar(chr)==1){
#   chr <- c(chr,"z") # for find_lines_GFF_Human
#}else{
#   chr <- strsplit(chr,split="")[[1]]
#}

if(extract.gene.names){

 region      <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start       <- region[1]
 end         <- region[2]
 gff.table   <- read.table(gff.file,sep="\t",colClasses=c("NULL","NULL","character",rep("NULL",5),"character"),
                           skip = start - 1, nrows = end - start + 1)
 ids       <- which(gff.table[,1]=="gene")
 gff.table <- gff.table[ids, ]
 
return(gff.table[,2])

}


if(feature[1]!=FALSE){

 region      <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start       <- region[1]
 end         <- region[2]
 gff.table   <- read.table(gff.file,sep="\t",colClasses=c(rep("NULL",2),"character",rep("integer",2),rep("NULL",4)),
                           skip = start - 1, nrows = end - start + 1)
 featids     <- which(gff.table[,1]==feature)
 poslist     <- apply(gff.table[featids,2:3],1,function(x){return(x[1]:x[2])})

return(poslist)
}

if(is(object)=="logical"){
 chr    <- as.character(chr)
 region <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start  <- region[1]
 end    <- region[2]
 info   <- .Call("get_gff_info_C",start,end,gff.file,position)
return(info)
}

if(is(object)=="GENOME"){

 chr    <- as.character(chr)
 region <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start  <- region[1]
 end    <- region[2]

RET.INFO <- vector("list",length(position))

for(xx in 1:length(position)){

 region <- position[xx]
 region <- object@region.data@biallelic.sites[[region]]
 info   <- sapply(region,function(x){return(.Call("get_gff_info_C",start,end,gff.file,x))})
 names(info)    <- region
 RET.INFO[[xx]] <- info

}

RET.INFO           <- as.matrix(RET.INFO)
rownames(RET.INFO) <- position

return(RET.INFO)
}


}
