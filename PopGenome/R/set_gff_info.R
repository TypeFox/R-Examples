set_gff_info <- function(object,gff.file, chr){

# reading in the GFF file

my_pos      <- .Call("find_lines_GFF_Human2", gff.file, chr)
gff.table   <- read.table(gff.file,sep="\t",colClasses=c(rep("character",3),rep("integer",2),rep("character",2),"character","NULL"),
                             skip = my_pos[1], nrows = my_pos[2] - my_pos[1]
                         )

#return(gff.table)

change      <- object@region.data
print("1")
features    <- as.character(gff.table[,3])
exon.ids    <- which(features=="exon")
gene.ids    <- which(features=="gene")
coding.ids  <- which(features=="CDS")
print(coding.ids)
print("2")
ffc <- as.numeric(gff.table[coding.ids,4:5])
ffg <- as.numeric(gff.table[gene.ids,4:5])
ffe <- as.numeric(gff.table[exon.ids,4:5])
print("3")
change@Coding.matrix[[1]] <- ff(ffc,dim=c(length(coding.ids),2))
change@Gene.matrix[[1]]   <- ff(ffg,dim=c(length(gene.ids),2))
change@Exon.matrix[[1]]   <- ff(ffe,dim=c(length(exon.ids),2))
print("4")
object@region.data <- change
print("5")
return(object)


}
