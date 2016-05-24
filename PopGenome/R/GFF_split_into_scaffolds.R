GFF_split_into_scaffolds <- function(GFF.file,output.folder){


dir.create(output.folder)

# Output Folder
if (.Platform$OS.type == "unix") {
 output.folder <- paste(output.folder,"/", sep="")
}else{
output.folder <- paste(output.folder,"\\", sep="")
}

# This Function can be used as well
out <- .Call("split_VCF_scaffolds",GFF.file,output.folder)

return(TRUE)
}
