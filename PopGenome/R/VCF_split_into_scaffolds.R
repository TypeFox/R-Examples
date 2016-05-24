VCF_split_into_scaffolds <- function(VCF.file,output.folder){


dir.create(output.folder)

# Output Folder
if (.Platform$OS.type == "unix") {
 output.folder <- paste(output.folder,"/", sep="")
}else{
output.folder <- paste(output.folder,"\\", sep="")
}


out <- .Call("split_VCF_scaffolds",VCF.file,output.folder)

return(TRUE)
}
