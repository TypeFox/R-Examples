get.ext <- function(file.name, file.has.ext=TRUE) {
# Returns the filename extension, if any, and filename without extension.
# file.name: the string representing one filename.
# file.has.ext: boolean stating whether or not extension is used.
#    If TRUE, then file.name is assumed to be of the format <anything>.<ext>
#    If FALSE, then file.name is assumed to be <anything>
#    
# Returns:
# - ans$part1 - the part of file.name without extension
# - ans$ext - the filename extension of file.name;
#    if file.has.ext=TRUE then it will be the extension with dot: (ex ".txt", ".ped")
#    if file.has.ext=FALSE then it will be an empty string ""
#

name.part1 <- file.name
name.ext <- ""
if(file.has.ext == TRUE){
        name.shred <- unlist(strsplit(file.name, ".", fixed=TRUE))
        name.part1 <- paste(name.shred[1:(length(name.shred)-1)], collapse=".")
        name.ext <- paste(".", name.shred[length(name.shred)], sep="")
}

return(list(part1=name.part1, ext=name.ext))

}
