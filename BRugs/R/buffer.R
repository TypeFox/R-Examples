buffer <- function(){
    buffer <- file.path(tempdir(), "buffer.txt")
    message(readLines(buffer))
}
