#función leer.mpa, lee un archivo de texto de palabras clave según el formato indicado
#parámetros: el nombre del archivo y su extensión txt

leer.mpa<-function (textfile = "",encoding="unknown")
{
    if (textfile == "")
        textfile <- file.choose(new = FALSE)
    if (!file.exists(textfile))
        stop(paste("file", textfile, "not found"))
    mpa.file.connector <- file(textfile, "rt")
    mpa.text.lines <- readLines(mpa.file.connector, -1,encoding=encoding)
    close(mpa.file.connector)
    mpa.text.lines <- paste(mpa.text.lines, sep = "")
    return(mpa.text.lines)
}
