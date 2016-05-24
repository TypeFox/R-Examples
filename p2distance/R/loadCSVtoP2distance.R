loadCSVtoP2distance <- function (path, header=TRUE, sep="\t", dec=".", quote="\"", na.strings="NA", fileEncoding = "", 
    encoding = "unknown") 
{
	matriz <- read.table(path, header = header, sep = sep, dec = dec, quote = quote, na.strings = na.strings, fileEncoding = fileEncoding, encoding = encoding) 

	# Coger los nombres de los municipios 
	nombres <- matriz[,1]
	# Eliminar la columna de municipios y ponerlos como nombres de filas. 
	matriz.datos <- matriz[-1]
	row.names(matriz.datos) <- nombres
	
	return(as.matrix(matriz.datos))
}
