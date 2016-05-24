#' Strips a file extension from a filename.
#' @title remove_file_extension
#' @param filename Character. The input filename.
#' @param extension_delimiter Character. The extension or extension delimiter (default ".") to remove.
#' @author Jonathan A. Greenberg \email{spatial.tools@@estarcion.net}
#' @examples
#' myfilename="my.file.gri"
#' remove_file_extension(myfilename,".")
#' remove_file_extension(myfilename,".file.gri")
#' @export


remove_file_extension=function(filename,extension_delimiter=".")
{
	split_filename=unlist(strsplit(filename,extension_delimiter,fixed=TRUE))
	split_filename_length=length(split_filename)
	if(split_filename_length==1)
	{
		return(split_filename[1])
	} else
	{
		return(paste(as.character(split_filename)[1:(split_filename_length-1)],collapse=extension_delimiter))	
	}
}
