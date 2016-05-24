## Create an iframe to represent all icons in a category in Komodo
## iconpicker from a list of URIs
makeIconGallery <- function (flist)
{
	flist <- as.character(flist)[1]
	if (!file.exists(flist))
		stop("'flist' file not found")
	## Read the list
	icns <- readLines(flist)
	## Eliminate empty lines
	icns <- icns[icns != ""]
	if (length(icns) < 1)
		stop("Nothing in the 'flist' file!")
	## Create the iframe
	iframe <- sub("\\.txt$", ".html", flist)
	if (iframe == flist)
		iframe <- paste(flist, "html", sep =".")
	head <- '<html>
<body>
<style>img:hover { border-color: black; }</style>
<style>img { border-color: white; }</style>

'
	tail <- '</body>\n</html>'
	item <- '<img border="1"
	ondblclick="parent.ValidatedPickIcon(<<<uri>>>);"
	onclick="parent.Pick_Icon(\'<<<uri>>>\');"
    src="<<<uri>>>"
    alt="<<<icn>>>"
    style="padding: 4px;"/>
	 
'
	cat(head, file = iframe)
	for (i in 1:length(icns)) {
		uri <- icns[i]
		icn <- basename(uri)
		itm <- gsub("<<<uri>>>", uri, item)
		itm <- sub("<<<icn>>>", icn, itm)
		cat(itm, file = iframe, append = TRUE)
	}
	cat(tail, file = iframe, append = TRUE)
	## Check if the file exists
	return(invisible(file.exists(iframe)))
}
