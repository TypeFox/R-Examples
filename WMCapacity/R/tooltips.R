# Taken from rattle
.loadTooltips <-function()
{
  if (! packageIsAvailable("XML", "load GUI tooltips"))
  {
    warning("The XML package is not available. Tooltips will not be available.")
    return(FALSE)
  }

  #require(XML, quietly=TRUE)

  filename = getpackagefile("tooltips.xml")
  doc <- xmlTreeParse(filename, useInternalNodes=TRUE)

  for (tt in getNodeSet(doc, "//tooltip"))
  {
    # format the tooltip. blank lines are retained, but other
    # line breaks are ignored.
    tip <- gsub("XoX", "\\\n\\\n",
                gsub("\n *", " ",
                     gsub("\n *\n *", "XoX", xmlValue(tt))))
    wd <- theWidget(xmlGetAttr(tt, 'widget'))
    wd["tooltip-text"] <- Rtxt (tip) #  Space after Rtxt is intentional.
	#print(tip)
    # The MS/Windows RGtk2 is compiled with an older GTK, even
    # though a user might have GTK 2.12.8 installed.  Thus,
    # setTooltipText is not avilable and so we use the above setting
    # of the property rather than using the function.
    
    # theWidget(xmlGetAttr(tt, 'widget'))$setTooltipText(xmlValue(tt))
  }
}

# Taken from rattle
.populateTextviews <- function()
{
  # Reset all text views to default content, as when Rattle
  # starts up or the user has selected New Project, or has loaded a
  # new dataset. The text for the texviews come from textviews.xml.
    
  if (! packageIsAvailable("XML", "load textview texts"))
  {
    warning("The XML package is not available. Textview texts will not be available.")
    return(FALSE)
  }

  #require(XML, quietly=TRUE)
  filename = getpackagefile("textviews.xml")

  doc <- xmlTreeParse(filename, useInternalNodes=TRUE)

  sapply(getNodeSet(doc, "//textview"),
           function(tt)
           {
             wd <- xmlGetAttr(tt, 'widget')
             tv = theWidget(wd)
			 buffer = gtkTextViewGetBuffer(tv)
			 gtkTextBufferSetText(buffer,Rtxt(xmlValue(tt)))
           })

	invisible()
}



