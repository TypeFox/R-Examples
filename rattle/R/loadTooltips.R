loadTooltips <- function()
{
  if (! packageIsAvailable("XML", "load GUI tooltips"))
  {
    warning("The XML package is not available. Tooltips will not be available.")
    return(FALSE)
  }

  result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                silent=TRUE)
  if (inherits(result, "try-error"))
    doc <- XML::xmlTreeParse("tooltips.xml", useInternalNodes=TRUE)
  else
    doc <- XML::xmlTreeParse(file.path(etc, "tooltips.xml"), useInternalNodes=TRUE)

  for (tt in XML::getNodeSet(doc, "//tooltip"))
  {
    # 100110 format the tooltip. blank lines are retained, but other
    # line breaks are ignored.

    tip <- gsub("XoX", "\\\n\\\n",
                gsub("\n *", " ",
                     gsub("\n *\n *", "XoX", XML::xmlValue(tt))))
    wd <- theWidget(XML::xmlGetAttr(tt, 'widget'))
    wd["tooltip-text"] <- Rtxt (tip) # 100408 Space after Rtxt is intentional.

  }
}
