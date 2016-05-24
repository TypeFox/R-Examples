library(XML)

# Doesn't release the doc. Still has one reference.
xmlSize(xmlRoot(xmlParse("~/Data/GovTrack/112/bills/h1.xml")))
gc()


doc = xmlParse("~/Data/GovTrack/112/bills/h1.xml")
xmlSize(xmlRoot(doc))
rm(doc)
gc()


doc = xmlParse("~/Data/GovTrack/112/bills/h1.xml")
ti = getNodeSet(doc, "//titles")
length(ti)
rm(ti, doc)
gc()

doc = xmlParse("~/Data/GovTrack/112/bills/h1.xml")
rm(doc)
gc()
