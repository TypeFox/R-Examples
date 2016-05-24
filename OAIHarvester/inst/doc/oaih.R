### R code from vignette source 'oaih.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: oaih.Rnw:193-195
###################################################
library("OAIHarvester")
baseurl <- "http://epub.wu.ac.at/cgi/oai2"


###################################################
### code chunk number 2: oaih.Rnw:199-202
###################################################
require("XML")
options(warnPartialMatchArgs = FALSE)
options(width = 80)


###################################################
### code chunk number 3: oaih.Rnw:207-209
###################################################
x <- oaih_identify(baseurl)
rbind(x, deparse.level = 0)


###################################################
### code chunk number 4: oaih.Rnw:215-216
###################################################
sapply(x$description, xmlName)


###################################################
### code chunk number 5: oaih.Rnw:220-221
###################################################
oaih_transform(x$description[[1L]])


###################################################
### code chunk number 6: oaih.Rnw:226-229
###################################################
oaih_list_metadata_formats(baseurl)
sets <- oaih_list_sets(baseurl)
rbind(head(sets, 3L), tail(sets, 3L))


###################################################
### code chunk number 7: oaih.Rnw:235-236
###################################################
spec <- unlist(sets[sets[, "setName"] == "Type = Thesis", "setSpec"])


###################################################
### code chunk number 8: oaih.Rnw:239-240
###################################################
x <- oaih_list_records(baseurl, set = spec)


###################################################
### code chunk number 9: oaih.Rnw:244-246
###################################################
dim(x)
colnames(x)


###################################################
### code chunk number 10: oaih.Rnw:251-254
###################################################
m <- x[, "metadata"]
m <- oaih_transform(m[sapply(m, length) > 0L])
dim(m)


###################################################
### code chunk number 11: oaih.Rnw:257-258
###################################################
colnames(m)


###################################################
### code chunk number 12: oaih.Rnw:310-311
###################################################
m[c(1L, 6L, 7L), "subject"]


###################################################
### code chunk number 13: oaih.Rnw:315-320
###################################################
sep <- "[[:space:]]*/[[:space:]]*"
keywords_by_thesis <-
    strsplit(unlist(lapply(m[, "subject"],  paste, collapse = " / ")),
             sep)
keywords <- unlist(keywords_by_thesis)


###################################################
### code chunk number 14: oaih.Rnw:324-326
###################################################
counts <- table(keywords)
table(counts)


###################################################
### code chunk number 15: oaih.Rnw:330-331
###################################################
sort(counts[counts >= 3L], decreasing = TRUE)


###################################################
### code chunk number 16: oaih.Rnw:335-336
###################################################
counts["R"]


###################################################
### code chunk number 17: oaih.Rnw:339-342
###################################################
lapply(m[sapply(keywords_by_thesis, function(kw) any(kw == "R")),
         c("title", "creator")],
       strwrap)


###################################################
### code chunk number 18: oaih.Rnw:347-349
###################################################
m[grep("^Feinerer", unlist(m[, "creator"])),
  c("title", "creator", "subject")]


