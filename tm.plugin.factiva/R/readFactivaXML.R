readFactivaXML <- readXML(
    spec = list(Author = list("function", function(node)
                toupper(gsub("^\\s+|\\s+$", "",
                             gsub("\n|\\s+", " ",
                                  sapply(getNodeSet(node, "/article/byline"), xmlValue))))),
    content = list("function", function(node)
                   c(sapply(getNodeSet(node, "/article/headline"), xmlValue),
                         sapply(getNodeSet(node, "/article/leadParagraph"), xmlValue),
                         sapply(getNodeSet(node, "/article/tailParagraphs/paragraph"), xmlValue))),
    datetimestamp = list("function", function(node)
                         strptime(sapply(getNodeSet(node, "/article/publicationDate/date"), xmlValue),
                                  format="%Y-%m-%d")),
    heading = list("node", "/article/headline"),
    id = list("function", function(node) {
              str <- gsub("^distdoc:archive/ArchiveDoc::Article/", "",
                           sapply(getNodeSet(node, "/article/reference"), xmlValue))
              # Extract useful information: origin, date, and two last characters to avoid collisions
              m <- regmatches(str, regexec("^([A-Za-z]+)0*[1-9][0-9]([0-9][0-9][0-3][0-9][0-3][0-9]).*([A-Za-z0-9]{2})$", str))[[1]]
              # If extraction failed for some reason, make sure we return a unique identifier
              if(length(m) == 4)
                  paste(toupper(m[2]), "-", m[3], "-", m[4], sep="")
              else
                  paste(sample(LETTERS, 10), collapse="")
    }),
    origin = list("node", "/article/sourceName"),
    language = list("function", function(node)
                    tolower(sapply(getNodeSet(node, "/article/baseLanguage"), xmlValue))),
    edition = list("node", "/article/edition"),
    section = list("node", "/article/sectionName"),
    subject = list("node", "/article/newsSubject/name"),
    coverage = list("node", "/article/region/name"),
    company = list("node", "/article/company/name"),
    industry = list("node", "/article/industry/name"),
    infocode = list("node", "/article/descField[@code!='ipd']"),
    infodesc = list("function", function(node) {
                    str <- sapply(getNodeSet(node, "/article/descField[@code='ipd']"), xmlValue)
                    if(length(str) > 0)
                        strsplit(str, "( +\\| +| +-+ +| +--+|--+ +|\\._)")[[1]]
                    else
                        character(0)
    }),
    wordcount = list("function", function(node)
                     as.numeric(sapply(getNodeSet(node, "/article/wordCount"), xmlValue))),
    publisher = list("node", "/article/publisherName"),
    rights = list("function", function(node)
                  gsub("^\\s+|\\s+$", "",
                       gsub("\n|\\s+", " ",
                            sapply(getNodeSet(node, "/article/copyright"), xmlValue))))),
    doc = PlainTextDocument())

