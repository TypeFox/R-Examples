######################################################################
# Rjpstatdb: Interface to ``Gateway to Advanced and User-friendly
#            Statistics Service''
#
# Author: Kiwamu Ishikura
######################################################################

### required packages
library(methods)
library(RCurl)
library(XML)


### Common parameters
.DBURL <- "http://statdb.nstac.go.jp/api/1.0b/app/"
.appId <- "03f6a0a1f01ec96fdaa8920c51bebbafcab71927"
#.lang <- ifelse(length(grep("ja_JP",
#                            unlist(strsplit(Sys.getlocale(), "/")))) > 0,
#                "J", "E")
.lang <- "J"


### Definition of class
setClass(Class = "jpstat",
         representation(data = "list", id = "character",
                        stat.name = "character", gov = "character",
                        statistics.name = "character",
                        title = "character", survey.date = "character"),
         prototype = list(
             data = list(data.frame(0)), id = "", stat.name = "",
             gov = "", statistics.name = "", title = "",
             survey.date = ""))

print.jpstat <- function(x, ...) {
    cat("ID: ", x@id, "\n")
    cat("Stat name: ", x@stat.name, "\n")
    cat("Government: ", x@gov, "\n")
    cat("Statistics name: ", x@statistics.name, "\n")
    cat("Title: ", x@title, "\n")
    cat("Survey date: ", x@survey.date, "\n")
    cat("Tables:\n")
    for (i in seq(x@data)) {
        cat("  Name: ", names(x@data)[i], "\n")
        print(head(x@data[[i]]), ...)
        cat("  Dimension: ", dim(x@data[[i]]), "\n")
    }
    return(invisible(NULL))
}

setMethod("show", "jpstat",
          function(object) print.jpstat(x = object))



### Get statistics list
getStatsList <- function(searchWord = "", surveyYears = "",
                         openYears = "", statsField = NULL,
                         statsCode = NULL, searchKind = 1,
                         statsNameList = "") {
    ## set locale temporarily
    locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_ALL", "ja_JP.UTF-8")

    gf <- match.call(expand.dots = FALSE)
    m <- match(c("searchWord", "surveyYears", "openYears", "statsField",
                 "statsCode", "searchKind", "statsNameList"),
               names(gf), 0L)
    gf <- gf[c(1L, 1L, m)]
    gf[[1L]] <- as.name("getForm")
    names(gf)[2L] <- "uri"
    gf[[2L]] <- paste(.DBURL, "getStatsList", sep = "")
    gf$appId <- .appId
    gf$lang <- .lang
    root <- xmlRoot(xmlTreeParse(eval(gf)))

    status <- as.numeric(xmlValue(root[[1L]]["STATUS"][[1L]]))
    if (status >= 100)
        stop(xmlValue(root[[1L]]["ERROR_MSG"][[1L]]))
    else if (status == 1) {
        cat(xmlValue(root[[1L]]["ERROR_MSG"][[1L]]), "\n")
        return(invisible(NULL))
    }

    ResList <- root[[3L]][names(root[[3L]]) == "LIST_INF"]
    res <- data.frame(
        Data_Set_ID = sapply(ResList, xmlGetAttr, "id"),
        Stat_Name = sapply(ResList,
                           function(x)
                               paste(substr(xmlValue(x["STATISTICS_NAME"][[1L]]),
                                            1, 10), "...", sep="")),
        Stat_Code = sapply(ResList,
                           function(x) xmlGetAttr(x["STAT_NAME"][[1L]],
                                                  "code")),
        Org = sapply(ResList,
                     function(x) xmlValue(x["GOV_ORG"][[1L]])),
        Survey = sapply(ResList,
                        function(x) xmlValue(x["SURVEY_DATE"][[1L]])),
        Open = sapply(ResList,
                      function(x) xmlValue(x["OPEN_DATE"][[1L]])))

    ## reset locale
    Sys.setlocale("LC_ALL", locale)

    res
}


### Get data
getStatsData <- function(statsDataId = NULL, dataSetId = NULL,
                         limit = NULL, lvTab = "", cdTab = NULL,
                         lvTime = "", cdTime = NULL,
                         lvArea = "", cdArea = NULL) {
    if (is.null(statsDataId) & is.null(dataSetId))
        stop("Either statsDataId or dataSetId should be specified.")

    ## set locale temporarily
    locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_ALL", "ja_JP.UTF-8")

    gf <- match.call(expand.dots = FALSE)
    m <- match(c("statsDataId", "dataSetId", "limit", "lvTab",
                 "cdTab", "lvTime", "cdTime", "lvArea", "cdArea"),
               names(gf), 0L)
    gf <- gf[c(1L, 1L, m)]
    gf[[1L]] <- as.name("getForm")
    names(gf)[2L] <- "uri"
    gf[[2L]] <- paste(.DBURL, "getStatsData", sep = "")

    if (length(gf$cdTab) > 1)
        gf$cdTab <- paste(gf$cdTab, collapse = ",")
    if (length(gf$cdTime) > 1)
        gf$cdTime <- paste(gf$cdTime, collapse = ",")
    if (length(gf$cdArea) > 1)
        gf$cdArea <- paste(gf$cdArea, collapse = ",")

    gf$appId <- .appId
    gf$lang <- .lang
    root <- xmlRoot(xmlTreeParse(eval(gf)))

    status <- as.numeric(xmlValue(root[[1L]]["STATUS"][[1L]]))
    if (status >= 100)
        stop(xmlValue(root[[1L]]["ERROR_MSG"][[1L]]))
    else if (status == 1) {
        cat(xmlValue(root[[1L]]["ERROR_MSG"][[1L]]), "\n")
        return(invisible(NULL))
    }

    ## Tables
    tables <- getNodeSet(root[[3L]][[2L]],
                         "//CLASS_INF//CLASS_OBJ[@id='tab']")[[1L]]
    name.tables <- xmlSApply(tables, xmlGetAttr, "name")
    names(name.tables) <- xmlSApply(tables, xmlGetAttr, "code")

    ## Categories
    cats <- getNodeSet(root[[3L]][[2L]],
                       "//CLASS_INF//CLASS_OBJ[contains(@id, 'cat')]")
    name.cats <- level.cats <- vector(length(cats), mode = "list")
    for (i in seq(cats)) {
        nc <- xmlSApply(cats[[i]], xmlGetAttr, "name")
        lc <- xmlSApply(cats[[i]], xmlGetAttr, "level")
        names(nc) <- names(lc) <-
            xmlSApply(cats[[i]], xmlGetAttr, "code")
        name.cats[[i]] <- nc
        level.cats[[i]] <- lc
    }
    names(name.cats) <- sapply(cats, xmlGetAttr, "id")
    names(level.cats) <- sapply(cats, xmlGetAttr, "id")

    ## Area
    areas <- getNodeSet(root[[3L]][[2L]],
                        "//CLASS_INF//CLASS_OBJ[@id='area']")[[1L]]
    name.areas <- xmlSApply(areas, xmlGetAttr, "name")
    level.areas <- xmlSApply(areas, xmlGetAttr, "level")
    names(name.areas) <- names(level.areas) <-
        xmlSApply(areas, xmlGetAttr, "code")

    ## Time
    times <- getNodeSet(root[[3L]][[2L]],
                       "//CLASS_INF//CLASS_OBJ[@id='time']")[[1L]]
    name.times <- xmlSApply(times, xmlGetAttr, "name")
    level.times <- xmlSApply(times, xmlGetAttr, "level")
    names(name.times) <- names(level.times) <-
        xmlSApply(times, xmlGetAttr, "code")

    res.data <- vector(length(tables), mode = "list")
    for (i in 1:length(tables)) {
        data <- getNodeSet(root[[3L]][[3L]],
                           paste("//DATA_INF//VALUE[@tab='",
                                 names(name.tables[i]),
                                 "']", sep = ""))
        res.i <- data.frame(
            value = as.numeric(sapply(data, xmlValue)),
            area.code = sapply(data, xmlGetAttr, "area"),
            time.code = sapply(data, xmlGetAttr, "time"))
        res.i <- transform(res.i,
                           area = name.areas[res.i$area.code],
                           area.level = level.areas[res.i$area.code],
                           time = name.times[res.i$time.code],
                           time.level = level.times[res.i$time.code])
        for (nc in names(name.cats)) {
           nc.full <- sapply(data, xmlGetAttr, nc)
           res.i[paste(nc, ".code", sep="")] <- nc.full
           res.i[nc] <- name.cats[nc][[1L]][nc.full]
           res.i[paste(nc, ".level", sep="")] <- level.cats[nc][[1L]][nc.full]
        }

        res.data[[i]] <- res.i
    }
    names(res.data) <- name.tables

    res <- new(
        "jpstat",
        data = res.data,
        id = xmlValue(getNodeSet(root[[2L]],
                      "//PARAMETER//STATS_DATA_ID")[[1L]]),
        stat.name = xmlValue(getNodeSet(root[[3L]][[1L]],
                             "//TABLE_INF//STAT_NAME")[[1L]]),
        gov = xmlValue(getNodeSet(root[[3L]][[1L]],
                                  "//TABLE_INF//GOV_ORG")[[1L]]),
        statistics.name = xmlValue(getNodeSet(root[[3L]][[1L]],
                                              "//TABLE_INF//STATISTICS_NAME")[[1L]]),
        title = xmlValue(getNodeSet(root[[3L]][[1L]],
                                    "//TABLE_INF//TITLE")[[1L]]),
        survey.date = xmlValue(getNodeSet(root[[3L]][[1L]],
                                          "//TABLE_INF//SURVEY_DATE")[[1L]]))

    ## reset locale
    Sys.setlocale("LC_ALL", locale)

    res
}
