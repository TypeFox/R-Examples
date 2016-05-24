readAssayTable <-
    function (file = "AssayTable.txt",
              fun = NULL,
              log = FALSE,
              rows = "",
              columns = "",
              replicateNames = c("Repetition", "Replicate", "Block",
                  "Plate", "Animal", "Group", "Litter"),
              doseNames = c("Dilution", "Dose"),
              sampleNames = c("Smp", "Sample"),
              otherNames = c("Outer", "Inner", "Assay"),
              responseName = "Response",
              combinedTreatment = FALSE,
              echoTitle = TRUE)
{
    options(SweaveHooks = list(fig = function()
        par(bg = "white", fg = "black")))
    color       <- "black"
    font        <- "\\\\rm"
    design      <- "";
    Data        <- NULL
    RawTable    <- matrix(nrow = 0, ncol = 0)
    useFill     <- FALSE
    Description <- NULL
    Name        <- "";
    firstColumn <- NA;
    ColumnNames <- NULL;
    try(Data <- read.table(file, header = TRUE))
    if (is.null(Data)) {
        warning("Dataformat is not correct\n")
        cat("Warning: Dataformat is not correct\n")
        useFill  <- TRUE
        options(SweaveHooks = list(fig = function()
            par(bg = "red", fg = "blue")))
        try(Data <- read.table(file, header = TRUE, fill = TRUE))
        color <- "red"
        font <- "\\\\it"
        if (!is.null(Data)) {
            warning(paste("Datafile with variable length of line read",
                          "- You should correct the datafile!"),
                    call. = FALSE)
            cat("Warning: Datafile with variable length of line read\n")
            cat("Warning: You should correct the datafile!\n")
        }
    }
    if (is.null(Data)) {
        warning("Warning: Giving up!", call. = FALSE)
        cat("Warning: Giving up!\n")
    }
    else {
        Header   <- dimnames(Data)[[2]]
        Reserved <- c(replicateNames, sampleNames, doseNames,
                      otherNames, "Row", "Column", responseName)
        Reserved <- responseName
        matching <- match(Reserved, Header)
        if (length(Header[matching][!is.na(matching)]) > 0) {
            readingTable      <- FALSE
            ## warning("You are readig a data.frame, not a table!")
            Plates            <- substr(unlist(Data[1]), 1, 1) != "#"
            ## RawTable       <- toTable(Data[Plates,])
            RawTable          <- Data[Plates,]
            design            <- "lsd"; ## !!!
            ## Table          <- Data[Plates,]
            mcm <- function(x, table = Header) {
                m <- charmatch(tolower(x), tolower(table))
                m <- m[!is.na(m)]
                return(unique(sort(m)))
            }
            indexHeader <- list(
                Response  = mcm(responseName),
                Label     = mcm(c("Label")),
                Inner     = mcm(c("Inner")),
                Sample    = mcm(sampleNames),
                Dose      = mcm(doseNames),
                Replicate = mcm(replicateNames),
                Row       = mcm(c("Row")),
                Column    = mcm(c("Column")),
                Outer     = mcm(c("Outer")),
                Assay     = mcm(c("Assay")))
            ## Patch for extra factors is just theses
            ## four plus the below two lines:
            restNames <- otherNames[is.na(match(otherNames,
                                                c("Outer", "Inner", "Assay")))]
            if (length(restNames) > 0)
                for (i in 1:length(restNames)) {
                    indexHeader <- append(indexHeader, mcm(restNames[i]))
                    names(indexHeader)[length(indexHeader)] <- restNames[i]
                }
            missingFactors <- lapply(indexHeader[4:6], length) == 0
            Description <- list(indexHeader = indexHeader,
                                missingFactors = missingFactors)
        } else {
            Factors           <- c(sampleNames, doseNames, replicateNames,
                                   otherNames, "Label")
            Factors           <- c(Factors, tolower(Factors))
            Factors           <- c(Factors, substr(Factors, 1, 3),
                                   substr(Factors, 1, 1))
            Factors           <- unique(sort(Factors))
            mtch              <- match(Factors, Header, nomatch = 0)
            oneCharMtch       <- (mtch > 1) & (nchar(Factors) == 1) &
                (Factors != "R") & (Factors != "D") & (Factors != "S")
            ## That is, 'R', 'D' and 'S' is still 'ok'.
            oneCharMtch       <- oneCharMtch |
                ((mtch > 2) & (nchar(Factors) == 1))
            ## Now 'R', 'D' and 'S' are also 'deleted'
            ## for columns after number 2.
            mtch[oneCharMtch] <- 0
            indexMatch        <- mtch[mtch > 0]
            HeaderMatch       <- Header[mtch]
            FactorsMatch      <- Factors[mtch > 0]
            Name              <- Header[1]
            ## It would be nice to not have the e.g. 'Replicate-' part:
            Plates            <- substr(unlist(Data[min(indexMatch)]), 1,
                                        nchar(Name) + 1) ==
                                            paste(Name, "-", sep = "")
            readingTable      <- TRUE
            RawTable          <- Data[Plates, -indexMatch]
            TableFirst        <- Data[Plates, 1]
            firstColumn       <- substr(TableFirst, nchar(Name)+2,
                                        nchar(as.character(TableFirst)))
            TableHead         <- as.character(Data[Plates, dim(Data)[2]])
            ColumnNames       <- as.character(Data[Plates, max(indexMatch)])
            if (length(indexMatch[FactorsMatch != "Label"]) > 1)
                RowNames   <- apply(Data[Plates,
                                         indexMatch[FactorsMatch !=
                                                        "Label"]], 1,
                                    FUN = function(i)
                                        paste0(i, collapse = ":"))
            else
                RowNames   <- Data[Plates,
                                   indexMatch[FactorsMatch != "Label"]]
            dimnames(RawTable)[[1]]   <- RowNames
            TableMatrix       <- matrix(as.numeric(as.matrix(RawTable)),
                                        ncol = dim(RawTable)[2])
            dimnames(TableMatrix)[[1]] <- RowNames
            dimnames(TableMatrix)[[2]] <- dimnames(RawTable)[[2]]
            invalidCells      <- is.na(TableMatrix) & !is.na(RawTable)
            if (any(invalidCells)) {
                invalidCellsT <- ifelse(invalidCells,
                                        RawTable[invalidCells], NA)
                cat("Warning: Invalid cells in the table:\n")
                print(invalidCellsT, na.print = "")
                cat("\n")
            }
            if (length(RowNames) == 0)
                warning("No data (levels of first factor invalid)")
            else
                Description <- .findArrayDescription(
                    Data, Plates, RawTable, HeaderMatch,
                    FactorsMatch, indexMatch, rows, columns,
                    replicateNames = replicateNames, doseNames = doseNames,
                    sampleNames = sampleNames, otherNames = otherNames,
                    combinedTreatment = combinedTreatment)
        }
    }
    if (!is.null(Data)) {
        result <- .insertTable(
            file, fun, log, rows, columns, echoTitle,
            useFill, Data, readingTable, Plates, Name, firstColumn,
            ColumnNames, RawTable, Description, design, color, font)
    }
    return(result)
}
