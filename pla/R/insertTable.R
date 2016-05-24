.insertTable <-
    function (file = "AssayTable.txt",
              fun = NULL,
              log = FALSE,
              rows = "Replicates",
              columns = "Samples & Dilutions",
              echoTitle = TRUE,
              useFill,
              Data,
              readingTable,
              Plates,
              firstName,
              firstColumn,
              columnNames,
              TableChar,
              description,
              design,
              color,
              font
              )
{
    myLine <- function(char = "=", length = 80, start = "", end = "\n") {
        cat(start)
        cat(paste(rep(char, 80), collapse = ""))
        cat(end)
    }
    status        <- ""
    Labels        <- ""
    model         <- ""
    Anova         <- matrix(nrow = 0, ncol = 0)
    Potency       <- matrix(nrow = 0, ncol = 0)
    Estimates     <- matrix(nrow = 0, ncol = 0)
    inv           <- function (x) x
    factor        <- 1
    dilutionRatio <- 2
    dfAdjustment  <- 0
    Project       <- ""
    Assay         <- ""
    Comment       <- ""
    Description   <- ""
    Resume        <- ""
    Date          <- ""
    Operator      <- ""
    try(DataLabels <- read.table(file, comment.char = "",
                                 header = TRUE, fill = useFill))
    if (!is.null(DataLabels)) {
        Data     <- DataLabels
        Header   <- dimnames(Data)[[2]]
        if (readingTable)
            Plates   <- substr(unlist(Data[1]), 1,
                               nchar(firstName) + 1) == paste(firstName,
                                  "-", sep = "")
        else
            Plates   <- substr(unlist(Data[1]), 1, 1) != "#"
    }
    if (!all(Plates)) {
        Row1              <- as.character(Data[!Plates, ][[1]])
        Row2              <- as.character(Data[!Plates, ][[2]])
        Row3              <- as.character(Data[!Plates, ][[3]])
        RowN              <- as.character(Data[!Plates, ][[dim(Data)[2]]])
        Comments          <- any(Row1 == "#")
        if (Comments) {
            rowLabelNames <- 2
            LabelNames    <- Row2
        } else {
            rowLabelNames <- 1
            LabelNames    <- Row1
        }
        IsValues          <- substr(LabelNames, 3, 3) == "@"
        expectedDots      <- Data[!Plates,
                                  -c(1:rowLabelNames,
                                     dim(Data)[2])][!IsValues, ]
        invalidComments   <- apply(expectedDots != ".", 1, any)

        NonPlateTails     <- Data[!Plates, -c(1:(rowLabelNames+1))]
        if (is.null(dim(NonPlateTails)))
            LabelTails    <- NonPlateTails[!IsValues]
        else LabelTails   <- NonPlateTails[!IsValues, ]
        blankEnds         <- is.na(LabelTails) | LabelTails == ""
        if (all(blankEnds)) {
            cat("Note: You omitted dots! \n")
            NamesAndLabels  <- Data[!Plates,
                                    rowLabelNames + c(0, 1)][!IsValues, ]
            invalidComments <- FALSE
        }
        else NamesAndLabels <- Data[!Plates, c(rowLabelNames,
                                               dim(Data)[2])][!IsValues, ]
        NamesLabels   <- as.character(unlist(NamesAndLabels[, 1]))
        Labels        <- as.character(unlist(NamesAndLabels[, 2]))
        names(Labels) <- NamesLabels
        Assay         <- Labels[NamesLabels == "ResultOfAssay"]
        Project       <- Labels[NamesLabels == "Project"]
        Comment       <- Labels[NamesLabels == "Comment"]
        Description   <- Labels[NamesLabels == "Description"]
        Resume        <- Labels[NamesLabels == "Resume"]
        Operator      <- Labels[NamesLabels == "Operator"]
        Date          <- Labels[NamesLabels == "Date"]
        factor <- 1
        if (any(names(Labels) == "Factor"))
            factor <- unlist(lapply(as.list(strsplit(Labels["Factor"],
                                                     ",")[[1]]),
                                    FUN = function(x)
                                        eval(parse(text = x))))
        if (any(names(Labels) == "FactorA")) {
            factor <- unlist(lapply(as.list(strsplit(Labels["FactorA"],
                                                     ",")[[1]]),
                                    FUN = function(x)
                                        eval(parse(text = x))))
            if (any(names(Labels) == "FactorB"))
                factor <- unlist(lapply(as.list(strsplit(Labels["FactorB"],
                                                         ",")[[1]]),
                                        FUN = function(y, x = factor)
                                            eval(parse(text = y))))
        }
        dilutionRatio <- 1
        if (any(names(Labels) == "DilutionRatio"))
            dilutionRatio <- eval(parse(text = Labels["DilutionRatio"]))
        dfAdjustment <- 1
        if (any(names(Labels) == "DfAdjustment"))
            dfAdjustment <- eval(parse(text = Labels["DfAdjustment"]))
        model <- ""
        if (any(names(Labels) == "Model"))
            model <- Labels["Model"]
        design <- ""
        if (any(names(Labels) == "Design"))
            design <- Labels["Design"]
        if (any(invalidComments)) {
            x <- which(invalidComments)
            cat("\n")
            for (i in x) {
                y <- expectedDots[i, ] != "."
                cat("Warning: Problem with label number  ")
                cat(i)
                cat(" , ")
                cat(NamesLabels[i])
                cat(" : ")
                cat(expectedDots[i, ][y])
                cat("\n")
            }
            cat("\n")
        }
        if (Comments) {
            Results <- Data[!Plates, 2:dim(Data)[2]][IsValues, ]
        } else {
            Results <- Data[!Plates, ][IsValues, ]
        }
        Char4 <- substr(Results[, 1], 4, 4)
        if (any(Char4 == "<")) {
            Result1 <- Results[Char4 == "<", ]
            Result2 <- Results[Char4 == ">", -1]
            Results <- cbind(Result1, Result2)
        }
        else {
            Char2 <- substr(Results[, 1], 2, 2)
            if (any(Char2 == "2")) {
                Result1 <- Results[Char2 == "1", ]
                Result2 <- Results[Char2 == "2", -1]
                Results <- cbind(Result1, Result2)
            }
        }
        Char <- substr(Results[, 1], 1, 1)
        AnovaC <- Results[Char == "A", ]
        if (dim(AnovaC)[1] > 1) {
            if (dim(AnovaC)[2] > 6)
                AnovaC <- AnovaC[, 1:6]
            Names <- AnovaC[, 1]
            Anova <- suppressWarnings(matrix(
                as.numeric(as.character(unlist(AnovaC[, -1]))),
                ncol = 5))
            dimnames(Anova)[[1]] <- substr(as.character(Names), 4, 40)
            dimnames(Anova)[[2]] <- c("DF", "SumOfSquares", "MeanSquare",
                                      "Fratio", "P-value")
        }
        else Anova <- matrix(nrow = 0, ncol = 0)
        PotencyC <- Results[Char == "P", c(1, 3:5)]
        if (dim(PotencyC)[1] > 1) {
            Names <- PotencyC[, 1]
            Potency <- suppressWarnings(matrix(as.numeric(
                as.character(unlist(PotencyC[, -1]))), ncol = 3))
            dimnames(Potency)[[1]] <- substr(as.character(Names), 4, 40)
            dimnames(Potency)[[2]] <- c("Lower", "Estimate", "Upper")
        }
        else Potency <- matrix(nrow = 0, ncol = 0)
        EstimatesC <- Results[Char == "C", c(1, 3)]
        if (dim(EstimatesC)[1] > 1) {
            Names <- EstimatesC[, 1]
            Estimates <- suppressWarnings(matrix(
                as.numeric(as.character(unlist(EstimatesC[, -1]))),
                ncol = 1))
            dimnames(Estimates)[[1]] <- substr(as.character(Names), 4, 40)
            dimnames(Estimates)[[2]] <- c("Estimate")
        }
        else Estimates <- matrix(nrow = 0, ncol = 0)
        if (echoTitle) {
            myLine()
            cat("Project: ")
            cat(Project)
            cat("\nAssay: ")
            cat(Assay)
            myLine(, start = "\n")
        }
    }
    if (TRUE | !all(Plates)) {
        if (is.null(fun)) {
            fun <- function(x) x
            if (any(names(Labels) == "FUN"))
                fun <- eval(parse(text = Labels["FUN"]))
        }
        inv <- function(x) x
        if (any(names(Labels) == "INV"))
            inv <- eval(parse(text = Labels["INV"]))
    }
    if (log == TRUE)
        log <- exp(1)
    else
        if (log == FALSE)
            log <- 0
    if (is.null(fun))
        fun <- function(x) x
    ## assayData <- list(
    if (readingTable)
        class <- "assayTable"
    else
        class <- "assayFrame"
    assayData <- new(class,
        internals = list(Source = file,
            firstName = firstName,
            firstColumn = firstColumn,
            columnNames = columnNames,
            description = description,
            Rows = rows, Columns = columns,
            Color = color, Font = font),
        status = "", labels = Labels,
        anova = Anova,
        potency = Potency,
        estimates = Estimates,
        tableRaw = TableChar,
        log = log, fun = fun, inv = inv,
        factor = factor, dilutionRatio = dilutionRatio,
        dfAdjustment = dfAdjustment,
        projectTitle = Project, assayTitle = Assay,
        comment = Comment, description = Description,
        resume = Resume, operator = Operator,
        date = Date, model = model, design = design
        )
    return(assayData)
}
