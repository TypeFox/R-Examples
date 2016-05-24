
EViewsHeader <- function() {
    memFormat(firstline=memBlock(80),
              headersize=integer8,
              unknown=memBlock(26),
              numvblesplusone=integer4,
              date=vectorBlock(ASCIIchar, 4),
              unkown=memBlock(2),
              datafreq=integer2,
              startperiod=integer2,
              startobs=integer4,
              unkown=memBlock(8),
              numobs=integer4)
}

EViewsVbleInfo <- function() {
    memFormat(unknown=memBlock(6),
              recsize=integer4,
              memsize=integer4,
              ptrtodata=integer8,
              vblename=vectorBlock(ASCIIchar, 32),
              ptrtohistory=integer8,
              vbletype=integer2,
              unknown=memBlock(6))
}

EViewsVbleData <- function(numObs) {
    memFormat(numobs=integer4,
              startobs=integer4,
              unknown=memBlock(8),
              endobs=integer4,
              unknown=memBlock(2),
              values=vectorBlock(real8, 
                numObs))
}

readEViews <- function(filename, as.data.frame=TRUE) {
    # Start with the header information
    header <- readFormat(filename, EViewsHeader())
    vbleStart <- blockValue(header$blocks$headersize) + 24 + 2
    # Going to ignore "C" and "RESID"
    numVbles <- blockValue(header$blocks$numvblesplusone) - 1 - 2
    numObs <- blockValue(header$blocks$numobs)
    # Empty container for data set
    Names <- vector("character", numVbles)
    Data <- vector("list", numVbles)
    vbleCount <- 1
    for (i in 1:(numVbles + 2)) {
        vbleInfo <- readFormat(filename,
                               offset=vbleStart + (i - 1)*70,
                               EViewsVbleInfo())
        vbleName <- blockString(vbleInfo$blocks$vblename)
        if (vbleName == "C" || vbleName == "RESID") {
            warning("Skipping boilerplate variable\n")
        } else {
            Names[vbleCount] <- vbleName
            dataLoc <- blockValue(vbleInfo$blocks$ptrtodata)
            vbleValues <- readFormat(filename,
                                     offset=dataLoc,
                                     EViewsVbleData(numObs))
            Data[[vbleCount]] <- blockValue(vbleValues$blocks$values)
            vbleCount <- vbleCount + 1
        }
    }
    names(Data) <- Names
    if (as.data.frame)
        as.data.frame(Data)
    else
        Data
}
