`report.mefa` <-
function (x, filename, segment = FALSE, n = NULL, by.taxa = TRUE,
    samp.var = NULL, drop.redundant = NULL, collapse = TRUE,
    taxa.name = NULL, author.name = NULL, taxa.order = NULL,
    grouping = FALSE, tex = FALSE, binary = FALSE,
    tex.control = list(ital.taxa = TRUE, noindent = TRUE,
    bold.sect = TRUE, bold.1st = TRUE, vspace1 = 0.5, vspace2 = 0.2),
    sep = c(",", ":", "(", ":", ",", ")", ";"), dir = getwd(), ...)
{
current.dir <- getwd()
setwd(dir)
on.exit(setwd(current.dir))

mf <- x[drop = TRUE]
if (is.null(taxa.name)) {
    taxa.name <- 1
    mf$taxa <- data.frame(rownames(mf$taxa), mf$taxa)}

# test
    if(!is.mefa(mf))
        stop("object is not of class 'mefa'")
    if (is.null(mf$samp) || is.null(mf$taxa))
        stop("report needs both '$samp' and '$taxa'")
    if (length(sep) != 7)
        stop("specify exactly 7 'sep' values")
    if (is.null(taxa.order))
        taxa.order <- taxa.name
    if (is.character(n))
        n <- which(names(mf$segm) %in% n)
    length.n <- 1
    if (!is.null(n) && is.null(mf$segm))
        stop("no segments found")
    if (segment && is.null(mf$segm))
        stop("no segments found")
    if (segment && is.null(n))
        length.n <- dim(mf)[3]
    if (segment && !is.null(n))
        length.n <- length(n)

# ordering sample attributes
    if (is.null(samp.var))
        loca <- mf$samp else loca <- mf$samp[, samp.var]
    nnloca <- loca
    names(nnloca) <- NULL
    ord_loca <- do.call(order, nnloca)
    loc <- loca[ord_loca, ]

    if (!is.null(drop.redundant) & length(drop.redundant) >= ncol(loca))
        stop("'drop.redundant' should be smaller than length of 'samp.var'")

# total count data
    if (!segment) {
        mfdata <- mf$xtab
        mfd <- list(mf$xtab)
        }

    if (segment) {
        if (is.null(n)) {
            mfdata <- mf$xtab
            mfd <- list()
            for (i in 1:dim(mf)[3]) {
                runname <- paste(dimnames(mf)$segm[i])
                mfd[[runname]] <- mf$segm[[i]]
                }
            }
        if (!is.null(n)) {
            if (length(n) == 1) mfdata <- mf$segm[[n]]
            if (length(n) > 1) {
                mfdata <- mf$segm[[n[1]]]
            for (i in 2:length(n)) mfdata <- mfdata + mf$segm[[n[i]]]
            }
            mfd <- list()
            for (i in 1:length(n)) {
                runname <- paste(dimnames(mf)$segm[n[i]])
                mfd[[runname]] <- mf$segm[[n[i]]]
                }
        }
    }

# ordering counts
    xcr <- mfdata[ord_loca, order(mf$taxa[, taxa.order])]
    mfdl <- list()
    for (i in 1:length.n) {
        mfdl[[i]] <- mfd[[i]][ord_loca, order(mf$taxa[, taxa.order])]
        }
    names(mfdl) <- names(mfd)

# species names
    nam <- as.vector(mf$taxa)[, taxa.name][order(mf$taxa[, taxa.order])]
    if (!is.null(author.name))
        autv <- as.vector(mf$taxa)[, author.name][order(mf$taxa[, taxa.order])]

# formatting
    if (tex & tex.control$ital.taxa) ti <- "\\textit{" else ti <- ""
    if (tex & tex.control$bold.1st) tb <- "\\textbf{" else tb <- ""
    if (tex & tex.control$bold.sect) tb2 <- "\\textbf{" else tb2 <- ""
    if (tex) noin <- "\\noindent " else noin <- ""
    if (tex & tex.control$noindent) noin1 <- "\\noindent " else noin1 <- ""
    if (tex & grouping & tex.control$noindent) noin2 <- "\\noindent " else noin2 <- ""
    if (tex & tex.control$ital.taxa) clbr <- "}" else clbr <- ""
    if (tex & tex.control$bold.1st) clbr2 <- "}" else clbr2 <- ""
    if (tex & tex.control$bold.sect) clbr3 <- "}" else clbr3 <- ""
    vspace1 <- tex.control$vspace1
    vspace2 <- tex.control$vspace2
    calling <- deparse(match.call())
    calling <- gsub("  ", "", calling)

# START ordering=species
    if (by.taxa) {

        zz <- file(filename, "w")
        cat("%% Start writing data from a 'mefa' object sorted by species into file \"",
            filename, "\" on ", date(), ".\n%% Call: ", calling, "\n", file = zz, sep = "")

# start of SPEC loop
    for (spec in 1:length(nam)) {

# start of IF >0
    if (sum(xcr[, spec]) > 0) {

# specname cat
            if (is.null(author.name)) {
                spprint <- paste(noin,tb2,ti,nam[spec],clbr,clbr3,"\n\n",sep="")
                } else{
                spprint <- paste(noin,tb2,ti,nam[spec],clbr," ",autv[spec],clbr3,"\n\n",sep="")
                }
            cat("\n\n", file = zz, sep = "")
            if (tex) cat(paste("\\vspace{",vspace1,"cm} ",sep=""), file = zz, sep = "")
            cat(spprint, file = zz, sep = "")
            if (tex & !grouping) cat(paste("\\vspace{",vspace2,"cm} ",sep=""), file = zz, sep = "")

# nonzero count subsets
            loc.sub <- as.matrix(subset(loc, xcr[, spec] > 0))
            lev.sub <- as.factor(loc.sub[,1])
            xcr.sub <- subset(xcr[, spec], xcr[, spec] > 0)
            mfdl.sub <- matrix(NA,length(xcr.sub),length.n)
            colnames(mfdl.sub) <- names(mfdl)
            for (i in 1:length.n) {
               mfdl.sub[,i] <- as.vector(subset(mfdl[[i]][, spec], xcr[, spec] > 0))
               }

# collapse - exclude and aggregate
    leave <- rep(1,nrow(loc.sub))
    aggr <- c(1:nrow(loc.sub))
    if (collapse) {
        for (i in 1:nrow(loc.sub)){
            if (i > 1) if (sum(loc.sub[i,] == loc.sub[(i-1),]) == ncol(loc.sub)) leave[i] <- 0
            aggr[i] <- paste(loc.sub[i,],collapse="")
            }
        aggr <- as.numeric(as.factor(aggr))
        loc.sub <- subset(loc.sub, leave==1)
        lev.sub <- subset(lev.sub, leave==1)
        xcr.sub <- aggregate(xcr.sub,list(aggr),sum)[,2]
        mfdl.sub <- aggregate(mfdl.sub,list(aggr),sum)[,-1]
        }

    if (!is.null(drop.redundant)) {
        kloc.sub <- as.matrix(loc.sub)
    if (nrow(kloc.sub) > 1) {
        for (col in 1:drop.redundant) {
            for (row in 1:(nrow(loc.sub)-1)) {
                if (loc.sub[row, col] == loc.sub[(row+1), col]) kloc.sub[(row+1), col] <- ""
            }
        }
    if (drop.redundant > 1) {
        for (col in 2:drop.redundant) {
            for (row in 1:nrow(loc.sub)) {
                if (kloc.sub[row, (col-1)] != "") kloc.sub[row, col] <- loc.sub[row, col]
            }
        }}
        for (col in 1:drop.redundant) {
            for (row in 1:(nrow(loc.sub)-1)) {
                if (kloc.sub[row, col] != "" & kloc.sub[(row+1), col] == "")
                    kloc.sub[row, col] <- paste(sep[2],sep[2],kloc.sub[row, col],sep[2]," ",sep="")
            }
        }
    } # end IF nrow > 1
    } else kloc.sub <- loc.sub

    for (j in 1:ncol(kloc.sub)) {
        for (i in 1:nrow(kloc.sub)) {
            test.1st <- substr(kloc.sub[i,j],1,2) == paste(sep[2],sep[2],sep="")
            if (kloc.sub[i,j] != "" & !test.1st & j!=ncol(kloc.sub))
                kloc.sub[i,j] <- paste(kloc.sub[i,j],sep[1]," ",sep="")
            if (kloc.sub[i,j] != "" & !test.1st & j==ncol(kloc.sub))
                kloc.sub[i,j] <- paste(kloc.sub[i,j]," ",sep="")
            if (test.1st) kloc.sub[i,j] <- gsub(paste(sep[2],sep[2],sep=""),"",kloc.sub[i,j])
        }
    }

if (!grouping) cat(noin1, file = zz, sep = "")

# loop for first column levels (grouping also)
        for (lev in 1:nlevels(lev.sub)) {
            if (tex & grouping) cat(paste("\\vspace{",vspace2,"cm} ",sep=""), file = zz, sep = "")
            xcr.sub2 <- subset(xcr.sub, lev.sub == levels(lev.sub)[lev])
            loc.sub2 <- subset(loc.sub, lev.sub == levels(lev.sub)[lev])
            kloc.sub2 <- subset(kloc.sub, lev.sub == levels(lev.sub)[lev])
            printcount <- rep("",length(lev.sub[lev.sub==levels(lev.sub)[lev]]))

if (segment){
        mfdl.sub2 <- subset(mfdl.sub, lev.sub == levels(lev.sub)[lev])
        mfdl.sub2[mfdl.sub2==0]  <- paste("DELETEME",sep[5],sep="")
        for (i in 1:nrow(mfdl.sub2)) {
            for (j in 1:ncol(mfdl.sub2)) {
                if (mfdl.sub2[i,j] != paste("DELETEME",sep[5],sep=""))
                    if (binary) {
                        mfdl.sub2[i,j] <- paste(colnames(mfdl.sub2)[j],sep="")
                            } else {
                        mfdl.sub2[i,j] <- paste(colnames(mfdl.sub2)[j],sep[4]," ",mfdl.sub2[i,j],sep="")}
            }
            printcount[i] <- paste(mfdl.sub2[i,],collapse=paste(sep[5]," ",sep=""))
            printcount[i] <- gsub(paste("DELETEME",sep[5],sep[5]," ",sep=""),"", printcount[i])
            printcount[i] <- gsub(paste("DELETEME",sep[5],sep=""),"", printcount[i])
        }
    }
if (!segment) {
    for (i in 1:length(printcount))
        if (binary) {
            printcount[i] <- paste("",sep="")
                } else {
            printcount[i] <- paste(xcr.sub2[i],sep="")}
    }

if (binary & !segment) {
    brace1 <- rep("",length(lev.sub[lev.sub==levels(lev.sub)[lev]]))
    brace2 <- rep("",length(lev.sub[lev.sub==levels(lev.sub)[lev]]))
        } else {
    brace1 <- rep(sep[3],length(lev.sub[lev.sub==levels(lev.sub)[lev]]))
    brace2 <- rep(sep[6],length(lev.sub[lev.sub==levels(lev.sub)[lev]]))
    }

ending <- c(rep(paste(sep[7]," ",sep=""),length(lev.sub[lev.sub==levels(lev.sub)[lev]])-1),". ")

printout0 <- data.frame(kloc.sub2, brace1, printcount, brace2, ending)
colnames(printout0) <- letters[1:ncol(printout0)]
printout <- as.matrix(printout0)
paragraph <- ""

for (i in 1:nrow(printout)) {
    for (j in 1:ncol(printout)) {
    if (i==1 & j==1) printout[i,j] <- paste(noin2,tb,as.character(printout[i,j]),clbr2,sep="")
        else printout[i,j] <- paste(as.character(printout[i,j]),sep="")
    paragraph <- paste(paragraph,printout[i,j],sep="")
    }}

    paragraph <- gsub("  "," ", paragraph)
    paragraph <- gsub(paste(sep[5]," ",sep[6],sep=""),paste(sep[6],sep=""), paragraph)
    paragraph <- gsub(paste(" ",sep[7]," ",sep=""),paste(sep[7]," ",sep=""), paragraph)
    paragraph <- gsub(" . ",". ", paragraph)
    paragraph <- gsub(paste(sep[5],". ",sep=""),". ", paragraph)

    cat(paragraph, file = zz, sep = "")

if (grouping) cat("\n\n", file = zz, sep = "")
} # end of LEV loop
} # end of IF >0
} # end of SPEC loop

        cat("\n\n%% End of output.\n", file = zz, sep = "")
        close(zz)
    } #END of species ordering

if (!by.taxa) stop("'by.taxa = FALSE' is not yet implemented\n")

invisible()
} # end of function


