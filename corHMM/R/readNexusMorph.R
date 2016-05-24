
######################################################################################################################################
######################################################################################################################################
### Example run
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu and John Nylander

#This is a hacked version of apes "read.nexus.data". I rewrote parts to deal with missing, complete and incomplete ambiguous characters, and so that format of the data fits perfectly into corHMM.

readNexusMorph <- function(file) {
    find.ntax <- function(x){
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        ntax
    }
    find.nchar <- function(x){
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        nchar
    }
    find.matrix.line <- function(x){
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
                matrix.line <- as.numeric(i)
                break
            }
        }
        matrix.line
    }
    trim.whitespace <- function(x){
        gsub("\\s+", "", x)
    }
    trim.semicolon <- function(x){
        gsub(";", "", x)
    }
    X <- scan(file = file, what = character(), sep = "\n",
    quiet = TRUE, comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0
    seq.mat <- matrix(0, ntax, nchar+1)
    count <- 1
    for (j in start.reading:NROW(X)) {
        Xj <- trim.semicolon(X[j])
        if(Xj == "") {
            break
        }
        if(any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE))) {
            break
        }
        ts <- unlist(strsplit(Xj, "(\\s+)", perl = TRUE))
        if (length(ts) > 2) {
            stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
        }
        if (length(ts) !=2) {
            stop("nexus parser failed to read the sequences (ts!=2)")
        }
        Name <- trim.whitespace(ts[1])
        Seq <- trim.whitespace(ts[2])
        Seq <- strsplit(Seq, NULL)[[1]]
        #Find gap characters and replace with "?":
        Seq[which(Seq=="-")] = "?"
        #Find ambiguous characters:
        open.bracket.index <- which(Seq=="{")
        close.bracket.index <- which(Seq=="}")
        new.Seq <-c()
        for(index in sequence(length(Seq)-sum(close.bracket.index - open.bracket.index))){
            open.bracket.index <- which(Seq=="{")
            close.bracket.index <- which(Seq=="}")
            if(any(index == open.bracket.index)){
                open.bracket <- open.bracket.index[which(index == open.bracket.index)]
                closed.bracket <- close.bracket.index[which(index == open.bracket.index)]
                new.char <- Seq[open.bracket+1]
                for(k in (open.bracket+2):(closed.bracket-1)){
                    new.char<-paste(new.char, Seq[k], sep="&")
                }
                Seq <- Seq[!sequence(length(Seq)) %in% c((open.bracket+1):closed.bracket)]
                Seq[index] <- new.char
            }
        }
        seq.mat[count,1] <- Name
        seq.mat[count,2:dim(seq.mat)[2]] = Seq
        count <- count+1
    }
    seq.mat <- as.data.frame(seq.mat)
    seq.final <- seq.mat[-1]
    rownames(seq.final) <- seq.mat[,1]
    colnames(seq.final) <- sequence(dim(seq.final)[2])
    return(seq.final)
}

