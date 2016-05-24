# Read in a file with multiple fasta records. Return a vector of
# individual fasta record sequences, each names with the fasta
# identifier as name.
read.fasta <- function(file.name) {
    x <- scan(file.name, what=character(), sep='\n', quiet=T)

    id.pos <- grep('>', x, fixed=T)
    stopifnot(id.pos[1] == 1)

    ids <- x[id.pos]
    # Trim the > and \n from each id
    ids <- gsub('>|\n', '', ids)

    zz <- character(length(ids))
    names(zz) <- ids

    starts <- id.pos+1
    ends <- id.pos-1
    ends <- ends[-1]
    ends[length(ids)] <- length(x)
    for (i in 1:length(ids)) {
        zz[i] <- paste(x[starts[i]:ends[i]], collapse='')
        # Remove white spaces useing perl regular expression syntax
        # (\W == whitespace)
        zz[i] <- gsub('\\W', '', zz[i], perl=T)
        zz[i] <- toupper(zz[i])
    }
    zz
}

# Find the mononucleotide repeats in text and return a data frame
# with <chr>, <start>, <stop> in each row.
#
# The input, text, is a vector of character strings; each character
# string should consist only of the nucleotides in a DNA sequence.
# The name of each string in text is taken as the chromosome name
# for the output data frame.
find.mono.repeats <- function(text, min.len=5) {
    pattern <- paste(paste(c('A', 'C', 'G', 'T'), '{', min.len, ',}', sep=''), collapse='|')
    mono.repeats <- gregexpr(pattern=pattern, text=text, perl=F)
    stopifnot(length(mono.repeats) == length(names(text)))
    output.table <- data.frame(chr=character(0), start=numeric(0), stop=numeric(0))
    for (i in 1:length(mono.repeats)) {
        if (mono.repeats[[i]][1] == -1) next # -1 if there was no match at all
		stop <- mono.repeats[[i]] + attr(mono.repeats[[i]], 'match.length') - 1
        tt <-data.frame(chr=names(text)[i], start=mono.repeats[[i]], stop=stop)
        output.table <- rbind(output.table, tt)
    }
    output.table
}

