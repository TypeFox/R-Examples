readFastaFile=function (fileName, sep=" "){
    fastaFile=file(fileName, "r")  
    sequences=list()
    line.text=readLines(fastaFile,1)
    name=NULL
    while (length(line.text)>0) {
        if (substr(line.text, 1,1)==">") {        
            name=line.text
            name=substr(strsplit (name, sep)[[1]][1], 2, 1000)            
            temp.seq=""    
            line.text=readLines(fastaFile,1)
            while ( length(line.text)>0 ) {
                if (substr(line.text, 0, 1) == ">") break;
                temp.seq=temp.seq %+% line.text
                line.text=readLines(fastaFile,1)
            }    
            sequences[[name]]=temp.seq
        } else {
            print ("sth wrong")
        }
    }    
    close(fastaFile)
    sequences
}

writeFastaFile=function (seqList, fileName) {
    outFile=file (fileName, open="w")
    for (i in 1:length(seqList)){
        write (file=outFile, ">"%+%names(seqList)[i], append=T)
        write (file=outFile, seqList[[i]], append=T)
    }
    close(outFile)
}

aaList=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")

aa2arabic=function (seq1) {
    temp=strsplit(seq1, "|")[[1]] # separate every character
    out=sapply(temp, function (i) {
        tmp=which (aaList==toupper(i)); 
        ifelse(length(tmp)==0, 21, tmp)
    } )
    names(out)=NULL
    out
}

string2arabic=function (seqList) {
# if out is not what you expect, maybe the strings don't have the same length
    out=mysapply(seqList, aa2arabic)
    dimnames(out)[[2]]=1:ncol(out)
    out
}

fastaFile2arabicFile=function(fastaFile, arabicFile, removeGapMajor=FALSE){
    stringList=readFastaFile(fastaFile)
    stringList2arabicFile(stringList, arabicFile, removeGapMajor)
}

selexFile2arabicFile=function(selexFile, arabicFile, removeGapMajor=FALSE){
    stringList=readSelexFile(selexFile)
    stringList2arabicFile(stringList, arabicFile, removeGapMajor)
}

stringList2arabicFile=function (seqList, arabicFile, removeGapMajor=FALSE) {
    alignment=string2arabic(seqList)
    # remove columns that are predominantly gaps
    if (removeGapMajor) {
        tmpCount=alignment2count(alignment, level=21)
        alignment=alignment[,tmpCount[,21]<10]
    }
    arabic2arabicFile (alignment, arabicFile)
}

arabic2arabicFile=function (alignment, arabicFile) {
    # adjust to 0 based index before writing to file 
    alignment=alignment-1 
    n=nrow(alignment)
    p=ncol(alignment)
    m=max(alignment)
    write (c(n, p, m), file=arabicFile)
    write (t(alignment), file=arabicFile, ncolumns=p, append=T)
    invisible(alignment+1)
}

# return a list of strings
readSelexFile=function (fileName) {
    fileCon=file(fileName, open="r")
    line1=readLines(fileCon, 1)
    while (startsWith(line1, "#") | line1=="") line1=readLines(fileCon, 1)
    seqs=list()
    while(length(line1)>0){
        tmp=strsplit(line1, " ")[[1]]
        seqs[[tmp[1]]]=tmp[length(tmp)]
        line1=readLines(fileCon, 1)
    }    
    close (fileCon)
    seqs
}

readSelexAsMatrix=function (fileName) {
    fileCon=file(fileName, open="r")
    outM=NULL
    names=NULL
    while(T){
        seqName=scan (fileCon, what="character", n=1)
        seq=scan (fileCon, what="character", n=1)
        if (length(seqName)>0) {
            outM=rbind(outM, strsplit(seq,"")[[1]])
            names=c(names, seqName)
        }else{
            break;
        }
    }
    close (fileCon)
    dimnames(outM)=list(names, NULL)
    outM
}

arabic2fastaFile=function (alignment, fileName){
    outFile=file (fileName, open="w")
    name=rownames(alignment)
    if (is.null(name)) name=1:nrow(alignment)
    for (i in 1:nrow(alignment)){
        write (file=outFile, ">"%+%name[i], append=T)
        write (file=outFile, concatList(aaList[alignment[i,]]), append=T)
    }
    close(outFile)
}

readArabicFile=function (fileName) {
    infile=file (fileName, open="r")
    n=scan (infile, n=1, quiet=T)
    p=scan (infile, n=1, quiet=T)
    m=scan (infile, n=1, quiet=T)
    out=mysapply (1:n, function (i) as.numeric( strsplit(readLines(infile, n=1), " ")[[1]] ) )
    close(infile)
    invisible (out+1)
}

readBlockFile=function (fileName) {
    file1= file (fileName, "r")
    n=scan (file1, n=1, quiet=T)
    p=scan (file1, n=1, quiet=T)
    M=scan (file1, n=1, quiet=T)
    out=matrix(scan (file1, quiet=T), nrow=n, byrow=T)+1
    close(file1)
    out
}
