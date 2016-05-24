read.fasta <- function
### Read sequences in FASTA format into a named character vector
(infile
### Name of the sequence file
 ){
  tmp <- sapply(strsplit(strsplit(paste('\n\n\n',paste(readLines(infile),collapse='\n'),sep=''),split='\n>')[[1]],'\n'),function(v)c(v[1],paste(v[2:length(v)],collapse='')))
  seqs <- tmp[2,]
  names(seqs) <- tmp[1,]
  blank <- names(seqs)==''
  names(seqs)[blank] <- seqs[blank]
  names(seqs) <- gsub(' ','',names(seqs))
  seqs[seqs!='']
}
##debug(read.fasta)

### Letters in the DNA alphabet, used to auto-detect sequence type
dna.letters <- c("*","A","T","G","C")
### DNA identity substitution matrix.
dna.identity <- matrix(0,nrow=length(dna.letters),ncol=length(dna.letters),
                       dimnames=list(dna.letters,dna.letters))
diag(dna.identity) <- 1
dna.identity['*','*'] <- 0
seqs.to.mat <- function
### Calculate pairwise differences between sequences using a
### substitution matrix.
(seq.vec,
### DNA or protein sequences.
 subs.mat=NULL){
### Substitution matrix with dimnames that match the letters used in
### the sequence data, or a character vector that specifies a common
### substitution matrix (as defined in the Biostrings package). NULL
### specifies that we will guess a suitable substitution matrix to
### match your input sequences (DNA=>identity, protein=>BLOSUM62).
  if(is.null(names(seq.vec)))names(seq.vec) <- seq.vec
  chars <- sapply(seq.vec,nchar)
  seqsum <- table(chars)
  if(length(seqsum)>1){
    print(as.data.frame(seqsum))
    i <- which.max(seqsum)
    cat("Sequences not of length ",names(seqsum)[i],":\n",sep="")
    rows <- as.integer(names(seqsum)[i])!=chars
    print(data.frame(chars,name=names(seq.vec),row.names=NULL)[rows,])
    stop("All input sequences must be of the same length.")
  }
  d <- toupper(gsub('[- .]',"*",seq.vec[unique(names(seq.vec))]))
  print(d)
  letters <- unique(c(unlist(strsplit(d,split='')),dna.letters))
  ##if dna alignment use simple identity matrix
  looks.like.dna <- identical(sort(letters),sort(dna.letters))
  if(is.null(subs.mat))
    subs.mat <- if(looks.like.dna)dna.identity else "BLOSUM62"
  if(mode(subs.mat)=="character"){
    ex <- substitute(data(M,package="Biostrings"),list(M=subs.mat))
    eval(ex)
    subs.mat <- get(subs.mat)
  }
  print(subs.mat)
  N <- length(d)
  m <- matrix(0,nrow=N,ncol=N,dimnames=list(names(d),names(d)))
  for(i in 1:N)for(j in 1:i){
    seqs <- sapply(strsplit(c(d[i],d[j]),split=''),c)
    entry <- try(apply(seqs,1,function(x)subs.mat[x[1],x[2]]))
    if(class(entry)=="try-error"){
      print(seqs)
      stop("Sequence difference matrix construction failed.")
    }
    ## subscript out of bounds here usually means bad matrix
    m[i,j] <- m[j,i] <- -sum(entry)
  }
  m <- m-min(m)
  attr(m,'seqs') <- d
  print(m[1:min(5,nrow(m)),1:min(5,ncol(m))])
  m
### The matrix of distances between each input sequence, with dimnames
### corresponding to either the sequences, or the sequence names (if
### they exist)
}
##debug(seqs.to.mat)

make.logo.ps <- function
### Create a logo using weblogo, then read it in using grImport
(helices,
### Sequences to plot in the logo
 psbase
### Base filename for the logo postscript and xml files, should be the
### full path name except for the trailing .ps
 ){
  psfile <- paste(psbase,'ps',sep='.')
  xmlfile <- paste(psfile,'xml',sep='.')
  seq.text <- paste(paste('>',helices,'\n',helices,sep=''),collapse='\n')
  write(seq.text,psbase)
  execdir <- system.file("exec",package="sublogo")
  cmd <- paste("PATH=",execdir,
               ":$PATH seqlogo -c -F EPS -f ",psbase,
               "|sed 's/^EndLine/%EndLine/'|sed 's/^EndLogo/%EndLogo/' >",
               psfile,sep="")
  cat(cmd,'\n')
  system(cmd)
  owd <- setwd(tempdir())
  PostScriptTrace(psfile,xmlfile)
  setwd(owd)
  pic <- readPicture(xmlfile)
  pic
### Grid picture grob as read using readPicture
}
##debug(make.logo.ps)

sublogo.dendrogram <- function
### Main function for drawing sublogo dendrograms.
(M,
### difference matrix as constructed by seqs.to.mat (although in
### principle any object with a valid as.dist method could be used)
 main='',
### plot title
 subtit=NULL,
### plot subtitle
 base=tempfile(),
### base file name for temporary logo files
 cutline=150,
### Distance for cutting the tree. Draw a sublogo for each
### leaf. Normally you will plot once, then inspect the dendrogram to
### determine which is a good value for cutline, then plot again using
### your chosen cutline.
 dend.width=30,
### Percent of the plot to be occupied by the dendrogram. The logos
### will be displayed with equal widths on either side.
 cex=1
### character expansion factor for the dendrogram
 ){
  hc <- hclust(as.dist(M),method="average")
  dend <- as.dendrogram(hc)
  fam <- cutree(hc,h=cutline)[labels(dend)] # order by plotting method
  famids <- unique(fam)
  famtab <- data.frame(fam,y=1:length(fam))
  famtab$seq <- attr(M,'seqs')[rownames(famtab)]
  fam.nontriv <- sapply(famids,function(i)sum(famtab$fam==i))>1
  names(fam.nontriv) <- famids
  if(is.null(subtit))
    subtit <- paste(nrow(M),"sequences,",sum(fam.nontriv),"families")
  xrange <- c(0,1)
  yrange <- c(0,max(famtab[,'y'])+1)
  draw.box <- function(i){ # replace with logos
    subtab <- famtab[famtab[,'fam']==i,]
    grprange <- range(subtab[,'y'])
    rect(xrange[1],grprange[1]-0.25,xrange[2],grprange[2]+0.25)
  }
  draw.logo <- function(i){
    ## do not draw sublogo for a family of trivial size
    if(fam.nontriv[as.character(i)]){ 
      subtab <- famtab[famtab[,'fam']==i,]
      grprange <- range(subtab[,'y'])
      tmpfile <- paste(base,i,sep='.')
      logo <- make.logo.ps(subtab$seq,tmpfile)
      grid.picture(logo,
                   x=xrange[1],
                   y=unit(grprange[1]+0.25,'native'),
                   height=unit(diff(grprange)-0.5,'native'),
                   distort=TRUE,
                   just=c(0,0),
                   fillText=TRUE)
    }
  }

  bottomspace <- 0.6
  topspace <- 0.4
  side.percents <- (100-dend.width)/2
  layout(matrix(1:3,ncol=3),c(side.percents,dend.width,side.percents))
  par(mai=c(bottomspace,0,topspace,0),cex=1)
  
  ## Big summary logo on left
  plot(xrange,yrange,
       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')
  biglogo <- make.logo.ps(famtab$seq,paste(base,'0',sep='.'))
  vps <- baseViewports()
  pushViewport(vps$inner,vps$figure,vps$plot)
  grid.picture(biglogo,
               y=unit(0,'native'),
               height=unit(length(fam),'native'),
               just=c(0.5,0),
               exp=0,
               distort=TRUE,
               fillText=TRUE)
  popViewport(3)
  
  ## Dendrogram in middle
  par(mai=c(bottomspace,0,topspace,
        max(strwidth(colnames(M),'inches',
                     cex))/6*5),
      family='mono')
  par(xpd=NA)
  plot(dend,h=TRUE,edgePar=list(lwd=2),nodePar=list(lab.cex=cex,pch=""))
  par(family="")
  segments(cutline,1,cutline,length(fam))
  ##axis(3,cutline,lty=0,line=0)
  
  ## Title in the middle
  title(main,line=0.5)
  mtext(subtit,line=-0.7,cex=1)

  # Sublogos on the right side
  par(mai=c(bottomspace,0,topspace,0))
  plot(xrange,yrange,
       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')
  vps <- baseViewports()
  pushViewport(vps$inner,vps$figure,vps$plot)
  sapply(famids,draw.logo)
  popViewport(3)
  
  hc
### The dendrogram from the call to hclust
}
##debug(sublogo.dendrogram)

sublogo <- function
### Draw a sublogo dendrogram for a sequence motif.
(seqs,
### Character vector of DNA or protein sequences (optionally named
### with sequence labels).
 mat=NULL,
### Substitution matrix passed to seqs.to.mat.
 ...
### Other arguments to pass to sublogo.dendrogram (see that help page
### for a full description).
 ){
  sublogo.dendrogram(seqs.to.mat(seqs,mat),...)
}

.result <-
list(dna.identity = list(definition = "dna.identity <- matrix(0,nrow=length(dna.letters),ncol=length(dna.letters),\n                       dimnames=list(dna.letters,dna.letters))", format="",title="dna identity",
    description = "DNA identity substitution matrix."),
     dna.letters = list(format="",title="dna letters",
    definition = "dna.letters <- c(\"*\",\"A\",\"T\",\"G\",\"C\")", 
    description = "Letters in the DNA alphabet, used to auto-detect sequence type"), 
    make.logo.ps = list(definition = "make.logo.ps <- function\n### Create a logo using weblogo, then read it in using grImport\n(helices,\n### Sequences to plot in the logo\n psbase\n### Base filename for the logo postscript and xml files, should be the\n### full path name except for the trailing .ps\n ){\n  psfile <- paste(psbase,'ps',sep='.')\n  xmlfile <- paste(psfile,'xml',sep='.')\n  seq.text <- paste(paste('>',helices,'\\n',helices,sep=''),collapse='\\n')\n  write(seq.text,psbase)\n  execdir <- system.file(\"exec\",package=\"sublogo\")\n  cmd <- paste(\"PATH=\",execdir,\n               \":$PATH seqlogo -c -F EPS -f \",psbase,\n               \"|sed 's/^EndLine/%EndLine/'|sed 's/^EndLogo/%EndLogo/' >\",\n               psfile,sep=\"\")\n  cat(cmd,'\\n')\n  system(cmd)\n  owd <- setwd(tempdir())\n  PostScriptTrace(psfile,xmlfile)\n  setwd(owd)\n  pic <- readPicture(xmlfile)\n  pic\n### Grid picture grob as read using readPicture\n}",
      format="",title="make logo ps",
        description = "Create a logo using weblogo, then read it in using grImport", 
        `item{helices}` = "Sequences to plot in the logo", `item{psbase}` = "Base filename for the logo postscript and xml files, should be the\nfull path name except for the trailing .ps", 
        value = "Grid picture grob as read using readPicture"), 
    read.fasta = list(definition = "read.fasta <- function\n### Read sequences in FASTA format into a named character vector\n(infile\n### Name of the sequence file\n ){\n  tmp <- sapply(strsplit(strsplit(paste('\\n\\n\\n',paste(readLines(infile),collapse='\\n'),sep=''),split='\\n>')[[1]],'\\n'),function(v)c(v[1],paste(v[2:length(v)],collapse='')))\n  seqs <- tmp[2,]\n  names(seqs) <- tmp[1,]\n  blank <- names(seqs)==''\n  names(seqs)[blank] <- seqs[blank]\n  names(seqs) <- gsub(' ','',names(seqs))\n  seqs[seqs!='']\n}",
      format="",title="read fasta",
        description = "Read sequences in FASTA format into a named character vector", 
        `item{infile}` = "Name of the sequence file"),
     seqs.to.mat = list(format="",title="seqs to mat",
        definition = "seqs.to.mat <- function\n### Calculate pairwise differences between sequences using a\n### substitution matrix.\n(seq.vec,\n### DNA or protein sequences.\n subs.mat=NULL){\n### Substitution matrix with dimnames that match the letters used in\n### the sequence data, or a character vector that specifies a common\n### substitution matrix (as defined in the Biostrings package). NULL\n### specifies that we will guess a suitable substitution matrix to\n### match your input sequences (DNA=>identity, protein=>BLOSUM62).\n  if(is.null(names(seq.vec)))names(seq.vec) <- seq.vec\n  chars <- sapply(seq.vec,nchar)\n  seqsum <- table(chars)\n  if(length(seqsum)>1){\n    print(as.data.frame(seqsum))\n    i <- which.max(seqsum)\n    cat(\"Sequences not of length \",names(seqsum)[i],\":\\n\",sep=\"\")\n    rows <- as.integer(names(seqsum)[i])!=chars\n    print(data.frame(chars,name=names(seq.vec),row.names=NULL)[rows,])\n    stop(\"All input sequences must be of the same length.\")\n  }\n  d <- toupper(gsub('[- .]',\"*\",seq.vec[unique(names(seq.vec))]))\n  print(d)\n  letters <- unique(c(unlist(strsplit(d,split='')),dna.letters))\n  ##if dna alignment use simple identity matrix\n  looks.like.dna <- identical(sort(letters),sort(dna.letters))\n  if(is.null(subs.mat))\n    subs.mat <- if(looks.like.dna)dna.identity else \"BLOSUM62\"\n  if(mode(subs.mat)==\"character\"){\n    ex <- substitute(data(M,package=\"Biostrings\"),list(M=subs.mat))\n    eval(ex)\n    subs.mat <- get(subs.mat)\n  }\n  print(subs.mat)\n  N <- length(d)\n  m <- matrix(0,nrow=N,ncol=N,dimnames=list(names(d),names(d)))\n  for(i in 1:N)for(j in 1:i){\n    seqs <- sapply(strsplit(c(d[i],d[j]),split=''),c)\n    entry <- try(apply(seqs,1,function(x)subs.mat[x[1],x[2]]))\n    if(class(entry)==\"try-error\"){\n      print(seqs)\n      stop(\"Sequence difference matrix construction failed.\")\n    }\n    ## subscript out of bounds here usually means bad matrix\n    m[i,j] <- m[j,i] <- -sum(entry)\n  }\n  m <- m-min(m)\n  attr(m,'seqs') <- d\n  print(m[1:min(5,nrow(m)),1:min(5,ncol(m))])\n  m\n### The matrix of distances between each input sequence, with dimnames\n### corresponding to either the sequences, or the sequence names (if\n### they exist)\n}", 
        description = "Calculate pairwise differences between sequences using a\nsubstitution matrix.", 
        `item{seq.vec}` = "DNA or protein sequences.", `item{subs.mat}` = "Substitution matrix with dimnames that match the letters used in\nthe sequence data, or a character vector that specifies a common\nsubstitution matrix (as defined in the Biostrings package). NULL\nspecifies that we will guess a suitable substitution matrix to\nmatch your input sequences (DNA=>identity, protein=>BLOSUM62).", 
        value = "The matrix of distances between each input sequence, with dimnames\ncorresponding to either the sequences, or the sequence names (if\nthey exist)"), 
    sublogo = list(definition = "sublogo <- function\n### Draw a sublogo dendrogram for a sequence motif.\n(seqs,\n### Character vector of DNA or protein sequences (optionally named\n### with sequence labels).\n mat=NULL,\n### Substitution matrix passed to seqs.to.mat.\n ...\n### Other arguments to pass to sublogo.dendrogram (see that help page\n### for a full description).\n ){\n  sublogo.dendrogram(seqs.to.mat(seqs,mat),...)\n}",
      format="",title="sublogo",
        description = "Draw a sublogo dendrogram for a sequence motif.", 
        `item{seqs}` = "Character vector of DNA or protein sequences (optionally named\nwith sequence labels).", 
        `item{mat}` = "Substitution matrix passed to seqs.to.mat.",
        `item{\\dots}` = "Other arguments to pass to sublogo.dendrogram (see that help page\nfor a full description)."), 
    sublogo.dendrogram = list(definition = "sublogo.dendrogram <- function\n### Main function for drawing sublogo dendrograms.\n(M,\n### difference matrix as constructed by seqs.to.mat (although in\n### principle any object with a valid as.dist method could be used)\n main='',\n### plot title\n subtit=NULL,\n### plot subtitle\n base=tempfile(),\n### base file name for temporary logo files\n cutline=150,\n### Distance for cutting the tree. Draw a sublogo for each\n### leaf. Normally you will plot once, then inspect the dendrogram to\n### determine which is a good value for cutline, then plot again using\n### your chosen cutline.\n dend.width=30,\n### Percent of the plot to be occupied by the dendrogram. The logos\n### will be displayed with equal widths on either side.\n cex=1\n### character expansion factor for the dendrogram\n ){\n  hc <- hclust(as.dist(M),method=\"average\")\n  dend <- as.dendrogram(hc)\n  fam <- cutree(hc,h=cutline)[labels(dend)] # order by plotting method\n  famids <- unique(fam)\n  famtab <- data.frame(fam,y=1:length(fam))\n  famtab$seq <- attr(M,'seqs')[rownames(famtab)]\n  fam.nontriv <- sapply(famids,function(i)sum(famtab$fam==i))>1\n  names(fam.nontriv) <- famids\n  if(is.null(subtit))\n    subtit <- paste(nrow(M),\"sequences,\",sum(fam.nontriv),\"families\")\n  xrange <- c(0,1)\n  yrange <- c(0,max(famtab[,'y'])+1)\n  draw.box <- function(i){ # replace with logos\n    subtab <- famtab[famtab[,'fam']==i,]\n    grprange <- range(subtab[,'y'])\n    rect(xrange[1],grprange[1]-0.25,xrange[2],grprange[2]+0.25)\n  }\n  draw.logo <- function(i){\n    ## do not draw sublogo for a family of trivial size\n    if(fam.nontriv[as.character(i)]){ \n      subtab <- famtab[famtab[,'fam']==i,]\n      grprange <- range(subtab[,'y'])\n      tmpfile <- paste(base,i,sep='.')\n      logo <- make.logo.ps(subtab$seq,tmpfile)\n      grid.picture(logo,\n                   x=xrange[1],\n                   y=unit(grprange[1]+0.25,'native'),\n                   height=unit(diff(grprange)-0.5,'native'),\n                   distort=TRUE,\n                   just=c(0,0),\n                   fillText=TRUE)\n    }\n  }\n\n  bottomspace <- 0.6\n  topspace <- 0.4\n  side.percents <- (100-dend.width)/2\n  layout(matrix(1:3,ncol=3),c(side.percents,dend.width,side.percents))\n  par(mai=c(bottomspace,0,topspace,0),cex=1)\n  \n  ## Big summary logo on left\n  plot(xrange,yrange,\n       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')\n  biglogo <- make.logo.ps(famtab$seq,paste(base,'0',sep='.'))\n  vps <- baseViewports()\n  pushViewport(vps$inner,vps$figure,vps$plot)\n  grid.picture(biglogo,\n               y=unit(0,'native'),\n               height=unit(length(fam),'native'),\n               just=c(0.5,0),\n               exp=0,\n               distort=TRUE,\n               fillText=TRUE)\n  popViewport(3)\n  \n  ## Dendrogram in middle\n  par(mai=c(bottomspace,0,topspace,\n        max(strwidth(colnames(M),'inches',\n                     cex))/6*5),\n      family='mono')\n  par(xpd=NA)\n  plot(dend,h=TRUE,edgePar=list(lwd=2),nodePar=list(lab.cex=cex,pch=\"\"))\n  par(family=\"\")\n  segments(cutline,1,cutline,length(fam))\n  ##axis(3,cutline,lty=0,line=0)\n  \n  ## Title in the middle\n  title(main,line=0.5)\n  mtext(subtit,line=-0.7,cex=1)\n\n  # Sublogos on the right side\n  par(mai=c(bottomspace,0,topspace,0))\n  plot(xrange,yrange,\n       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')\n  vps <- baseViewports()\n  pushViewport(vps$inner,vps$figure,vps$plot)\n  sapply(famids,draw.logo)\n  popViewport(3)\n  \n  hc\n### The dendrogram from the call to hclust\n}",
      format="",title="sublogo dendrogram",
        description = "Main function for drawing sublogo dendrograms.", 
        `item{M}` = "difference matrix as constructed by seqs.to.mat (although in\nprinciple any object with a valid as.dist method could be used)", 
        `item{main}` = "plot title", `item{subtit}` = "plot subtitle", 
        `item{base}` = "base file name for temporary logo files", 
        `item{cutline}` = "Distance for cutting the tree. Draw a sublogo for each\nleaf. Normally you will plot once, then inspect the dendrogram to\ndetermine which is a good value for cutline, then plot again using\nyour chosen cutline.", 
        `item{dend.width}` = "Percent of the plot to be occupied by the dendrogram. The logos\nwill be displayed with equal widths on either side.", 
        `item{cex}` = "character expansion factor for the dendrogram", 
        value = "The dendrogram from the call to hclust"))
