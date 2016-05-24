print.qrtpcr <- function(x, ...){
  print(summary(as.data.frame(x[c("RQ")])),...)
  invisible(x)
}

.normalize.qrtpcr <- function(x,...){
  treats = as.numeric(names(table(x$Treat)))
  samples = as.numeric(names(table(x$Sample)))
  for(i in treats){
    if(i>0){
      for(j in samples){
        sele = x$Treat == i & x$Sample == j
        y = x$ddCt[sele]
        y = y[!is.na(y)]
        fy = density(y)
        ny = which(fy$y == max(fy$y))[1]
        x$ddCt[sele] = x$ddCt[sele] - fy$x[ny]
      }
    }
  }
  x$RQ = 2^x$ddCt
  x
}

.denplot <- function(x){
  x = x[!is.na(x)]
  plot(density(x))
  abline(v=0, col='gray')
}

.plot.qrtpcr <- function(x, ...){
  INDEX = paste("T",x$Treat,"-S", x$Sample,sep='')
  tapply(x$ddCt, INDEX, FUN=.denplot)
  cat("\nPlots in order:", levels(as.factor(INDEX)),"\n")
  invisible(x)
}

.validate <- function(x,y,treat=1,sample='all', ...)
UseMethod("validate")

.validate.default <- function(x,y,treat=1,sample='all', ...){
  ##  to look up the qRT-PCR results for one given gene
  stopifnot(is.character(x))
  stopifnot(length(x)==1)
  if(missing(y)) {
    ##    data(qrtpcr);
    ##    if(!is.null(qpcr)) y=normalize(qpcr);
    warning("Default qRT-PCR results in 'data(qrtpct)' are used. ")
  }
  stopifnot(class(y)=='qrtpcr')
  ## get simple gene names from a possibly compound gene names such as
  ## in Exiqon microarrays.  Then find the matchs.
  gname = x
  rtnames = y$Gene
  if(is.null(rtnames))
    stop("qRT-PCR data file not found!")
  out = .Exiqon.GeneName(gname)
  ##  if(!is.null(out$parameter)){
  out2 = match(rtnames, out$Genes)
  sele0 = which(!is.na(out2))
  sele1 = y$Treat == treat
  if(missing(sample)) sele1 = sele1
  if(tolower(sample=='all')){
    sele1 = sele1
  }else{
    sele1 = sele1 & (y$Sample==sample)
  }
  sele2 = rep(FALSE, length(sele1))
  sele2[sele0] = TRUE
  sele = sele1 & sele2
  if(sum(sele)>0){
    Sample=y$Sample[sele]
    Treat = y$Treat[sele]
    dCt=y$dCt[sele]
    RQ=y$RQ[sele]
    pv = t.test(log2(RQ))$p.value
  }else{
    Sample = NA
    Treat = NA
    dCt=NA
    RQ=NA
    pv=NA
  }
  structure(list(
                 Gene = gname,
                 Sample = Sample,
                 Treat = Treat,
                 dCt = dCt,
                 RQ = RQ,
                 p.value=pv),
            class='rtprofile')
}

  
.print.rtprofile <- function(x,...){
  cat("\nTreatment: ", x$Treat[1],"\t")
  cat("\tp-value: ", x$p.value,"\n\n")
  y = data.frame(Sample=x$Sample, dCt=x$dCt,RQ=x$RQ)
  print(y)
  invisible(x)
}

.boxplot.rtprofile <- function(x,...){
  y = x$RQ
  boxplot(y,...)
}

.Exiqon.GeneName <- function(gname){
  ##  Designed for Exiqon data.  It can be used for any string(s).
  gname=as.character(gname);
  n=nchar(gname);
  out = NULL # to record the positions of the first "-" and all "/"s
  for(i in 1:n){
    pos=n+1-i;
    if(substr(gname,pos,pos)=="/"){
      out = c(out, pos)
    }
  }
  genes = gname;
  if(!is.null(out)){
    out = rev(out)
    for(i in 1:out[1]){
      pos=out[1]+1-i;
      if(substr(gname,pos,pos)=="-"){
        out = c(pos, out)
        break;
      }
    }
    genes = substr(gname,1,out[2]-1);
    chead = substr(gname,1,out[1]);
    j = length(out)
    out2 = c(out, n+1)
    for(i in 2:j){
      tmp = substr(gname, out2[i]+1, out2[i+1]-1)
      genes = c(genes, paste(chead, tmp,sep=''))
    }
  }
  list(Original=gname, parameters=out, Genes = genes);
}


