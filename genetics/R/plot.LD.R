# $Id: plot.LD.R 150 2003-06-04 21:22:57Z warnesgr $

plot.LD.data.frame <- function(x,
                               digits=3,

                               colorcut=c(0,0.01, 0.025, 0.5, 0.1, 1),
                               colors=heat.colors(length(colorcut)),
                               textcol="black",

                               marker,
                               which="D'",
                               distance,
                               ...)
  {
    oldpar <- par("mfrow")
    
    par(mfrow=c(1,2))

    LDtable(x, digits=digits, colorcut=colorcut, colors=colors,
            textcol=textcol, ...)
    LDplot(x, marker=marker, which=which, distance=distance, ...)
    
    par(mfrow=oldpar)
    invisible()
  }


LDtable <- function(x, 
                    colorcut=c(0,0.01, 0.025, 0.5, 0.1, 1),
                    colors=heat.colors(length(colorcut)),
                    textcol="black",
                    digits=3,
                    show.all=FALSE,
                    which=c("D", "D'", "r", "X^2", "P-value", "n"),
                    colorize="P-value",
                    cex,
                    ...)
  {
    if(! colorize %in% names(x))
      stop(colorize, " not an element of ", deparse(substitute(x)) )

    datatab <- summary(x)
    
    missmatch <- which[!(which %in% names(x))]
    if(length(missmatch)>0)
      stop(missmatch, " not an element of ", deparse(substitute(x)) )

    matform <- function( value, template )
      {
        dim(value) <- dim(template)
        dimnames(value) <- dimnames(template)
        value
      }
    
    tmp <- cut(x[[colorize]], colorcut, include.lowest=TRUE)
    colormat <- matform(as.numeric(tmp), x[[colorize]] )
    n <- matform( paste("(",x$n,")",sep="" ), x$n)

    if(!show.all)
      { # remove blank row/column
        colormat <- colormat[-nrow(colormat),-1, drop=FALSE]
        n <- n[-nrow(n),-1, drop=FALSE]
      }

    #
    # color coded frame boxes
    #
    image(x=1:ncol(colormat), y=1:ncol(colormat),
          z=t(colormat[nrow(colormat):1,]),
          col=colors, xlab="Marker 2\n\n", ylab="Marker 1",
          xaxt="n", yaxt="n",...)
    
    abline(v=-0.5 + 1:(ncol(colormat)+1))
    abline(h=-0.5 + 1:(nrow(colormat)+1))
    
    axis(3, 1:ncol(colormat), colnames(colormat) )
    axis(2, 1:nrow(colormat), rev(rownames(colormat)) )
    
    #
    # text in boxes
    #
    cex.old <- par("cex")

    if(missing(cex))
      cex <-min( c(1/10, 1/(length(which)+1 ) ) /
                 c(strwidth("W"), strheight("W")*1.5))
    
    par(cex=cex)

    lineheight <- strheight("W")*1.5
    center <- lineheight * length(which)/2
    
    for(i in 1:length(which))
      {
        displaymat <- x[[which[i]]]

        if(!show.all)
          displaymat <- displaymat[-nrow(displaymat),-1, drop=FALSE]
        
        if( which[i]=="P-value" )
          displaymat <- format.pval(displaymat, digits=digits)
        else if (which[i]!="n")
          displaymat <- format(displaymat, digits=digits)

        displaymat[] <- gsub("NA.*", "", as.character(displaymat))
        
        text(x=col(colormat),
             y=nrow(colormat) - row(colormat)+ 1 + center - lineheight*(i-1),
             displaymat,
             col=textcol,
             adj=c(0.5, 1)
             )
      }

    text(x=1, y=1, paste(which, collapse="\n"), adj=c(0.5,0.5) )

    par(cex=cex.old)
    
    #
    # title
    #
    title(main="Linkage Disequilibrium\n")

    invisible(colormat)
  }



LDplot <- function(x, 
                   digits=3,
                   marker,
                   distance,
                   which=c("D", "D'", "r", "X^2", "P-value", "n", " "),
                   ...)
{
  which = match.arg(which)
  
  if(missing(marker))
    marker <- colnames(x[[which]])
  else if (is.numeric(marker))
    marker <- colnames(x[[which]])[marker]
  
  datamat <- ifelse( is.na(x[[which]]), t(x[[which]]), x[[which]])

  if(which %in% c("D'","r") )
    diag(datamat) <- 1.0
  else if (which=="P-value")
    diag(datamat) <- 0.0
  
  dimnames(datamat) <- dimnames(x[[which]])
  
  if(missing(distance)) distance <- 1:ncol(datamat)
  distance <- matrix(distance, ncol=ncol(datamat), nrow=nrow(datamat),
                     byrow=TRUE)
  dimnames(distance) <- dimnames(datamat)
  
  matplot(x=t(distance[marker,,drop=FALSE]),
          t(datamat[marker,,drop=FALSE]),
          type="b", 
          xlab="Marker",
          ylab=paste("Linkage Disequilibrium: ", which, sep=""),
          xaxt="n",
          ... )
  
  axis(1, distance[1,], paste(1:ncol(datamat), colnames(datamat), sep=": " ))
  
  title("Pairwise Disequilibrium Plot")

  invisible()
}
