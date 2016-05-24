
htmlify <- function (x) # Taken from the tools package
{
    fsub <- function(pattern, replacement, x) 
       gsub(pattern, replacement, x, fixed=TRUE, useBytes=TRUE)
   
    x <- fsub("&", "&amp;", x)
    x <- fsub("---", "&mdash;", x)
    x <- fsub("--", "&ndash;", x)
    x <- fsub("``", "&ldquo;", x)
    x <- fsub("''", "&rdquo;", x)
    x <- gsub("`([^']+)'", "&lsquo;\\1&rsquo;", x, perl=TRUE, useBytes=TRUE)
    x <- fsub("`", "'", x)
    x <- fsub("<", "&lt;", x)
    x <- fsub(">", "&gt;", x)
    x <- fsub("\"\\{\"", "\"{\"", x)
    x <- fsub("\"", "&quot;", x)
    x
}

htmlNumeric <- function(chars, minus=TRUE, leftpad=TRUE, rightpad=TRUE) {
    regexp <- "^( *)([-]?)([^ -][^ ]*)( *)$"
    leadin <- sub(regexp, "\\1", chars)
    sign <- sub(regexp, "\\2", chars)
    rest <- sub(regexp, "\\3", chars)
    tail <- sub(regexp, "\\4", chars)
    
    figurespace <- "&#x2007;"
    minussign <- "&minus;"
    
    if (minus && any(neg <- sign == "-")) {
    	if (any(leadin[!neg] == ""))
    	    leadin <- sub("^", " ", leadin)
    	leadin[!neg] <- sub(" ", "", leadin[!neg])
    	sign[!neg] <- figurespace
    	sign[neg] <- minussign
    }
    if (leftpad && any(ind <- leadin != "")) 
    	leadin[ind] <- gsub(" ", figurespace, leadin[ind])
    	
    if (rightpad && any(ind <- tail != ""))
    	tail[ind] <- gsub(" ", figurespace, tail[ind])
    	
    paste(leadin, sign, rest, tail, sep="")
}

CSSclassname <- function(just) 
    ifelse(just == "l", "left",
    ifelse(just == "c", "center",
    ifelse(just == "r", "right", just)))

html.tabular <- function(object, file="", 
                         options=NULL, id=NULL, ...) {
    if (is.character(file)) {
	if (file == "")
    	    out <- ""
    	else {
    	    out <- file(file, open="wt")
    	    on.exit(close(out))
    	}
    } else
    	out <- file
    	
    if (!is.null(options)) {
    	saveopts <- do.call(table_options, options)
    	on.exit(table_options(saveopts), add=TRUE)
    }
    opts <- table_options()
    
    mycat <- function(...) cat(..., file=out)
    
    defjust <- opts$justification
    blankhead <- "  <th>&nbsp;</th>\n"

    classes <- chars <- format(object, html = TRUE, minus = opts$HTMLminus, 
                               leftpad = opts$HTMLleftpad, 
                               rightpad = opts$HTMLrightpad, ...) # format without justification
    classes[] <- ""
    
    vjust <- attr(object, "justification")
    vjust[is.na(vjust)] <- defjust
    ind <- vjust != defjust
    classes[ind] <- sprintf(' class="%s"', CSSclassname(vjust[ind]))
    chars[chars == ""] <- "&nbsp;"
    chars[] <- sprintf("  <td%s>%s</td>\n", classes, chars)
    
    rowClasses <- rowLabels <- attr(object, "rowLabels")
    rowClasses[] <- ""
    nleading <- ncol(rowLabels)
    rowLabels[is.na(rowLabels)] <- "&nbsp;"
    rjust <- attr(rowLabels, "justification")
    rjust[is.na(rjust)] <- opts$rowlabeljustification
    ind <- rjust != defjust
    rowClasses[ind] <- sprintf(' class="%s"', CSSclassname(rjust[ind]))
    rowLabels[] <- sprintf( "  <th%s>%s</th>\n", rowClasses, rowLabels)
    
    colnamejust <- attr(rowLabels, "colnamejust")
    colnamejust <- rep(colnamejust, length.out=nleading)
    colnameClasses <- colnames(rowLabels)
    colnameClasses[] <- ""
    ind <- is.na(colnamejust)
    colnamejust[ind] <- defjust
    ind <- colnamejust != defjust
    colnameClasses[ind] <- sprintf(' class="%s"', CSSclassname(colnamejust[ind]))
    colnames(rowLabels) <- sprintf("  <th%s>%s</th>\n", colnameClasses, colnames(rowLabels))
    
    clabels <- attr(object, "colLabels")
    cjust <- attr(clabels, "justification")
    ind <- is.na(cjust)
    cjust[ind] <- defjust

    multi <- matrix(0, nrow(clabels), ncol(clabels))
    prevmulti <- rep(0, nrow(multi))
    for (i in rev(seq_len(ncol(multi)))) {
    	ind <- is.na(clabels[,i])
    	multi[!ind, i] <- 1 + prevmulti[!ind]
    	prevmulti[ind] <- 1 + prevmulti[ind]
    	prevmulti[!ind] <- 0
    }
    colspan <- ifelse(multi < 2, "", sprintf(' colspan="%d"', multi))
    class <-   ifelse(cjust == defjust || multi == 0, "", sprintf(' class="%s"', CSSclassname(cjust)))
    clabels[clabels == ""] <- "&nbsp;"
    clabels <- ifelse(multi == 0, "", sprintf('  <th%s%s>%s</th>\n', colspan, class, clabels))
    
    rowLabelHeadings <- matrix(blankhead, nrow(clabels), ncol(rowLabels))
    rowLabelHeadings[nrow(clabels),] <- colnames(rowLabels)
    
    if (opts$doHTMLheader) {
    	head <- sub("CHARSET", localeToCharset(), opts$HTMLhead, fixed=TRUE)
     	mycat(head)
    }
     	
    if (opts$doCSS) {
        if (is.null(id)) 
            css <- gsub("#ID ", "", opts$CSS, fixed=TRUE)
        else
            css <- gsub("#ID", paste0("#", id), opts$CSS, fixed=TRUE)
     	mycat(css)
    }
    
    if (opts$doHTMLbody)
     	mycat(opts$HTMLbody)
   
    if (opts$doBegin) {
    	if (is.null(id)) 
    	    id <- ""
        else
    	    id <- sprintf(' id="%s"', id)
        mycat(sprintf('<table%s %s>\n', id, opts$HTMLattributes))
    }
    if (!is.null(opts$HTMLcaption))
        mycat(sprintf('<caption>%s</caption>\n', opts$HTMLcaption))
	
    if (opts$doHeader) {
        rows <- apply(cbind(rowLabelHeadings, clabels), 1, paste0, collapse="")
        mycat('<thead>\n')
	mycat(sprintf('<tr class="%s">\n%s</tr>\n', CSSclassname(defjust), rows))
	mycat('</thead>\n')
    }
    if (opts$doFooter && !is.null(opts$HTMLfooter)) {
        mycat('<tfoot>\n')
        mycat(opts$HTMLfooter)
        mycat('</tfoot>\n')
    }
    if (opts$doBody) {
	rows <- apply(cbind(rowLabels, chars), 1, paste0, collapse="")
	mycat('<tbody>\n')
	mycat(sprintf('<tr class="%s">\n%s</tr>\n', CSSclassname(defjust), rows))
	mycat('</tbody>\n')
    }
    if (opts$doEnd)
    	mycat("</table>\n")
    invisible(structure(list(file=file), class="html"))
}
