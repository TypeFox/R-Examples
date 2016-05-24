## environemnt to store result and headers in
.e <- new.env(parent = emptyenv())

out <- function(..., sep='', eol='\n')
    .e$out <- c(.e$out, paste(..., sep=sep, collapse=eol))

oclear <- function(output=TRUE, headers=FALSE) {
  if (output) .e$out <- character(0)
  if (headers) .e$headers <- NULL
}

otable <- function(..., tab='', tr='', cs='</td><td>') {
  a <- list(...)
  if (length(a)==1 && is.list(a[[1]])) a <- a[[1]]
  ml <- max(unlist(lapply(a,length)))
  m <- matrix(unlist(lapply(a,function(x) rep(as.character(x),length.out=ml))),ml)
  tout <- unlist(lapply(1:ml, function(x) paste("<tr",tr,"><td>",paste(m[x,],collapse=cs),"</td></tr>",sep='')))
  .e$out <- c(.e$out, paste("<table ",tab,">\n",sep=''), tout, '</table>')
}

ohead <- function(..., level=3)
  .e$out <- c(.e$out, paste("<h",level,">",paste(...,sep=''),"</h",level,">",sep=''))

oprint <- function(..., sep='\n')
  .e$out <- c(.e$out, paste("<pre>",paste(capture.output(print(...)),collapse=sep),"</pre>",sep=''))

.opts <- function(..., disabled=FALSE) {
  l <- list(...)
  disabled <- if(isTRUE(disabled)) " disabled" else ""
  paste(
        if (length(l)) {
          n <- names(l)
          if (is.null(n) || any(n=="")) stop("Invalid unnamed argument")
          paste(unlist(lapply(seq.int(l), function(i) paste(" ",n[i],"=\"",gsub("\"","&quot;",as.character(l[[i]])[1],fixed=TRUE),"\"",sep=''))), collapse='')
        } else ""
        , disabled, sep='')
}

oselection <- function(name, text, values=text, sel.index, sel.value, size, ...) {
  if (!missing(sel.index) && !missing(sel.value)) stop("only one of 'sel.index' and 'sel.value' can be specified")
  if (!length(name)) stop("element name must be not be empty")
  name <- as.character(name)[1]
  if (missing(sel.index) && missing(sel.value)) sel.index <- integer(0)
  if (missing(sel.index)) sel.index <- !is.na(match(values, sel.value))
  size <- if (missing(size)) '' else paste(" size=\"",as.character(size)[1],"\"",sep='')
  if (!is.logical(sel.index)) sel.index <- !is.na(match(seq.int(values), sel.index))
  name <- as.character(name)[1]
  .e$out <- c(.e$out,
              paste("<select name=\"", name, "\"", size, .opts(...), ">",sep=''),
              paste("<option value=\"",gsub('"','&quot;',values,fixed=TRUE),"\"",c(""," selected")[as.integer(sel.index) + 1L],">",text,"</option>",sep='',collapse='\n'),
              "</select>")
}

oinput <- function(name, value, size, type="text", checked=FALSE, ...) {
  if (!length(name)) stop("element name must be not be empty")
  name <- as.character(name)[1]
  size <- if (missing(size)) '' else paste(" size='",as.character(size)[1],"'",sep='')
  value <- if (missing(value)) '' else paste(" value=\"",gsub('"','&quot;',as.character(value)[1]),'"',sep='')
  checked <- if (isTRUE(checked)) " checked" else ""
  .e$out <- c(.e$out, paste("<input type=\"", as.character(type)[1], "\" name=\"",name,"\"", value, size, checked, .opts(...), ">", sep=''))
}

osubmit <- function(name="submit", ...) oinput(name=name, type="submit", ...)


arequest <- function(txt, target, where, ..., attr='') {
     if (length(list(...)))
          paste("<a href='javascript:void(0);' onclick=\"javascript:req('",target,"','",where,"','",gsub("'","&#39;", paste(..., sep=''), fixed=TRUE),"');\"",attr,">",txt,"</a>",sep='')
     else
          paste("<a href='javascript:void(0);' onclick=\"javascript:req('",target,"','",where,"');\"",attr,">",txt,"</a>",sep='')
}

add.header <- function(txt) {
    if (!is.null(.e$headers))
      .e$headers <- as.character(txt)
    else
          .e$headers <- c(.e$headers, as.character(txt))
    invisible(.e$headers)
}

#link <- function(url,target,par,...) paste("<a href=# onclick=\"req('",where,"','",target,"','",par,"');\">",...,"</a>",sep='',co

# note: .e$headers are added by WebResult automatically
done <- function(..., cmd="html", type="text/html; charset=utf-8")
  WebResult(cmd, ifelse(length(list(...)), paste(as.character(.e$out),paste(...,sep='',collapse='\n'),sep='',collapse='\n'), paste(as.character(.e$out),collapse='\n')), type)

# create query string from 'pars' and merge in any additional parameters passed
.npar <- function(...) {
     # we have to use get(), otherwise codetools get confused...
     q <- if (exists("pars") && is.list(get("pars"))) get("pars") else list()
     l <- list(...)
     for (i in names(l)) q[[i]] <- if (l[[i]] == "") NULL else l[[i]]
     if (!length(q)) return("")
     ml <- max(unlist(lapply(q,length)))
     for (i in names(q)) q[[i]] <- rep(paste(i,"=",as.character(q[[i]]),sep=''),length.out=ml)
     unlist(lapply(1:ml, function(x) paste(unlist(lapply(q, function(a) a[[x]])),sep='',collapse='&')))
}

#alink <- function(text, href, ...) {
#  if (missing(href)) href <- 'javascript:void(0);'
#  a <- list(...)
#}

#arequest <- function(what, where, ...) {
#  a <- list(...)
#  xp <- ''
#  if (length(a)) paste(names(a),as.character(a),sep='=',collapse='&'
#}


## general tools

# a fast way to convert a list into a dataframe
.df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.null(x)) return(NULL)
  if (!is.list(x)) x <- as.list(x)
  attr(x,"row.names") <- c(NA_integer_, -length(x[[1]]))
  class(x) <- "data.frame"
  x
}
