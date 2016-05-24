## page related functions

openPage=function(filename, dirname=NULL, title=filename, link.javascript=NULL,
  link.css=NULL, css=NULL, head=NULL, charset="utf-8", lang="en",
  head.attributes=NULL, body.attributes=NULL) {
  if (!is.null(dirname)) {
    if (!file.exists(dirname)) dir.create(dirname, recursive=TRUE, showWarnings=FALSE)
    filename = file.path(dirname, filename)
  }
  page = file(filename,'wt')
  doctype = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'
  meta = hmakeTag('meta',NULL,'http-equiv'='Content-Type',content=paste("text/html; charset=", charset, sep=''), newline=FALSE)
  
  if (!is.null(link.javascript)) link.javascript = paste(hmakeTag('script', language='JavaScript', src=link.javascript), collapse='\n')
  if (!is.null(link.css)) link.css = paste(hmakeTag('link', rel='stylesheet', type='text/css', href=link.css), collapse='\n')
  if (!is.null(css)) css = paste(hmakeTag('style', css), collapse='\n')
  
  head = paste(meta, hmakeTag('title',title), head, link.javascript, link.css, css, sep='\n')
  head = do.call(hmakeTag, c(list('head', head, newline=TRUE), head.attributes))
  bodyStart = do.call(hmakeTag, c(list('body', NULL), body.attributes))
  bodyStart = substr(bodyStart, 1, regexpr('</body>', bodyStart)-1)
  hwrite(paste(doctype, "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='", lang, "' lang='", lang, "'>", head, bodyStart, sep=''), page)
  page
}

getHwriterVersion=function() {
  (sessionInfo()$otherPkgs)[['hwriter']]$Version
}

closePage=function(page, splash=TRUE) {
  hwriterlink = hwrite('hwriter', link='http://www.embl.de/~gpau/hwriter/index.html')
  if (splash) hwrite(paste('\n<br/><br/><font size=\"-2\">(Page generated on ', date(), ' by ', hwriterlink, ' ', getHwriterVersion(), ')</font>', sep=''), page, br=TRUE)
  else hwrite('\n<br/><br/>', page, br=TRUE)
  hwrite('</body></html>', page, br=FALSE)
  close(page)
}
