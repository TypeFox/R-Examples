print.fedreg_document <- function(x, ...){
    
    cat(x$type)
    if(!is.null(x$subtype))
        cat(', ', x$subtype)
    if(!is.null(x$toc_doc))
        cat(': ', x$toc_doc, sep='')
    else if(!is.null(x$title))
        cat(': ', x$title, sep='')
    cat('\n')
    cat('FR Document Number:',x$document_number)
    if(!is.null(x$citation))
        cat(' Citation: ', x$citation)
    cat('\n')
    if(!is.null(x$publication_date))
    cat('Date:', x$publication_date, '\n')
    if(!is.null(x$page_length))
        cat(x$page_length, 'Pages:', x$start_page, '-', x$end_page, '\n')
    
    if(!is.null(x$agencies)){
        cat('Agencies mentioned:\n')
        lapply(x$agencies, function(a) {
            if(is.list(a)) cat('  ', a$raw_name, '(',a$id,')\n',sep='')
        })
    }
    
    cat('URLs:\n')
    if(!is.null(x$raw_text_url))
        cat('  Raw text: ', x$raw_text_url, '\n')
    if(!is.null(x$body_html_url))
        cat('  HTML:     ', x$body_html_url, '\n')
    else if(!is.null(x$html_url))
        cat('  HTML:     ', x$html_url, '\n')
    if(!is.null(x$json_url))
        cat('  JSON:     ', x$json_url, '\n')
    cat('\n')
    invisible(x)
}
