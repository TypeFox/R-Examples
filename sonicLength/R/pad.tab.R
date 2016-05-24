pad.tab <- function(tab, start.at=NULL,end.at=NULL)
  {
    ## Purpose: create a padded table of quasi.counts
    ## ----------------------------------------------------------------------
    
    ## Arguments: tab - named vector of non-negative entries, with
    ##                    !any(is.na(as.numeric(as.character(names(tab)))))
    ##            start.at - pad zeroes from what (named) value
    ##            end.at - pad zeroes to what value
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 20 Feb 2011, 16:19
    
    if (is.null(start.at)) start.at <- 1L
    if (is.null(end.at)) end.at <- max(as.numeric(as.character(names(tab))))+10
    
    new.tab.names <- start.at:end.at
    new.tab <- 0.0*new.tab.names
    names(new.tab) <- new.tab.names
    new.tab[ names(tab) ] <- tab
    data.frame(y=unname(new.tab),x=new.tab.names,orig=names(new.tab)%in%names(tab))
  }
