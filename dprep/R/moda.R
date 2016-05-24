moda <-
structure(function (x, na.rm = TRUE) 
{
    if (na.rm == TRUE) 
        m1 = rev(sort(table(x[])))
    else m1 = rev(sort(table(x, exclude = NULL)))
    moda = names(m1[m1 == m1[1]])
    if (is(x, "numeric")) 
        moda = as.numeric(moda)
    return(moda)
}, source = c("function(x,na.rm=TRUE)", "{", "  ", "#Function that finds the mode of vector x", 
"", "  if(na.rm==TRUE) m1=rev(sort(table(x[])))", "    else m1=rev(sort(table(x,exclude=NULL)))", 
"  moda=names(m1[m1==m1[1]])", "  if (is(x,\"numeric\")) moda=as.numeric(moda)", 
"  return(moda)", "}"))
