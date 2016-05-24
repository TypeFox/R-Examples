month_cnv=function(monthinput, up=FALSE, short=FALSE){

  numelem = length(monthinput)
  monthnames = c('january', 'february', 'march', 'april', 'may', 'june',
    'july', 'august', 'september', 'october', 'november', 'december')
  monthshort = substr(monthnames,1,3)

  if (is.character(monthinput)){
    monthinput = gsub('(^ +)|( +$)','',monthinput)
    shortinput = tolower(substr(monthinput,1,3))
    return(which(shortinput==monthshort))
  }
  else {
    if(any(monthinput<1 | monthinput>12)) {
      stop("bad input values.  month numbers must be 1-12.")
    }
    else {
      result = monthnames[monthinput]
      if(short) result = substr(result,1,3)
      if(up) result = toupper(result)
    }
  }
  return(result)
}   #   function month_cnv

