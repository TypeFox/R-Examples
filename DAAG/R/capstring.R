"capstring" <-
function (names) 
{
paste(toupper(substring(names, first=1, last=1)),substring(names, first=2), sep="")
}

