"list2vec" <-
function(liste)
{
    # function that creates vector with numerical elements of list
    # first element of list needs to be numeric

    # error control
## Wird die folgende Abfrage gebraucht ?
##    if (is.null(liste)) return(NULL)
##
##
##    else 
##    {

    if (!is.list(liste) || !is.numeric(liste[[1]])) 
    stop(" the argument for list2vec needs to be a list with first element numeric.")
    ausgabe <- liste[[1]]
    if (length(liste) > 1)
    {
    for (a in 2:length(liste))
    {
    if (is.numeric(liste[[a]])) ausgabe<-c(ausgabe,as.numeric(liste[[a]]))
    }
    }
    return(ausgabe)
##   } ## falls die oben die vorgeschaltete if-Abfrage 
}

