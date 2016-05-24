document.lengths <-
function (docs) 
{
    sapply(docs, function(x) sum(x[2, ]))
}
