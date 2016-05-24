shift.word.indices <-
function (documents, amount) 
{
    lapply(documents, function(x) {
        x[1, ] <- x[1, ] + amount
        x
    })
}
