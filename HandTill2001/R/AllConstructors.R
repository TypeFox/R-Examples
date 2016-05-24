bincap <-
function (response, predicted, true = "1") 
{
    return(new(Class = "bincap", response = response, predicted = predicted, 
        true = true))
}


multcap <-
function (response, predicted) 
{
    return(new(Class = "multcap", response = response, predicted = predicted))
}
