extract = function (data, direction = c("time", "x"), timeorder, xorder) 
{
    direction = match.arg(direction)
    if (direction == "time") {
        output = extract.time(data = data, timeorder = timeorder)
    }
    else {
        output = extract.x(data = data, xorder = xorder)
    }
    return(output)
}