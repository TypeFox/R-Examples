`jadjust.length` <-
function (inputdata) 
{
    if (is.character(inputdata)) 
        s <- scan(inputdata)
    else s <- inputdata
    np <- length(s)
    pow <- 1
    while (2 * pow < np) pow <- 2 * pow
    new.np <- 2 * pow
    if (np == new.np) 
        list(signal = s, length = np)
    else {
        new.s <- 1:new.np
        new.s[1:new.np] <- 0
        new.s[1:np] <- s
        list(signal = new.s, length = new.np)
    }
}
