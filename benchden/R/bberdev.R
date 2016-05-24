ben01 <- function()  {return(c(0,1))}
ben02 <- function()  {return(0)}
ben03 <- function()  {return(0)}
ben04 <- function()  {return(0)}
ben05 <- function()  {return(NULL)}
ben06 <- function()  {return(NULL)}
ben07 <- function()  {return(NULL)}
ben08 <- function()  {return(c(0,1))}
ben09 <- function()  {return(1)} 
ben10 <- function()  {return(0)}
ben11 <- function()  {return(NULL)}
ben12 <- function()  {return(0)}
ben13 <- function()  {return(c(-5,-0.5,0.5,5))}
ben14 <- function()  {return(c(-1/exp(2),0,1/exp(2)))}
ben15 <- function()  {return(c(0,1))}
ben16 <- function()  {return(c(-1,0,1))}
ben17 <- function()  {return(c(0,1))}
ben18 <- function()  {return(0)}
ben19 <- function()  {return(NULL)}
ben20 <- function()  {return(0)}
ben21 <- function()  {return(NULL)}
ben22 <- function()  {return(NULL)}
ben23 <- function()  {return(NULL)}
ben24 <- function()  {return(NULL)}
ben25 <- function()  {return(c(-1.1,-0.1,0.1,1.1))}
ben26 <- function()  {return(c(-20.1,-20,-1,1,20,20.1))}
ben27 <- function()  {return(seq(-10,10))}
ben28 <- function()  {return(c(0,1))}

`bberdev` <- function(dnum=1) {
    if (is.nan(dnum) || ! dnum %in% 1:28)
        stop("dnum must be between 1 and 28")
    return (
            eval(               
                parse(text = sprintf("ben%02d()", dnum)) # evaluate "ben[dnum]()"-string
            )
        )
}
