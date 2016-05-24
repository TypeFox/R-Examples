ColList <- function(n)
{
    n<-round(n)
    if (n <= 15) {
        col_list <- c (
            "#ee3333",                     # red
            "#3366aa",                     # blue
            "#009872",                     # green
            "#982187",                     # purple
            "#faa200",                     # yellow
            "#267aa4",                     # teal
            "#910000",                     # maroon
            "#b56cfe",                     # violet
            "#00b7ec",                     # cyan
            "#f36a18",                     # orange
            "#534731",                     # brown
            "#fdb5da",                     # pink
            "#064650",                     # bluegreen
            "#b5dafe",                     # sky
            "#000000"                      # black
        )
    } else 
        col_list <- NULL
    
    
    return(col_list)
}
