chron_trans <-
function(format = "%Y-%m-%d", n = 5)
{
    breaks. <- function(x)
        chron(scales::pretty_breaks(n)(x))
    format. <- function(x)
        format(as.POSIXct(x, tz = "GMT"), format = format)
    scales::trans_new("chron",
                      transform = as.numeric, inverse = chron,
                      breaks = breaks., format = format.)
}

scale_x_chron <-
function(..., format = "%Y-%m-%d", n = 5)
{
    ggplot2::scale_x_continuous(..., trans = chron_trans(format, n))
}

scale_y_chron <-
function(..., format = "%Y-%m-%d", n = 5)
{
    ggplot2::scale_y_continuous(..., trans = chron_trans(format, n))
}
