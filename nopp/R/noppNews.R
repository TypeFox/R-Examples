noppNews <- function () 
{
    newsfile <- file.path(system.file(package = "nopp"), 
        "NEWS")
    file.show(newsfile)
}
