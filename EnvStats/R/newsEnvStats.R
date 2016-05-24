newsEnvStats <-
function () 
{
    newsfile <- file.path(system.file(package = "EnvStats"), 
        "NEWS")
    file.show(newsfile)
}
