TempEnv <- function ()
{
    pos <-  match("SciViews:TempEnv", search())
    if (is.na(pos)) { # Must create it
        `SciViews:TempEnv` <- list()
        Attach <- function (...) get("attach", mode = "function")(...)
        Attach(`SciViews:TempEnv`, pos = length(search()) - 1)
        rm(`SciViews:TempEnv`)
        pos <- match("SciViews:TempEnv", search())
    }
    pos.to.env(pos)
}
