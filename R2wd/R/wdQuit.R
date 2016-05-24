wdQuit <-
function (wdapp = .R2wd)
{
    wdapp$Quit()
    rm(.R2wd,pos=1)
}

