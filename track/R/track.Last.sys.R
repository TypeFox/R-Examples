## Doesn't work properly to have track.stop() in .Last.sys or
## .Last in the tracking package -- need to assign .Last in
## the global environment.
## If .Last or .Last.sys exist in the track package, they
## appear to get called after saveimage() in the
## Windows GUI (if at all), which results in all objects being saved
## in .RData as well as in the tracking DB, which results
## in problems on startup because of duplicated vars.

## .Last and .Last.sys are called from R_dot_Last() in src/main/main.c
## from R_CleanUp in src/gnuwin32/system.c
## .Last is called before save.image()
## Then R_RunExitFinalizers(); is called after that.
track.Last <- function() {
    if (!isTRUE(getOption("global.track.options")$inhibit.Last))
        track.stop(all=TRUE, sessionEnd=TRUE, callFrom=".Last()")
}

