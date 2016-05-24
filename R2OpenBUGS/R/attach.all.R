attach.all <- function(x, overwrite = NA, name = "attach.all"){
    rem <- names(x) %in% ls(.GlobalEnv)
    if(!any(rem)) overwrite <- FALSE
    rem <- names(x)[rem]
    if(is.na(overwrite)){
        question <- paste("The following objects in .GlobalEnv will mask\nobjects in the attached database:\n", 
                        paste(rem, collapse=", "), 
                        "\nRemove these objects from .GlobalEnv?", sep="")
        if(interactive()){
            if(.Platform$OS.type == "windows")
                overwrite <- "YES" == winDialog(type = "yesno", question)
            else
                overwrite <- 1 == menu(c("YES", "NO"), graphics = FALSE, title = question)
        }
        else overwrite <- FALSE
    }
    if(overwrite) remove(list=rem, envir=.GlobalEnv)
    attach(x, name=name)
}

attach.bugs <- function (x, overwrite = NA){
  if(class(x) != "bugs")
    stop("attach.all() requires a bugs object.")
  if("bugs.sims" %in% search()){
    detach("bugs.sims")}
  x$sims.list$n.sims <- x$n.sims   # put n.sims into sims.list for convenience
  r2 <- attach.all(x$sims.list, overwrite = overwrite, name = "bugs.sims")
  invisible (r2)
}

detach.all <- function(name = "attach.all"){
    do.call("detach", list(name=name))
}

detach.bugs <- function(){
    detach.all("bugs.sims")
}
