attach.all <- function(x, overwrite = NA, name = "attach.all"){
  if(is.R()){
    rem <- names(x) %in% ls(.GlobalEnv)
    if(!any(rem)) overwrite <- FALSE
    rem <- names(x)[rem]
    if(is.na(overwrite)){
        question <- paste("The following objects in .GlobalEnv will mask\nobjects in the attached database:\n", 
                        paste(rem, collapse=", "), 
                        "\nRemove these objects from .GlobalEnv?", sep="")
        if(interactive()){
            if(.Platform$OS.type == "windows")
                overwrite <- "YES" == utils::winDialog(type = "yesno", question)
            else
                overwrite <- 1 == menu(c("YES", "NO"), graphics = FALSE, title = question)
        }
        else overwrite <- FALSE
    }
    if(overwrite) remove(list=rem, envir=.GlobalEnv)
    attach(x, name=name)
  } else {
    ## next line is a dirty trick for R'd codetools check in R-2.5.0
    ## (should be removed after codetoold have been improved):
    attach.default <- get("attach.default") 
    attach.default(x, name = name)
  }
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
  if (is.R()){
    do.call("detach", list(name=name))
  } else {
    do.call("detach", list(what=name))
  }
}

detach.bugs <- function(){
  if (is.R()){
    detach.all("bugs.sims")
  } else {
    invisible(detach.all("bugs.sims"))
  }
}
