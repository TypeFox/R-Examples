printparms <-
function() {
curparms<-"Make me local" ## Probably only relevant to the checker, which hasn't checked the "load" file for what it contains, not real locality.
cfile<-paste0(tempdir(),"/","curparms")
trial<-try(load(cfile))
if (class(trial)=="try-error") {cat(paste0("No parameters saved yet -- presumably you haven't had a succesful call of phyreg() yet\n"))} else
print(curparms)
}
