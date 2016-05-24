library(R2admb)
## I'm not sure whether this problem *really* demonstrates a race condition;
## it might be a separate effect.  I do know that it's a problem.  I don't
## know if it's

load(system.file("admbtests","racecond.RData",package="R2admb"))
file.copy(system.file("admbtests","deltavar_prof.tpl",package="R2admb"),"deltavar_prof.tpl")
d1 <- try(do_admb("deltavar_prof",
              data=list(nt=length(dvar)+1,
                        dvar=dvar,max_var_g=10),
              profile=TRUE,
              run.opts=run.control(compile=FALSE,clean_files=FALSE),
              params=svals))
file.copy("deltavar_prof.pin","deltavar_prof1.pin",overwrite=TRUE)

## works
d1B <- try(do_admb("deltavar_prof",
              data=list(nt=length(dvar)+1,
                        dvar=dvar,max_var_g=10),
              profile=TRUE,
                   extra.args="-crit 1.e-7",
                   run.opts=run.control(compile=FALSE,clean_files=FALSE),
                   params=svals))
file.copy("deltavar_prof.pin","deltavar_prof1.pin",overwrite=TRUE)

svals2 <- list(rho=0.5,var_g=4)

d2 <- try(do_admb("deltavar_prof",
              data=list(nt=length(dvar)+1,
                        dvar=dvar,max_var_g=10),
              profile=TRUE,
              run.opts=run.control(compile=FALSE,clean_files=FALSE),
              params=svals2))
file.copy("deltavar_prof.pin","deltavar_prof2.pin",overwrite=TRUE)
