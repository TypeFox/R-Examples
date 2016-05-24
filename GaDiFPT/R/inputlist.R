inputlist <- 
  function(m,s,Sty,tini,xini,Tend,delta,Ntime,quadfl,RSfl)
  {
#   out <-  as.list(c(m,s,Sty,tini,xini,Tend,delta,Ntime,quadfl,RSfl))
    out <-  list(m,s,Sty,tini,xini,Tend,delta,Ntime,quadfl,RSfl)
   attr(out,"names") <- c("mu","sigma2","Stype","t0","x0","Tfin","deltat","M","quadflag","RStudioflag") 
   class(out) <- c("inputlist", "list")
   return(out)
  }



