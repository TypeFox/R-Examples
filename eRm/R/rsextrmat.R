"rsextrmat" <-
function(RSobj, mat.no = 1)
{
    obj.name <- deparse(substitute(RSobj))
    if (!(class(RSobj)=="RSmpl" || class(RSobj)=="RSmplext")){
         err.text<-paste(obj.name," is not a sample object - see help(\"rsextrobj\")",sep ="",collapse="")
         stop(err.text)
    }
    if(mat.no > RSobj$n_tot)
         stop("\n\tElement ",mat.no," not available (",obj.name," has ", RSobj$n_tot, " elements).")
    obj<-rsextrobj(RSobj, start = mat.no, end = mat.no)
    RET<-rstats(obj, function(x) matrix(x, nrow = obj$n))[[1]]
    RET
}

