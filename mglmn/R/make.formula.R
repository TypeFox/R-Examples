make.formula <-
function(lhs, var.list){
vars.vector <- var.list
if (length(vars.vector) < 2)
return(noquote(paste(lhs, '~',vars.vector)))

rhs <- vars.vector[1]
for (v in vars.vector[2:length(vars.vector)])
rhs <- paste(rhs, '+', v)
noquote(paste(lhs, '~', rhs))

}
