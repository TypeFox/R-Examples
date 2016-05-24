if(is.element("X",ls())) warning("X was replaced",call. = FALSE)
if(is.element("pdata",ls())) warning("pdata was replaced",call. = FALSE)
if(is.element("pos1.cell.counter",ls())) warning("pos1.cell.counter was replaced",call. = FALSE)

load("ACL394data.rda")
X$images$path<-factor(system.file('img', package='RcellData'))