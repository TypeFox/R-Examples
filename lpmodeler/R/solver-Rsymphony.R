
mipSolve.Rsymphony <- function(p, 
	verbosity=-1,
	time_limit=30,
	gap_limit=-1,
	first_feasible=F,
	write_mps=F,
	write_lp=F,
  ...
)
{
  require(Rsymphony)
  checkDims(p)
  Rsymphony_solve_LP(
    p$obj, p$mat, p$dir, p$rhs, types=p$types, max=p$max,
    verbosity=verbosity,
    time_limit=time_limit,
    gap_limit=gap_limit,
    first_feasible=first_feasible,
    write_mps=write_mps,
    write_lp=write_lp,
    ...)
}

getValue <- function(p,s,name)
  s$solution[get(name,envir=p$vars)]
