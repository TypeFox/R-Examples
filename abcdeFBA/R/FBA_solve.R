FBA_solve<-function (fba_object, precision = 6, verbosity = FALSE, maximize = TRUE) 
{
	if (maximize == TRUE)
	{
	fba_sol <- Rglpk_solve_LP(obj=fba_object$obj, mat=fba_object$mat,
				dir=fba_object$dir, rhs=fba_object$rhs,
				bounds=fba_object$bounds,types=fba_object$types, 
				max=fba_object$max,control = list("verbose"=verbosity))
	}

	if (maximize == FALSE)
	{
	fba_object$max = FALSE
	fba_sol <- Rglpk_solve_LP(obj=fba_object$obj, mat=fba_object$mat,
				dir=fba_object$dir, rhs=fba_object$rhs,
				bounds=fba_object$bounds, types=fba_object$types, 
				max=fba_object$max,control=list("verbose" =verbosity))
 	 }

	if (fba_sol$status == 1)
	{
	#if (verbosity == TRUE) {
	message("FBA solution Infeasible,Zeroing...")
	#}
	FBA_solution = list()
	FBA_solution$fluxes = rep(0, length(fba_object$obj))
	FBA_solution$objective = 0
	FBA_solution$status = 1
	return(FBA_solution)
	}

	if (fba_sol$status == 0)
	{
	fba_sol$solution <- round(fba_sol$solution, precision)
	fba_sol$optimum <- round(fba_sol$optimum, precision)
	FBA_solution = list()
	FBA_solution$fluxes = fba_sol$solution
	FBA_solution$objective = fba_sol$optimum
	FBA_solution$status = fba_sol$status
	return(FBA_solution)
	}
}
