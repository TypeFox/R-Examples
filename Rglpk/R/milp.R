## MILP constructor
## borrowed from KH's relations package

## objective: numeric vector
## constraints: list [[1]] constr_mat
##                   [[2]] constr_dir,
##                   [[3]] constr_rhs),
## integers:
## maximum: logical

MILP <-
function(objective, constraints, bounds = NULL, types = NULL,
         maximum = FALSE)
{
    ## In the simples case, 'constraints' is a (not necessarily named)
    ## list with mat, dir and rhs.  Advanced solvers might allow for
    ## more advanced constraints, but let's worry about this later (and
    ## maybe also a little MILP_constraints() wrapper ...).

    structure(list(objective = slam::as.simple_triplet_matrix(objective), constraints = constraints,
                   bounds = bounds, types = types, maximum = maximum),
              class = "MILP")
}
