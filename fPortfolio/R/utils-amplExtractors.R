
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:                     DESCRIPTION:
#  .amplObjval                   Extracts objective function value
#  .amplSolution                 Extracts solution vector
#  .amplModel                    Extracts model file information
#  .amplRun                      Extracts model file information
#  .amplSolver                   Extracts solver name
#  .amplVersion                  Extracts version number
#  .amplPresolve                 Extracts presolve information
###############################################################################


.amplObjval <- 
    function(ampl) 
{  
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts objective function value.
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Return Value:
    ampl$objective
}


# -----------------------------------------------------------------------------


.amplSolution <- 
    function(ampl) 
{
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts solution vector.
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Return Value:
    ampl$solution
}


# -----------------------------------------------------------------------------


.amplModel <- 
    function(ampl) 
{
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts model file information
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    model <- ampl$opt$model
    cat(model, sep="\n")
    
    # Return Value:
    invisible(model) 
}


# -----------------------------------------------------------------------------


.amplRun <- 
    function(ampl) 
{ 
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts run file information.
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Get run file:
    run <- ampl$opt$run
    cat(, sep="\n")
    
    # Return Value:
    invisible(run)
}


# -----------------------------------------------------------------------------


.amplSolver <- 
    function(ampl) 
{   
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts solver information
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Get solver vector:
    solver <- ampl$solver
    cat(solver, sep="\n")
    
    # Return Value:
    invisible(solver)
}


# -----------------------------------------------------------------------------


.amplVersion <- 
    function(ampl) 
{
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Extracts solver version.
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Get Version Number:
    version <- ampl$version
    cat(version, sep="\n")
    
    # Return Value:
    invisible(version)
}


# -----------------------------------------------------------------------------


.amplPresolve <- 
    function(ampl)
{
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #   Extracts presolver results
    
    # Arguments:
    #    ampl - an object as returned by the the rampl[OPT] solver
    #         functions.
    
    # FUNCTION:
    
    # Get presolve information:
    solve <- ampl$opt$solve
    Index <- grep("^Presolve", solve)
    solve <- solve[-(1:(Index-1))]
    Index <- grep("^$", solve)[1]
    solve <- solve[1:(Index-1)]
    cat(solve, sep="\n")
    
    # Return Value:
    invisible(solve)
}


###############################################################################

