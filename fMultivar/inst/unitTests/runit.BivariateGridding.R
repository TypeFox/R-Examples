
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA


################################################################################
# FUNCTION:           DESCRIPTION:
#  gridData            Generates a grid data set
#  persp.gridData      Generates a perspective plot from a grid data object
#  contour.gridData    Generates a contour plot from a grid data object
# REQUIREMENT:        DESCRIPTION:
#  akima               Package
#  spatial             Package
################################################################################


test.gridData = 
    function()
{
    #  gridData            Generates a grid data set
    #  persp.gridData      Generates a perspective plot from a grid data object
    #  contour.gridData    Generates a contour plot from a grid data object
     
    # Generate Grid Data:
    gD <- gridData()
    
    # Perspective Plot:
    persp(gD)
    
    # Contour Plot:
    contour(gD) 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.gridDataPlot = 
    function()
{   
    #  gridData            Generates a grid data set
    #  persp.gridData      Generates a perspective plot from a grid data object
    #  contour.gridData    Generates a contour plot from a grid data object
      
    # Generate Akima interpolated Grid Data:
    require(akima)
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = runif(999)-0.5
    y = runif(999)-0.5
    z = cos(2*pi*(x^2+y^2))
    ans = akimaInterp(x, y, z, extrap = FALSE)
    persp(ans)
    title(main = "Akima Interpolation") 
    contour(ans)
    title(main = "Akima Interpolation") 
    
    # Generate Kriged Grid Data:
    require(spatial)
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = runif(999)-0.5
    y = runif(999)-0.5
    z = cos(2*pi*(x^2+y^2))
    ans = krigeInterp(x, y, z, extrap = FALSE)
    persp(ans)
    title(main = "Kriging") 
    contour(ans)
    title(main = "Kriging") 
    
    # Return Value:
    return()    
}


################################################################################

