### pseudoinverse.R  (2004-09-25)
###
###    Computation of the Pseudoinverse of a Matrix
###
### Copyright 2003-04 Korbinian Strimmer
###
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



pseudoinverse = function (m, tol)
{
    msvd = fast.svd(m, tol)
    
    if (length(msvd$d) == 0)
    {
       return(
            array(0, dim(m)[2:1])
            )
    }
    else
    {
       return( 
            msvd$v %*% (1/msvd$d * t(msvd$u))
            )
     }    
}
