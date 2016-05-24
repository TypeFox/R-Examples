### cor0.test.R  (2004-01-15)
###
###    Tests of Vanishing Correlation
###    
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
###
### This file is part of the `GeneNet' library for R and related languages.
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


# Test of vanishing partial correlation
#  H0: rho == 0
#  HA: rho != 0
#
# Input:  observed partial correlation r 
#         degree of freedom kappa
# Output: p-value
cor0.test = function(r, kappa, method=c("student", "dcor0", "ztransform"))
{
  method = match.arg(method)

  if (method == "student") # exact method
  {
       # t is distributed with df=kappa-1 
       t = r*sqrt((kappa-1)/(1-r*r))
       
       # two-sided test around zero
       pval = 2*pt(-abs(t), df=kappa-1)
  }
 
  if (method == "dcor0") # exact method
  {
       # two-sided test around zero
       pval = 2*pcor0(-abs(r), kappa)
  }
  
  if (method == "ztransform") # approximate method
  {
    # apply Fisher's z-transform
    z = z.transform(r)
    
    # then use two-sided normal test around zero
    sd = 1/sqrt(kappa-2)
    pval = 2*pnorm(-abs(z), mean=0, sd=sd)
  }

  return(pval)
}
