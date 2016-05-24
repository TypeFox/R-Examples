# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 28-05-2014
# last modification: 10-10-2014
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.



# random.subsampling(Y)
# prediction.formula(data,H,uloadings,vloadings,CH,means.X,means.Y,sigma.X,sigma.Y)
# pre.screening(X,coeff)

# random subsamplings
grey_area=function()
{
    #source(paste(directory.working,"prediction.R",sep=""))
    #source(paste(directory.working,"CI.subsampling.prediction.R",sep=""))
    #X.learn=X
    #Y.learn=Y
    #result=prediction(X.learn,Y.learn,indice.X,H=H,X.test=X.learn)
    
    #quant.MSC=quantile(result$Y.hat.learn[Y.learn=="MSC",1,H],0.01)
    #quant.NonMSC=quantile(result$Y.hat.learn[Y.learn=="Non-MSC",1,H],0.99)
    #cat("Quantile 1% of the density of the MSC values on the learning set are lower than",quant.MSC ," \n")
    #cat("Quantile 1% of the density of the non-MSC values on the learning set are higher than",quant.NonMSC ," \n")
    #a=density(result$Y.hat.learn[Y.learn=="Non-MSC",1,H])
    #names(a)
    #b=density(result$Y.hat.learn[Y.learn=="MSC",1,H])
    #ind.a=which(a$x>quant.MSC & a$x<quant.NonMSC)
    #ind.b=which(b$x>quant.MSC & b$x<quant.NonMSC)

}

