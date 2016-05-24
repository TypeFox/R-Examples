################################################################################
#    This file is part of growthrate.                                          #
#                                                                              #
#    growthrate is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by      #
#    the Free Software Foundation, either version 3 of the License, or         #
#    (at your option) any later version.                                       #
#                                                                              #
#    growthrate is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with growthrate.  If not, see <http://www.gnu.org/licenses/>.       #
################################################################################

priormeanobs <- function(YI,tobs) {
    df = diff(tobs);
    n = length(tobs);

    N = dim(YI)[1];
    n = dim(YI)[2]+1;
    Xtilda = matrix(nrow=n,ncol=N);

    w = rep(0,n-2);
    for (i in 2:(n-1)) {
        w[i] = df[i]/(df[i-1]+df[i]);

        for (j in 1:N) {
            Xtilda[i,j] = w[i]*YI[j,i-1]+(1-w[i])*YI[j,i];
        }
    }

    for (j in 1:N) {
        Xtilda[1,j] = YI[j,1];
        Xtilda[n,j] = YI[j,n-1];
    }

    muprior = apply(Xtilda,1,mean);

    result <- list(muprior=muprior,Xtilda=Xtilda);
    return(result);
}

