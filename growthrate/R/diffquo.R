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

diffquo <- function(data,tobs) {
    n = length(tobs);
    dtobs = diff(tobs);
    data = as.matrix(data);

    N = dim(data)[1];
    n = dim(data)[2];
    YI = matrix(nrow=N,ncol=(n-1));
    diffdatak = rep(0,(n-1));

    for (k in 1:N) {
        YI[k,] = diff(data[k,])/dtobs;
    }

    Xtilda = matrix(nrow=n,ncol=N);

    w = rep(0,n-2);
    for (i in 2:(n-1)) {
        w[i] = dtobs[i]/(dtobs[i-1]+dtobs[i]);

        for (j in 1:N) {
            Xtilda[i,j] = w[i]*YI[j,i-1]+(1-w[i])*YI[j,i];
        }
    }

    for (j in 1:N) {
        Xtilda[1,j] = YI[j,1];
        Xtilda[n,j] = YI[j,n-1];
    }

    result = list(YI=YI,Xtilda=Xtilda);
    return(result);
}
