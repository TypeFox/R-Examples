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

posteriorobs <- function(Sigma0inv,sigma,muprior,Xtilda,tobs,YI) {
    df = diff(tobs);
    n = length(tobs);

    N = dim(YI)[1];
    n = dim(YI)[2]+1;

    ## Calculate the D, Q matrices and Xtilde in Theorem 2.1
    d = c(1/df,0)+c(0,1/df);
    D = (6/sigma^2)*diag(d);
    d1 = diag(1/df);

    Q = rbind(cbind(rep(0,n-1),d1),rep(0,n))+rbind(rep(0,n),cbind(d1,rep(0,n-1)))+diag(d);
    Q = (3/sigma^2)*Q;

    ## Calculate the estimated posterior mean and covariance
    Sigmahat = solve(Sigma0inv+Q);
    muhatMatrix = c();

    for (k in 1:N) {
        muhat = Sigmahat%*%((Sigma0inv)%*%muprior+D%*%Xtilda[,k]);
        muhatMatrix = cbind(muhatMatrix,muhat);
    }

    result = list(muhatMatrix=muhatMatrix,Sigmahat=Sigmahat);
    return(result);
}

