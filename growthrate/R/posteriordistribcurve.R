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

posteriordistribcurve <- function(muhatmatrix,Sigmahat,sigma,tobs,d,YI) {
    muhatmatrix = t(muhatmatrix);
    n = length(tobs);
    t0 = min(tobs);
    T = max(tobs);
    N = dim(muhatmatrix)[1];
    tg = seq(t0,T,length=d);
    ngrid = length(tg);

    muhatcurvematrix = matrix(nrow=N,ncol=ngrid);
    varhatcurvematrix = matrix(nrow=N,ncol=ngrid);

    dt = diff(tobs);

    for (k in 1:N) {
        muhat = muhatmatrix[k,];
        muhatcurve = muhat[1];
        varhatcurve = Sigmahat[1,1];
        ti <- tg[1];

        for (i in 1:(n-1)) {
            td = tg[tobs[i]<tg & tg<=tobs[i+1]];
            muhatcurvei = muhat[i]+(muhat[i+1]-muhat[i])*
                ((td-tobs[i])/dt[i])+6*(td-tobs[i])*(tobs[i+1]-td)*
                    (YI[k,i]-(muhat[i]+muhat[i+1])/2)/dt[i]^2;

            ai = 1-(td-tobs[i])/dt[i]-3*(td-tobs[i])*(tobs[i+1]-td)/dt[i]^2;
            bi = (td-tobs[i])/dt[i]-3*(td-tobs[i])*(tobs[i+1]-td)/dt[i]^2;
            Ktildett = (td-tobs[i])-(td-tobs[i])^2/dt[i]-3*(td-tobs[i])^2*(tobs[i+1]-td)^2/dt[i]^3;

            varhatcurvei = (sigma^2)*Ktildett+(ai^2)*Sigmahat[i,i]+2*ai*bi*Sigmahat[i,i+1]+bi^2*Sigmahat[i+1,i+1];

            muhatcurve = cbind(muhatcurve,t(muhatcurvei));
            varhatcurve = cbind(varhatcurve,t(varhatcurvei));

            ti <- cbind(ti,t(td));
        }

        muhatcurvematrix[k,] = muhatcurve;
        varhatcurvematrix[k,] = varhatcurve;
        tgrid = ti;
    }

    result = list(muhatcurvematrix=muhatcurvematrix,varhatcurvematrix=varhatcurvematrix,tgrid=tgrid);
    return(result);
}

