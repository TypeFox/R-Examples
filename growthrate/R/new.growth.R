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

new.growth <- function(newdata,newtobs,sigma,d,muhatcurve,Khat,tgrid) {
    p = length(newtobs);
    dim(muhatcurve);

    resul = diffquo(newdata,newtobs);

    YInew = resul$YI;
    Xtildanew = resul$Xtilda;
    newtobsgrid = rep(0,p);
    index = rep(0,p);

    ## Find indices of the newtobs in fine grid (tgrid)
    for (i in 1:p) {
        s = newtobs[i];
        tk1 = tgrid[max(which(s>=tgrid))];
        tk2 = tgrid[min(which(s<=tgrid))];

        dtk1 = s-tk1;
        dtk2 = tk2-s;

        if (dtk2 >= dtk1) {
            k = max(which(s>=tgrid));
        } else {
            k = min(which(s<=tgrid));
        }

        newtobsgrid[i] = tgrid[k];
        index[i] = k;
    }

    Khatind = Khat[index,index];
    muhatcurveaver = apply(muhatcurve,2,mean);
    muhatcurveaverind = muhatcurveaver[index];

    sigma0invnew = solve(Khatind);
    mupriornew = muhatcurveaverind;

    r = posteriorobs(sigma0invnew,sigma,mupriornew,Xtildanew,newtobs,YInew);

    muhatmatrix = r$muhatMatrix;
    muhatmatrix = as.matrix(muhatmatrix);
    dim(muhatmatrix);
    Sigmahat = r$Sigmahat;
    rc = posteriordistribcurve(muhatmatrix,Sigmahat,sigma,newtobs,d,YInew);

    Khat = Khatf(tgrid,newtobs,sigma,Sigmahat);

    result = list(muhatcurvenew=rc$muhatcurve,Khatnew=Khat);

    return(result);
}
