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

Khatf <- function(tgrid,tobs,sigma,Sigmahat) {
    ng = length(tgrid);
    nobs = length(tobs);

    ## Construct Ktilde
    A = 0;
    for (k in 1:(nobs-1)) {
        t1 = tobs[k];
        t2 = tobs[k+1];

        tgrids = tgrid[tgrid>=t1 & tgrid<t2];
        p1 = length(tgrids);
        Sigma = matrix(nrow=p1,ncol=p1);
        for (i in 1:p1) {
            for (j in 1:p1) {
                Sigma[i,j] = (min(tgrids[i],tgrids[j])-t1)-
                    (tgrids[j]-t1)*(tgrids[i]-t1)/(t2-t1)-
                        (3/(t2-t1)^3)*(tgrids[i]-t1)*(t2-tgrids[i])*(tgrids[j]-t1)*(t2-tgrids[j]);
            }
        }
        A = bdiag(A,Sigma);
    }

    A = as.matrix(A);
    Ktilde = A[-1,-1];
    Ktilde = bdiag(Ktilde,0);
    Ktilde = as.matrix(Ktilde);

    ## Now we construct Kstar
    nobs = length(tobs);
    tg = tgrid;
    ng = length(tg);
    Kstar = matrix(nrow=ng,ncol=ng);

    for (i in 1:(ng-1)) {
        for (j in 1:(ng-1)) {
            s = tg[i];
            t = tg[j];

            ## Find interval where s is in
            tk1 = tobs[max(which(s>=tobs))];
            tk2 = tobs[min(which(s<tobs))];
            k = max(which(s>=tobs)); # tells you which interval s is in

            ## Find interval where t is in
            tl1 = tobs[max(which(t>=tobs))];
            tl2 = tobs[min(which(t<tobs))];
            l = max(which(t>=tobs)); # tells you which interval t is in

            aks = 1-(s-tk1)/(tk2-tk1)-3*(s-tk1)*(tk2-s)/((tk2-tk1)^2);
            bks = (s-tk1)/(tk2-tk1)-3*(s-tk1)*(tk2-s)/((tk2-tk1)^2);
            alt = 1-(t-tl1)/(tl2-tl1)-3*(t-tl1)*(tl2-t)/((tl2-tl1)^2);
            blt = (t-tl1)/(tl2-tl1)-3*(t-tl1)*(tl2-t)/((tl2-tl1)^2);

            Kstar[i,j] = aks*Sigmahat[k,l]*alt+bks*Sigmahat[k+1,l]*alt+aks*Sigmahat[k,l+1]*blt+bks*Sigmahat[k+1,l+1]*blt;
        }
    }

    Kstar[ng,ng] = Sigmahat[nobs,nobs];
    for (i in 1:(ng-1)) {

        ## Find the interval where s is in
        s = tg[i];
        tk1 = tobs[max(which(s>=tobs))];
        tk2 = tobs[min(which(s<tobs))];
        k = max(which(s>=tobs)); # tells you which interval s is in

        aks = 1-(s-tk1)/(tk2-tk1)-3*(s-tk1)*(tk2-s)/((tk2-tk1)^2);
        bks = (s-tk1)/(tk2-tk1)-3*(s-tk1)*(tk2-s)/((tk2-tk1)^2);
        Kstar[i,ng] = aks*Sigmahat[k,nobs]+bks*Sigmahat[k+1,nobs];
        Kstar[ng,i] = Kstar[i,ng];
    }

    ## Construct Khat
    Khat = sigma^2*Ktilde+Kstar;
    return(Khat);
}

