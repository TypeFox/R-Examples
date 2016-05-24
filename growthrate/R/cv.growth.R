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

cv.growth <- function(data,tobs,d,K,a,b,r) {
    obsdata = data;

    ## The rows of the matrix obsdata correspond to the individuals or subjects
    N = dim(obsdata)[1];
    n = dim(obsdata)[2];

    ## Take differences of the observed time points
    dobspts = diff(tobs);

    ## First we construct the yij using all time points
    YI = matrix(nrow = N,ncol = (n-1));

    for (k in 1:N) {
        YI[k,] = diff(obsdata[k,])/dobspts;
    }

    ## Create the vector Xtilde (or ybold) based on YI
    Xtilda = matrix(nrow = n,ncol = N);

    w = rep(0,n-2);
    for (i in 2:(n-1)) {
        w[i] = dobspts[i]/(dobspts[i-1]+dobspts[i]);

        for (j in 1:N) {
            Xtilda[i,j] = w[i]*YI[j,i-1]+(1-w[i])*YI[j,i];
        }
    }

    for (j in 1:N){
        Xtilda[1,j] = YI[j,1];
        Xtilda[n,j] = YI[j,n-1];
    }

    ybar = t(Xtilda);

    ## Estimate the prior mean based on nonparametric approach
    muprior = apply(ybar,2,mean);
    muprior = t(muprior);

    ## Estimate the prior variance using CLIME
    re.clime = clime(ybar,perturb = TRUE);
    re.cv = cv.clime(re.clime,loss = c("likelihood", "tracel2"));
    re.clime.opt <- clime(ybar, re.cv$lambdaopt, perturb = TRUE);
    Sigma0inv <- re.clime.opt$Omegalist[[1]];

    ## Take out the rth observed time point
    tm = tobs[-r];
    dobsptsm = diff(tm);
    YIm = matrix(nrow = N,ncol = length(dobsptsm));
    obsdatam = obsdata[,-r];

    ## Construct YIm with one less time point
    for (k in 1:N) {
        YIm[k,] = diff(obsdatam[k,])/dobsptsm;
    }

    n = length(dobsptsm)+1;

    ## Create the vector Xtilde (or ybold in the paper) based on new YIm
    Xtilda = matrix(nrow = n,ncol = N);

    w = rep(0,n-2);
    for (i in 2:(n-1)) {
        w[i] = dobsptsm[i]/(dobsptsm[i-1]+dobsptsm[i]);
        for (j in 1:N) {
            Xtilda[i,j] = w[i]*YIm[j,i-1]+(1-w[i])*YIm[j,i];
        }
    }

    for (j in 1:N) {
        Xtilda[1,j] = YIm[j,1];
        Xtilda[n,j] = YIm[j,n-1];
    }

    ybar = t(Xtilda);

    ## Based on the estimated priors with the full data set select the prior mean and variance-covariance taking out the rth time point
    mupriorr = muprior[-r];
    Sigma0invr = Sigma0inv[-r,-r];

    ## Define the values of sigma we are considering and initialize the variables
    sigmavec = seq(a,b,length = K);
    CVe = rep(0,K);
    e = rep(0,N);

    for (k in 1:K) {
        sigma = sigmavec[k];

        ## Calculate the posterior mean and variance covariance at observation time points for each value of sigma
        resulpost = posteriorobs(Sigma0invr,sigma,mupriorr,Xtilda,tm,YIm);

        muhatMatrix = resulpost$muhatMatrix;
        Sigmahat = resulpost$Sigmahat;

        ## Calculate the posterior distribution at the fine grid for each subject
        result = posteriordistribcurve(muhatMatrix,Sigmahat,sigma,tm,d,YIm);

        muhatcurve = result$muhatcurvematrix;
        tgrid = result$tgrid;
        t = tobs;
        Khat = Khatf(tgrid,tm,sigma,Sigmahat);

        ## Reconstruct at time point t(r) the value of y for each subject
        tgrids = tgrid[tgrid >= t[r-1] & tgrid <= t[r]];
        itgrids = which(tgrid >= t[r-1] & tgrid <= t[r]);
        ngrids = length(tgrids);

        muhatcurve1 = muhatcurve[,itgrids];
        N = dim(muhatcurve1)[1];
        ng = dim(muhatcurve1)[2];

        Khat1 = Khat[itgrids,itgrids];

        ## Calculate the integral of muhatcurve1
        dtgrids = diff(tgrids);
        integ = rep(0,N);

        for (j in 1:N) {
            integ1 = sum(dtgrids*((muhatcurve1[j,1:(ng-1)]+muhatcurve1[j,2:ng])/2));
            integ[j] = integ1/(t[r]-t[r-1]);

            ## Calculate the error in the reconstruction for each subject
            e[j] = (YI[j,r-1]-integ[j])^2;
        }

        int1 = dtgrids%*%Khat1[-1,];
        int2 = int1[-1]%*%dtgrids;

        CVe[k] = mean(e) + (1/((t[r]-t[r-1])^2))*int2;
    }

    ## Calculate and return the cross-validation error removing observed time point r for the different values of sigma in sigmavec and the vector of sigma values in sigmavec
    CVer = CVe;

    result = list(CVer = CVer,sigmavec = sigmavec);
    return(result);
}
