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

growth <- function(data,tobs,sigma,d) {
    data = as.matrix(data);
    obsdata = data;
    dim(obsdata);

    ## obsdata is the data matrix of dimension N (number of subjects) times n (number of time points).
    N = dim(obsdata)[1];
    n = dim(obsdata)[2];

    ## calculate the difference of consecutive time points
    dobspts = diff(tobs);

    ## calculate the vector y in the paper
    YI = matrix(nrow=N,ncol=(n-1));
    dobsdatak = rep(0,(n-1));

    for (k in 1:N) {
        YI[k,] = diff(obsdata[k,])/dobspts;
    }

    ## Estimate the prior mean based on nonparametric approach and obtain Xtilda (ybar in the paper)
    resulmuprior = priormeanobs(YI,tobs);
    muprior = resulmuprior$muprior;
    Xtilda = resulmuprior$Xtilda;

    ## Estimate the prior precision based on CLIME
    re.clime = clime(t(Xtilda),perturb=TRUE);
    re.cv = cv.clime(re.clime,loss=c("likelihood", "tracel2"));
    re.clime.opt <- clime(t(Xtilda), re.cv$lambdaopt, perturb =TRUE);
    Sigma0inv <- re.clime.opt$Omegalist[[1]];

    ## Get the posterior mean and variance-covariance matrix at observed time points
    resulpost = posteriorobs(Sigma0inv,sigma,muprior,Xtilda,tobs,YI);
    muhatMatrix = resulpost$muhatMatrix;
    Sigmahat = resulpost$Sigmahat;

    ## The posterior mean curves (muhatcurve) and variance covariance operator (Khat) evaluated at the fine grid
    r = posteriordistribcurve(muhatMatrix,Sigmahat,sigma,tobs,d,YI);
    tgrid = r$tgrid;
    muhatcurve = r$muhatcurvematrix;
    Khat = Khatf(tgrid,tobs,sigma,Sigmahat);

    ## Return the muhat curves for each subject, variance-covariance operator and the fine grid.
    result = list(muhatcurve=muhatcurve,Khat=Khat,tgrid=tgrid);
    return(result);
}

