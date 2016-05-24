### $Id: linregEst.R 45 2006-08-15 13:11:29Z bhm $
# %=============== linregEst.m ====================
# %  [BetaU,VmodelDivS,VextraDivS1,msError,errorObs,Yhat] = linregEst(X,Y)
# %        performs multivariate multiple linear regression modelling: Y = XB + E
# %        A principal component regression approach is used (generalised inverse):
# %        X is decomposed using [U,S,V] = svd(X).
# %        A model is made by using r columns of U where r = rank of X.
# %             Umodel is the first r columns of U.
# %        BetaU refers to model with Umodel instead of X
# %        Prediction cannot (not estimable) be made if Xnew*VextraDivS1 is
# %        nonzero.
# %        -  Estimablity can be checked by
# %              is_estimable(Xnew,VextraDivS1)
# %        -  Predictions can be made by
# %              Unew = Xnew * VmodelDivS;
# %              Ypred = Unew*BetaU;
# %              stdYpred = sqrt(sum(Unew .^ 2,2)*msError);
# % Input:
# %       X(*,*) - Regressors
# %       Y(*,*) - Response
# %
# % Output:
# %     BetaU(*,*)       - Parameters in PCR model: Y=Umodel*BetaU + E
# %     VmodelDivS(*,*)  - Umodel = X*VmodelDivS
# %     VextraDivS1(*,*) - For checking estimability
# %     msError(1,*)     - msError for each response
# %     errorObs(*,*)    - error observations (can be used in multivariate testing)
# %     Yhat(*,*)        - fitted values
# %
# % Note:
# %     Only two lines of code:
# %             [Umodel,VmodelDivS,VextraDivS1] = linregStart(X);
# %             [BetaU,msError,errorObs,Yhat] = linregEnd(Umodel,Y);
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linregEst = function(X,Y){
linreg_Start = linregStart(X)
linreg_End =linregEnd(linreg_Start$Umodel,Y)
c(linreg_Start,linreg_End) # Umodel not needed as output
}# end linregEst
