### $Id: manova_5050.R 53 2007-04-20 12:05:00Z bhm $
# %=============== manova_5050.m ====================
# % function results = manova_5050(xObj,Y,stand)
# %    Takes a model object (created by x_Obj.m) togheter with a
# %     a matrix, Y, of responses and produces 50-50 MANOVA output.
# %     Use stand=1 to standardize responses (stand=0 otherwise).
# %
# % --OUTPUT-- results is a structure with fields:
# %     termNames: name of model terms (including "error").
# %       exVarSS: (Sum of SS for each response)/(Sum of total SS for each response).
# %            df: degrees of freedom - adjusted for other terms in model.
# %         df_om: degrees of freedom - adjusted for terms contained in actual term.
# %           nPC: number of principal components used for testing.
# %           nBU: number of principal components used as buffer components.
# %       exVarPC: variance explained by nPC components
# %       exVarBU: variance explained by (nPC+nBU) components
# %       pValues: 50-50 MANOVA p-values.
# %    outputText: 50-50 MANOVA results as text.
# %
# %  NOTE: The function can be called by
# %        manova_5050(x_Obj(Xinput,cova,model,xNames),Y,stand)
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function results = manova_5050(xObj,Y,stand)
# partBufDim = 0.5;        %%% !! hard coded constant !! %%%
# minBufDim = 0;           %%% !! hard coded constant !! %%%
# maxBufDim = 100000000;   %%% !! hard coded constant !! %%%
# minErrDf  = 3;           %%% !! hard coded constant !! %%%
# cp = -1;                 %%% !! hard coded constant !! %%%
# part1 = 0.9;             %%% !! hard coded constant !! %%%
# part2 = 0.5;             %%% !! hard coded constant !! %%%
# part  = [part1,part2]';  %%% !! hard coded constant !! %%%
# yNames=[];
# if(stand)
#    Y = stdStand(Y);
# end
# model = xObj.model;
# xyObj = xy_Obj(xObj,Y,yNames);
# nTerms = length(xyObj.xObj.df_D_test);
# %results.Yhat = xyObj.Yhat;
# results.termNames = xyObj.xObj.termNames;
# results.exVarSS = xyObj.ss / xyObj.ssTot;
# results.df = [xyObj.xObj.df_D_test xyObj.xObj.df_error];
# results.df_om = [xyObj.xObj.df_D_om xyObj.xObj.df_error];
# nPC = [];
# nBU = [];
# exVarPC = [];
# exVarBU = [];
# pValues = [];
# normY = norm(Y);
# %errorData = xyObj.errorObs;
# for i=1:nTerms
#    modelData = xyObj.hypObs{i};
#    %if(normY < 1e-250 | norm([errorData',modelData'])/normY < 1e-12)% Singularity problems
#    %   [exVar1_,exVar2_,dim_,dimX_,dimY_,bufferDim_,D_,E_,A_,M_,pD_,pE_,pA_,pM_] = ...
#    %      ffmanovatest(modelData(:,[]),errorData(:,[]),part,partBufDim,minBufDim, ...
#    %      maxBufDim,minErrDf,cp,stand);
#    %%%---%%%
#    if(iscell(xyObj.errorObs))
#        normTest =  norm([(xyObj.errorObs{1})',modelData']);
#        dfError = xyObj.errorObs{2};
#    else
#        normTest =  norm([xyObj.errorObs',modelData']);
#        dfError = size(xyObj.errorObs,1);
#    end
#    if(normY < 1e-250 | normTest/normY < 1e-12)% Singularity problems
#       [exVar1_,exVar2_,dim_,dimX_,dimY_,bufferDim_,D_,E_,A_,M_,pD_,pE_,pA_,pM_] = ...
#          ffmanovatest(modelData(:,[]),zeros(dfError,0),part,partBufDim,minBufDim, ...
#          maxBufDim,minErrDf,cp,stand);
#    else
#       [exVar1_,exVar2_,dim_,dimX_,dimY_,bufferDim_,D_,E_,A_,M_,pD_,pE_,pA_,pM_] = ...
#          ffmanovatest(modelData,xyObj.errorObs,part,partBufDim,minBufDim,...
#          maxBufDim,minErrDf,cp,stand);
#    end
#    nPC = [nPC dimY_];
#    nBU = [nBU bufferDim_];
#    exVarPC = [exVarPC exVar1_];
#    exVarBU = [exVarBU exVar2_];
#    pValues = [pValues pA_];
# end
# results.nPC = nPC;
# results.nBU = nBU;
# results.exVarPC = exVarPC;
# results.exVarBU = exVarBU;
# results.pValues = pValues;
# % Start making outputText
#     outputText=[];
#     outputText=outLine(outputText,sprintf('  --- 50-50 MANOVA Version 2.0 --- %d objects -- %d responses:',size(Y,1),size(Y,2)));
#     approx = 0;
#     %names = strvcat(strvcat(results.termNames),'Source'); % Changed - Octave
#     names = '';
#     for i=1:length(results.termNames)
#         names = strvcat(names,results.termNames{i});
#     end
#     names = strvcat(names,'Source');
#     outputText=outLine(outputText,sprintf('  %s  DF        exVarSS nPC nBu exVarPC exVarBU    p-Value ',names(size(model,1)+2,:)));
#     for i=1:(nTerms+1)
#        s1 = sprintf('  %s',names(i,:));
#        s2 = sprintf('%4d',results.df(i));
#        if(results.df(i)==results.df_om(i))
#           dfFull = '       ';
#       else
#           dfFull = sprintf('(%d)',results.df_om(i));
#       end
#        dfFull = strjust(sprintf('%7s',dfFull),'left');
#        s3 = sprintf('%s',dfFull(1:5));
#        s4 = sprintf('  %8.6f',results.exVarSS(i));
#        if(i <= nTerms)
#           s5 = sprintf(' %3d',nPC(i));
#           s6 = sprintf(' %3d ',nBU(i));
#           s7 = sprintf(' %5.3f  ',exVarPC(i));
#           s8 = sprintf(' %5.3f    ',exVarBU(i));
#           if(pValues(i)<2)
#              s9 = sprintf('%8.6f ',pValues(i));
#              if(nPC(i)>2 & results.df(i)>2)
#                 s10 =sprintf('x');
#                 approx = 1;
#              else
#                 s10 = ' ';
#              end
#           else
#              s9 =sprintf(' ....... ');
#              s10 = ' ';
#           end
#           outputText=outLine(outputText,sprintf('%s%s%s%s%s%s%s%s%s%s',s1,s2,s3,s4,s5,s6,s7,s8,s9,s10));
#        end
#     end
#     if(stand)
#        s5 = sprintf(' - STANDARDIZATION ON  ');
#     else
#        s5 = sprintf(' - STANDARDIZATION OFF ');
#     end
#     if(approx)
#        s6 = sprintf('- x Approx p');
#     else
#        s6 = sprintf('------------');
#     end
#     outputText=outLine(outputText,sprintf('%s%s%s%s%s%s%s%s%s%s',s1,s2,s3,s4,s5,s6));
# results.outputText = outputText;
#############################################################################
manova5050 = function(xyObj,stand){
#if(stand) Y = stdStand(Y)
#xyObj = xy_Obj(xObj,Y)
model = xyObj$xObj$model
nTerms = length(xyObj$xObj$df_D_test)
nPC = c()
nBU = c()
exVarPC = c()
exVarBU = c()
pValues = c()
normY = norm(xyObj$Y)
for(i in 1:nTerms){
   modelData = xyObj$hypObs[[i]]
   if(is.list(xyObj$errorObs)){
       normTest =  norm(rbind(xyObj$errorObs[[1]],modelData))
       dfError = xyObj$errorObs[[2]]
   }else{
       normTest =  norm(rbind(xyObj$errorObs,modelData))
       dfError = nrow(xyObj$errorObs)
   }#end
   if(normY < 1e-250 | normTest/normY < 1e-12){ #% Singularity problems
      res = ffmanovatest(modelData[,numeric(0)],matrix(0,dfError,0),stand)
   }else{
      res = ffmanovatest(modelData,xyObj$errorObs,stand)
   }#end
   nPC = c(nPC ,res$dimY)
   nBU = c(nBU ,res$bufferDim)
   exVarPC = c(exVarPC ,res$exVar1)
   exVarBU = c(exVarBU ,res$exVar2)
   pValues = c(pValues ,res$pA)
}#end
list(termNames=xyObj$xObj$termNames,
     exVarSS = xyObj$ss/xyObj$ssTot,
     df = c(xyObj$xObj$df_D_test, xyObj$xObj$df_error),
     df_om = c(xyObj$xObj$df_D_om, xyObj$xObj$df_error),
     nPC = nPC,
     nBU = nBU,
     exVarPC = exVarPC,
     exVarBU = exVarBU,
     pValues = pValues,
     stand = stand)
}

# % Subfunction
# function outputText=outLine(text,line)
# outputText=strvcat(text,line);
