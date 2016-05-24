### $Id: x_Obj.R 53 2007-04-20 12:05:00Z bhm $
# %=============== x_Obj.m ====================
# %  xObj = x_Obj(Xinput,cova,model,names)
# %    Transforms a matrix of factor setting inputs, Xinput{1,*},
# %         to a model matrix for a regression analysis and performs
# %         initial computations for prediction and testing.
# %    Categorical design variables can be represented by
# %       a numeric vector, a numeric matrix (each unique row is a group),
# %       a character matrix (each row representing a group name),
# %       or a cell array of strings stored as a column vector.
# %    Nonzero elements of cova(1,*) indicate cells of X that are covariates.
# %    Multiple columns (several DF's) of covariate model terms are allowed.
# %    names{1,*} contain the factor names.
# %    Model example: Xinput = {A B C}
# %         model = [ 0 0 0 ; 1 0 0; 0 1 0 ; 0 0 1; 2 0 0; 1 1 0; 1 0 1; 0 1 1; 3 0 0]
# %         ->   Constant + A + B + C + A^2 + A*B + A*C + B*C + A^3
# %
# %   Input:
# %         Xinput{1,*}: factor settings
# %           cova(1,*): covariates
# %          model(*,*): model matrix
# %          names{1,*}: factor names
# %
# %   Output: XObj is a structure with fields
# %         Xinput{1,*}: same as input
# %           cova(1,*): same as input
# %          model(*,*): same as input
# %          names{1,*}: same as input
# %      termNames{1,*}: names for all terms in model
# %            df_error: degrees of freedom for error
# %                nVar: number of input variables
# %            catNames: category variable level names
# %                   X: as X, but with dummy coded categorical variables
# %        X_norm_means: X_norm = normalize(X,X_norm_means,X_norm_stds)
# %         X_norm_stds: --------- // ---------
# %                   D: model matrix for regression analysis: D = my_x2fx(X_norm,model)
# %              D_test: as D, but with Type II adjusted model terms
# %                D_om: as D, but with OM-adjusted model terms
# %             df_D_om: degrees of freedom according to D_om
# %           df_D_test: degrees of freedom according to D_test
# %              Beta_D: regr model: D_om = D*Beta_D  (see linregEst)
# %        VmodelDivS_D: --------------- // ----------------------
# %       VextraDivS1_D: --------------- // ----------------------
# %              Umodel: regr model: Y (unknown) = D_om*Beta  (see linregStart)
# %          VmodelDivS: --------------- // ----------------------
# %         VextraDivS1: --------------- // ----------------------
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function xObj = x_Obj(Xinput,cova,model,names)
# xObj.Xinput = Xinput;
# xObj.cova = cova;
# xObj.model = model;
# xObj.names = names;
# xObj.termNames = makeTermNames(names,model);
# % Make X
# % - where cat variables are coded as contineous
# X=Xinput;
# nVar = size(X,2);
# catNames = cell(size(X));
# for i=1:nVar
#     if(cova(i)==0)
#         [G,GN] = my_grp2idx(X{i});
#         X{i}=dummy_var(G); %X{i}=dummyvar(G);
#         catNames{i}=GN;
#     end
# end
# % D = my_x2fx(X,model);
# % center and scale -> more stable computation
# if(min(sum(model'))>0) % no constant term in model
#     [X_norm X_norm_stds] = absStand(X);
#     X_norm_means = zeros(size(X_norm_stds));
# else
#     [X_norm X_norm_means X_norm_stds] = normalize(X);
# end
# D = my_x2fx(X_norm,model);
# % Make Type~II*-adjusted and OM-adjusted model matrix
# %D_test = orth_D(D,model,'test');
# D_om   = orth_D(D,model,'om');
# D_test = orth_D(D_om,model,'test');  % More stable computation with D_om as input ?
# [D_om_,df_D_om] = c2m(D_om);
# df_D_test = c2df(D_test);
# % Estimate model where X="model matrix" , Y="OM-adjusted model matrix"
# [Beta_D,VmodelDivS_D,VextraDivS1_D] = linregEst(c2m(D),D_om_);
# % Start estimating model where
# %  X = "OM-adjusted model matrix" ,
# %  Y = unknown reponse data
# %%%---%%% [Umodel,Uerror,VmodelDivS,VextraDivS1] = linregStart(D_om_);
# %%%---%%% xObj.df_error = size(Uerror,2);
# [Umodel,VmodelDivS,VextraDivS1] = linregStart(D_om_);
# xObj.df_error = size(Umodel,1) - size(Umodel,2);
# xObj.nVar = nVar;
# xObj.catNames = catNames;
# xObj.X = X;
# xObj.X_norm_means = X_norm_means;
# xObj.X_norm_stds = X_norm_stds;
# xObj.D = D;
# xObj.D_test = D_test;
# xObj.D_om = D_om;
# xObj.df_D_om = df_D_om;
# xObj.df_D_test = df_D_test;
# xObj.Beta_D = Beta_D;
# xObj.VmodelDivS_D = VmodelDivS_D;
# xObj.VextraDivS1_D = VextraDivS1_D;
# xObj.Umodel = Umodel;
# %%%---%%% xObj.Uerror = Uerror;
# xObj.VmodelDivS = VmodelDivS;
# xObj.VextraDivS1 = VextraDivS1;
##############################################################
# Lages (foreloepig) en funksjon som tar D som input (ikke Xinput)
x_Obj = function(D,model){
# % Make Type~II*-adjusted and OM-adjusted model matrix
# %D_test = orth_D(D,model,'test');
D_om   = orth_D(D,model,'om')
D_test = orth_D(D_om,model,'test'); #% More stable computation with D_om as input ?
D_om_  = c2m(D_om)
df_D_om = c2df(D_om)
df_D_test = c2df(D_test)
# % Estimate model where X="model matrix" , Y="OM-adjusted model matrix"
# [Beta_D,VmodelDivS_D,VextraDivS1_D] = linregEst(c2m(D),D_om_);
linregEst_ =linregEst(c2m(D),D_om_)
Beta_D = linregEst_$BetaU
VmodelDivS_D = linregEst_$VmodelDivS
VextraDivS1_D = linregEst_$VextraDivS1
# % Start estimating model where
# %  X = "OM-adjusted model matrix" ,
# %  Y = unknown reponse data
# %%%---%%% [Umodel,Uerror,VmodelDivS,VextraDivS1] = linregStart(D_om_);
# %%%---%%% xObj.df_error = size(Uerror,2);
# [Umodel,VmodelDivS,VextraDivS1] = linregStart(D_om_);
xObj2 = linregStart(D_om_)
# nVar ikke med paa lista
xObj1 = list(df_error = nrow(xObj2$Umodel) - ncol(xObj2$Umodel),
  D=D,
  D_test = D_test,
  D_om = D_om,
  df_D_om = df_D_om,
  df_D_test = df_D_test,
  Beta_D = Beta_D,
  VmodelDivS_D = VmodelDivS_D,
  VextraDivS1_D = VextraDivS1_D,
  termNames = attr(model,"dimnames")[[1]]
  )
c(xObj1,xObj2)
}



# % Dorth = orth_D(D,model,method)
# %    Adjusts/orthogonalizes model matrix D{1,*} according to
# %    method = 'test' (Type II* testing)
# %             'seq'  (Type I testing)
# %             'om'   (adjusts according to model hierarchy)
# %          'ssIII'   (Type III* testinng - dependes on parameterization )
# %   D and model is on the form described in x_Obj.m
# function Dorth = orth_D(D,model,method)
# Dorth = cell(1,size(model,1));
# if(length(D)~=size(model,1))
#     return;
# end
# for i=1:size(model,1)
#     d = D{i};
#     d_adj = d(:,[]);
#     for j=1:size(model,1)
#         switch lower(method)
#             case {'test'}
#                 if(min( model(j,:) - model(i,:)) < 0 )
#                     d_adj = [d_adj D{j}];
#                 end
#             case {'om'}
#                 if( (min( model(j,:) - model(i,:)) < 0)  &  (max( model(j,:) - model(i,:)) <= 0) )
#                     d_adj = [d_adj D{j}];
#                 end
#             case {'seq'}
#                 if(j<i)
#                     d_adj = [d_adj D{j}];
#                 end
#             case {'ssIII'}
#                 if(j~=i)
#                     d_adj = [d_adj D{j}];
#                 end
#         end
#     end
#     Dorth{i} =  adjust(d,d_adj);
# end
####################################################
orth_D = function(D,model,method){
Dorth = vector("list",nrow(model))
if(length(D)!= nrow(model))
    return(Dorth)
# end
for(i in 1:nrow(model)){
   d = D[[i]]
   d_adj = d[,numeric(0)]
   for(j in 1:nrow(model)){
      switch(method, # Bare lagt inn "test" og "om" forloepig
         test = {if(min( model[j,] - model[i,]) < 0 ) d_adj = cbind(d_adj,D[[j]])},
          om  = {if( (min( model[j,] - model[i,]) < 0 )  &  (max( model[j,] - model[i,]) <= 0) ) d_adj = cbind(d_adj,D[[j]])})
         #seq  =,
         #ssIII = )
      }#end
    Dorth[[i]] =  adjust(d,d_adj)
 }# end
Dorth
}# end orth_D
