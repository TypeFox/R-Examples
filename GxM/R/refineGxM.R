refineGxM <-
function(param, modelname) {

  paramnames = names(param);
  ind_aM = which(paramnames %in% 'aM');    sign_aM = sign(param[ind_aM]);    param[ind_aM] = abs(param[ind_aM]);
  ind_cM = which(paramnames %in% 'cM');    sign_cM = sign(param[ind_cM]);    param[ind_cM] = abs(param[ind_cM]);
  ind_eM = which(paramnames %in% 'eM');    sign_eM = sign(param[ind_eM]);    param[ind_eM] = abs(param[ind_eM]);

  ## model 3
  if(modelname == 'Chol') {
    ind_aC = which(paramnames %in% 'aC');    param[ind_aC] = sign_aM * param[ind_aC];
    ind_cC = which(paramnames %in% 'cC');    param[ind_cC] = sign_cM * param[ind_cC];
    ind_eC = which(paramnames %in% 'eC');    param[ind_eC] = sign_eM * param[ind_eC];
    ind_aU = which(paramnames %in% 'aU');    param[ind_aU] = abs(param[ind_aU]);
    ind_cU = which(paramnames %in% 'cU');    param[ind_cU] = abs(param[ind_cU]);
    ind_eU = which(paramnames %in% 'eU');    param[ind_eU] = abs(param[ind_eU]);
  }

  ## model 4
  if(modelname == 'CholGxM') {
    ind_aC = which(paramnames %in% 'aC');                param[ind_aC] = sign_aM * param[ind_aC];
    ind_alphaC = which(paramnames %in% 'alphaC');        param[ind_alphaC] = sign_aM * param[ind_alphaC];

    ind_cC = which(paramnames %in% 'cC');                param[ind_cC] = sign_cM * param[ind_cC];
    ind_kappaC = which(paramnames %in% 'kappaC');        param[ind_kappaC] = sign_cM * param[ind_kappaC];

    ind_eC = which(paramnames %in% 'eC');                param[ind_eC] = sign_eM * param[ind_eC]; 
    ind_epsilonC = which(paramnames %in% 'epsilonC');    param[ind_epsilonC] = sign_eM * param[ind_epsilonC];

    ind_aU = which(paramnames %in% 'aU');                sign_aU = sign(param[ind_aU]);    param[ind_aU] = abs(param[ind_aU]);
    ind_alphaU = which(paramnames %in% 'alphaU');        param[ind_alphaU] = sign_aU * param[ind_alphaU];

    ind_cU = which(paramnames %in% 'cU');                sign_cU = sign(param[ind_cU]);    param[ind_cU] = abs(param[ind_cU]);
    ind_kappaU = which(paramnames %in% 'kappaU');        param[ind_kappaU] = sign_cU * param[ind_kappaU];

    ind_eU = which(paramnames %in% 'eU');                sign_eU = sign(param[ind_eU]);    param[ind_eU] = abs(param[ind_eU]);
    ind_epsilonU = which(paramnames %in% 'epsilonU');    param[ind_epsilonU] = sign_eU * param[ind_epsilonU];
  }

  ## model 5
  if(modelname == 'NLMainGxM') {
    ind_aU = which(paramnames %in% 'aU');                sign_aU = sign(param[ind_aU]);    param[ind_aU] = abs(param[ind_aU]);
    ind_alphaU = which(paramnames %in% 'alphaU');        param[ind_alphaU] = sign_aU * param[ind_alphaU];

    ind_cU = which(paramnames %in% 'cU');                sign_cU = sign(param[ind_cU]);    param[ind_cU] = abs(param[ind_cU]);
    ind_kappaU = which(paramnames %in% 'kappaU');        param[ind_kappaU] = sign_cU * param[ind_kappaU];

    ind_eU = which(paramnames %in% 'eU');                sign_eU = sign(param[ind_eU]);    param[ind_eU] = abs(param[ind_eU]);
    ind_epsilonU = which(paramnames %in% 'epsilonU');    param[ind_epsilonU] = sign_eU * param[ind_epsilonU];
  }

  ## model 6
  if(modelname == 'CorrGxM') {
    ind_aP = which(paramnames %in% 'aP');                sign_aP = sign(param[ind_aP]);    param[ind_aP] = abs(param[ind_aP]);
    ind_alphaP = which(paramnames %in% 'alphaP');        param[ind_alphaP] = sign_aP * param[ind_alphaP];

    ind_cP = which(paramnames %in% 'cP');                sign_cP = sign(param[ind_cP]);    param[ind_cP] = abs(param[ind_cP]);
    ind_kappaP = which(paramnames %in% 'kappaP');        param[ind_kappaP] = sign_cP * param[ind_kappaP];

    ind_eP = which(paramnames %in% 'eP');                sign_eP = sign(param[ind_eP]);    param[ind_eP] = abs(param[ind_eP]);
    ind_epsilonP = which(paramnames %in% 'epsilonP');    param[ind_epsilonP] = sign_eP * param[ind_epsilonP];

    ind_rA = which(paramnames %in% 'rA');                param[ind_rA] = sign_aM * sign_aP * param[ind_rA];
    ind_rC = which(paramnames %in% 'rC');                param[ind_rC] = sign_cM * sign_cP * param[ind_rC];
    ind_rE = which(paramnames %in% 'rE');                param[ind_rE] = sign_eM * sign_eP * param[ind_rE];
  }

  ## model 7
  if(modelname == 'CholNonLin') {
    ind_aC = which(paramnames %in% 'aC');                param[ind_aC] = sign_aM * param[ind_aC];
    ind_cC = which(paramnames %in% 'cC');                param[ind_cC] = sign_cM * param[ind_cC];
    ind_eC = which(paramnames %in% 'eC');                param[ind_eC] = sign_eM * param[ind_eC]; 

    ind_aU = which(paramnames %in% 'aU');                sign_aU = sign(param[ind_aU]);    param[ind_aU] = abs(param[ind_aU]);
    ind_delta1 = which(paramnames %in% 'delta1');        param[ind_delta1] = sign_aM * sign_aU * param[ind_delta1];

    ind_cU = which(paramnames %in% 'cU');                sign_cU = sign(param[ind_cU]);    param[ind_cU] = abs(param[ind_cU]);
    ind_delta2 = which(paramnames %in% 'delta2');        param[ind_delta2] = sign_cM * sign_cU * param[ind_delta2];

    ind_eU = which(paramnames %in% 'eU');                sign_eU = sign(param[ind_eU]);    param[ind_eU] = abs(param[ind_eU]);
    ind_delta3 = which(paramnames %in% 'delta3');        param[ind_delta3] = sign_eM * sign_eU * param[ind_delta3];
  }

  ## model 8
  if(modelname == 'CorrNonLin') {
    ind_aP = which(paramnames %in% 'aP');                sign_aP = sign(param[ind_aP]);    param[ind_aP] = abs(param[ind_aP]);
    ind_lambda1 = which(paramnames %in% 'lambda1');      param[ind_lambda1] = sign_aM * sign_aP * param[ind_lambda1];

    ind_cP = which(paramnames %in% 'cP');                sign_cP = sign(param[ind_cP]);    param[ind_cP] = abs(param[ind_cP]);
    ind_lambda2 = which(paramnames %in% 'lambda2');      param[ind_lambda2] = sign_cM * sign_cP * param[ind_lambda2];

    ind_eP = which(paramnames %in% 'eP');                sign_eP = sign(param[ind_eP]);    param[ind_eP] = abs(param[ind_eP]);
    ind_lambda3 = which(paramnames %in% 'lambda3');      param[ind_lambda3] = sign_eM * sign_eP * param[ind_lambda3];

    ind_rA = which(paramnames %in% 'rA');                param[ind_rA] = sign_aM * sign_aP * param[ind_rA];
    ind_rC = which(paramnames %in% 'rC');                param[ind_rC] = sign_cM * sign_cP * param[ind_rC];
    ind_rE = which(paramnames %in% 'rE');                param[ind_rE] = sign_eM * sign_eP * param[ind_rE];
  }

  return(param);
}
