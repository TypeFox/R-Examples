learnIQ1var.formula <-
function (formula, data=list(), treatName,
  intNames, method="homo", cmObject, ...){
  if (method=="homo"){
    stop ("method must be hetero to use formula specification");
  }
  mf <- model.frame (formula=formula, data=data);
  x <- model.matrix (attr (mf, "terms"), data=mf);
  xNames = colnames (x);
  findTxt = xNames %in% treatName;
  txtCol = which (findTxt*1==1);
  findFirst = xNames %in% paste (treatName, xNames, sep=":");
  findLast = xNames %in% paste (xNames, treatName, sep=":");

  if (length (findFirst) > 0){
    txtFirst = which (findFirst*1==1);
  }
  else{
    txtFirst = NULL;
  }
  if (length (findLast) > 0){
    txtLast = which (findLast*1==1);
  }
  else{
    txtLast = NULL;
  }
  intOrder = c (txtCol, txtFirst, txtLast);
  numInt = length (intOrder);
  H1CVar = cbind (x[,-c (1,intOrder)]);
  A1 = x[,txtCol];
  mainNames = colnames (H1CVar);

  if (length (intNames) > 0){
    ints = which (mainNames %in% intNames);
  }
  else{
    ints = NULL;
  }
  
  est <- learnIQ1var.default (cmObject, H1CVar, ints, method, ...);
  est$call <- match.call ();
  est$formula <- formula;
  est
}
