learnIQ2.formula <-
function (formula, data=list(), treatName,
  intNames, ...){
  mf <- model.frame (formula=formula, data=data);
  Y <- model.response (mf);
  x <- model.matrix (attr (mf, "terms"), data=mf);
  xNames = colnames (x);
  findTxt = xNames %in% treatName;
  txtCol = which (findTxt*1==1);
  findFirst = xNames %in% paste (treatName, xNames, sep=":");
  txtFirst = which (findFirst*1==1);
  findLast = xNames %in% paste (xNames, treatName, sep=":");
  txtLast = which (findLast*1==1);
  intOrder = c (txtCol, txtFirst, txtLast);
  numInt = length (intOrder);
  H2 = cbind (x[,-c (1,intOrder)]);
  A2 = x[,txtCol];
  mainNames = colnames (H2);
  ints = which (mainNames %in% intNames);
  est <- learnIQ2.default (H2, Y, A2, ints, ...);
  est$call <- match.call ();
  est$formula <- formula;
  est
}
