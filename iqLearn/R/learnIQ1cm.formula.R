learnIQ1cm.formula <-
function (formula, data=list(), treatName,
                                intNames, s2object, ...){
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
  H1CMean = cbind (x[,-c (1,intOrder)]);
  A1 = x[,txtCol];
  mainNames = colnames (H1CMean);

  if (length (intNames) > 0){
    ints = which (mainNames %in% intNames);
  }
  else{
    ints = NULL;
  }
  est <- learnIQ1cm.default (s2object, H1CMean, A1, ints, ...);
  est$call <- match.call ();
  est$formula <- formula;
  est
}
