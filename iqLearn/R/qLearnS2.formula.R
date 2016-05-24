qLearnS2.formula <-
function (formula, data=list(), treatName,
                              intNames, ...){
  mf <- model.frame (formula=formula, data=data);
  Y <- model.response (mf);
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
  H2 = cbind (x[,-c (1,intOrder)]);
  A2 = x[,txtCol];
  mainNames = colnames (H2);
  
  if (length (intNames) > 0){
    ints = which (mainNames %in% intNames);
  }
  else{
    ints = NULL;
  }
      
  est <- qLearnS2.default (H2, Y, A2, ints, ...);
  est$call <- match.call ();
  est$formula <- formula;
  est
}
