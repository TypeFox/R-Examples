valueEst <-
function (d1, d2, Y, A1, A2){
  ind = (A1==d1)*(A2==d2);
  num = sum (Y*ind);
  denom = sum (ind);
  value = num/denom;

  indPosPos = (A1==1)*(A2==1);
  indPosNeg = (A1==1)*(A2==-1);
  indNegPos = (A1==-1)*(A2==1);
  indNegNeg = (A1==-1)*(A2==-1);

  valPosPos = sum (Y*indPosPos)/sum (indPosPos); 
  valPosNeg = sum (Y*indPosNeg)/sum (indPosNeg);
  valNegPos = sum (Y*indNegPos)/sum (indNegPos);
  valNegNeg = sum (Y*indNegNeg)/sum (indNegNeg);

  list ("value"=value, "valPosPos"=valPosPos, "valPosNeg"=valPosNeg,
        "valNegPos"=valNegPos, "valNegNeg"=valNegNeg);
}
