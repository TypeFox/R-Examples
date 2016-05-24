initialization_corr <-
function(dataset, rGvalue, fullparam, zeroset, manualinitial) {

  DataMZ = dataset[dataset$rG==rGvalue,]
  M = cbind(DataMZ$M1, DataMZ$M2); 
  nMZ = nrow(M); 
  muMMZ = mean(M);
  P = cbind(DataMZ$P1, DataMZ$P2);
  muPMZ = mean(P);
  Sx = cov(cbind(M,P));

  DataDZ = dataset[!(dataset$rG==rGvalue),]
  M = cbind(DataDZ$M1, DataDZ$M2); 
  nDZ = nrow(M);
  muMDZ = mean(M);
  P = cbind(DataDZ$P1, DataDZ$P2);
  muPDZ = mean(P);
  Sy = cov(cbind(M,P));

  if (nrow(DataMZ) == 0 || nrow(DataDZ) == 0) {
    n = nrow(dataset);
    Dataset = dataset[order(dataset$rG,decreasing=TRUE),];

    DataMZ = Dataset[1:round(n/2),];
    M = cbind(DataMZ$M1, DataMZ$M2); 
    nMZ = nrow(M); 
    muMMZ = mean(M);
    P = cbind(DataMZ$P1, DataMZ$P2);
    muPMZ = mean(P);
    Sx = cov(cbind(M,P));

    DataDZ = Dataset[(round(n/2)+1):n,]
    M = cbind(DataDZ$M1, DataDZ$M2); 
    nDZ = nrow(M);
    muMDZ = mean(M);
    P = cbind(DataDZ$P1, DataDZ$P2);
    muPDZ = mean(P);
    Sy = cov(cbind(M,P));
  }

  muM = (muMMZ*nMZ+muMDZ*nDZ)/(nMZ+nDZ);
  muP = (muPMZ*nMZ+muPDZ*nDZ)/(nMZ+nDZ);

  ACE = ((Sx[1,1]+Sx[2,2])*nMZ+(Sy[1,1]+Sy[2,2])*nDZ)/2/(nMZ+nDZ);
  AC = Sx[1,2];
  eM = sqrt(max(0,ACE-AC));
  aM = sqrt(max(0,2*(Sx[1,2]-Sy[1,2])));
  cM = sqrt(max(0,2*Sy[1,2]-Sx[1,2])); 

  ACEM = ((Sx[1,3]+Sx[2,4])*nMZ + (Sy[1,3]+Sy[2,4])*nDZ)/2/(nMZ+nDZ);
  ACM = (Sx[1,4]+Sx[2,3])/2;
  eC = ifelse(eM>0,(ACEM-ACM)/eM,0);
  aC = ifelse(aM>0,(Sx[1,4]+Sx[2,3]-Sy[1,4]-Sy[2,3])/aM,0);
  cC = ifelse(cM>0,(Sy[1,4]+Sy[2,3]-Sx[1,4]/2-Sx[2,3]/2)/cM,0);

  ACEACEU = ((Sx[3,3]+Sx[4,4])*nMZ + (Sy[3,3]+Sy[4,4])*nDZ)/2/(nMZ+nDZ);
  ACEU = max(0, ACEACEU-aC^2-cC^2-eC^2);
  ACU = Sx[3,4]-aC^2-cC^2;
  eU = sqrt(max(0, ACEU-ACU));
  halfACU = Sy[3,4]-aC^2/2-cC^2;
  aU = sqrt(max(0, 2*(ACU-halfACU)));
  cU = sqrt(max(0, 2*halfACU-ACU));

  aP = sqrt(aC^2 + aU^2);
  rA = ifelse(abs(aP)>1e-8, aC/aP, 0);  rA = min(rA,0.999); rA = max(rA,-0.999);
  cP = sqrt(cC^2 + cU^2);
  rC = ifelse(abs(cP)>1e-8, cC/cP, 0);  rC = min(rA,0.999); rC = max(rA,-0.999);
  eP = sqrt(eC^2 + eU^2);
  rE = ifelse(abs(eP)>1e-8, eC/eP, 0);  rE = min(rA,0.999); rE = max(rA,-0.999);

  initialvalues = list(muM=muM, aM=aM, cM=cM, eM=eM, muP=muP, aP=aP, cP=cP, eP=eP, rA=rA, rC=rC, rE=rE);
  usedparam = setdiff(fullparam, zeroset);

  set1 = c("muM","aM", "cM", "eM","muP","aP","cP","eP");
  usedset1 = set1[is.element(set1, usedparam)];
  d1 = length(usedset1);
  param1 = matrix(0, 3, d1);
  if(d1 > 0) {
    for (i in 1:d1) {
      ind = which(usedset1[i] == names(initialvalues));
      param1[1,i] = initialvalues[[ind]];
    }
    param1[2,] = param1[1,] - 10;
    param1[3,] = param1[1,] + 10;
    colnames(param1) = usedset1;
  }

  set2 = c("rA", "rC", "rE");
  usedset2 = set2[is.element(set2, usedparam)];
  d2 = length(usedset2);
  param2 = matrix(0, 3, d2)
  if(d2 > 0) {
    for (i in 1:d2) {
      ind = which(usedset2[i] == names(initialvalues));
      param2[1,i] = initialvalues[[ind]];
    }
    param2[2,] = -0.999;
    param2[3,] = 0.999;
    colnames(param2) = usedset2;
  }

  usedset3 = setdiff(setdiff(usedparam, usedset1),usedset2);
  d3 = length(usedset3);
  param3 = matrix(rep(c(0, -10, 10), d3), 3, d3);
  colnames(param3) = usedset3;
  param = cbind(param1, param2, param3);

  paramnames = colnames(param);
  initial=param[1,];
  if(length(manualinitial) > 0) {
    manualnames = names(manualinitial);
    for (i in 1:length(manualnames)) {
      localind = which(paramnames %in% manualnames[i]);
      if(length(localind)<1) {
        cat(paste(manualnames[i],'can not be manually initialized.\n'));
      } else {
        initial[localind] = manualinitial[[i]];
        param[2,localind] = initial[localind] - 10;
        param[3,localind] = initial[localind] + 10;
        if (is.element(manualnames[i], set2)) {
          param[2,localind] = -0.999;
          param[3,localind] = 0.999;
        }
      }
    }
  }

  return(list(paramnames=paramnames,initial=initial,lower = param[2,],upper=param[3,]));
}
