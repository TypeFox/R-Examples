.criteria <-
function(loglik,T,prms,n){
  K = prms$K
  p = prms$p
  comp = switch(prms$model,
                'DkBk' = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K^2*(K-1)/2 + K,
                'DkB'  = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K^2*(K-1)/2 + 1,
                'DBk'  = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K*(K-1)/2 + K,
                'DB'   = (K-1) + K*(K-1)+ (K-1)*(p-K/2) + K*(K-1)/2 + 1,
                'AkjBk'= (K-1) + K*(K-1) + (K-1)*(p-K/2) + K^2,
                'AkjB' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K*(K-1)+1,
                'AkBk' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + 2*K,
                'AkB'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K+1,
                'AjBk' = (K-1) + K*(K-1) + (K-1)*(p-K/2) + (K-1)+K,
                'AjB'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + (K-1)+1,
                'ABk'  = (K-1) + K*(K-1) + (K-1)*(p-K/2) + K+1,
                'AB'   = (K-1) + K*(K-1) + (K-1)*(p-K/2) + 2)
  aic = loglik - comp # AIC criterion
  bic = loglik - 1/2 * comp * log(n) # BIC criterion
  T[T<1e-6] = 1e-6
  icl = loglik - 1/2 *  comp * log(n) - sum(T*log(T)) # ICL criterion
  list(aic=aic,bic=bic,icl=icl,nbprm=comp)
}
.estep <-
function(prms,fd,U){
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  K = prms$K
  mu = prms$mean
  prop = prms$prop
  D = prms$D
  d = K-1
  
  QQ = matrix(NA,n,K)
  QQ2 = matrix(NA,n,K)
  T = matrix(NA,n,K)
  
  # Compute posterior probabilities
  for (k in 1:K){
    bk = D[k,p,p]
    mY = prms$my[k,]
    YY = Y-t(matrix(rep(mY,n),p,n)) 
    projYY = YY %*% U %*% t(U) # proj dans l'espace latent
    
    if (d==1){
      for (i in 1:n){QQ[i,k] =  1/D[k,1,1] * sum(projYY[i,]^2) + 1/D[k,p,p]*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)}
    }
    else{
      sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
      for (i in 1:n){QQ[i,k] =  projYY[i,] %*% sY %*% as.matrix(projYY[i, ],p,1) + 1/bk*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)}
    }
  }
  A = -1/2 * QQ
  loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
  for (k in 1:K) {
    # browser()
    T[,k] = 1 / rowSums(exp(0.5*(QQ[,k]*matrix(1,n,K)-QQ)))}
  list(T=T,loglik=loglik)
}
.fstep <-
function(fd,T,lambda){
  G = t(fd$coefs)
  n = nrow(G)
  p = ncol(G)
  d = ncol(T) - 1
  basisobj <- fd$basis
  #W <- (basisobj + t(basisobj))/2
  W <- inprod(basisobj,basisobj)
  Ttilde  = t(apply(T,1,'/',sqrt(colSums(T))))
  eig = svd(ginv(t(G)%*%G%*%W) %*% (t(G)%*%Ttilde%*%t(Ttilde)%*%G%*%W),nu=d,nv=0)
  U   = eig$u[,1:d]

  # Sparse version
  if (lambda>0){
    X = G %*% U
    Utilde = U
    for (i in 1:d){ x.predict = X[,i]
        res.enet = enet(G,x.predict,intercept=FALSE)
        coef     = predict.enet(res.enet,G,type="coefficients",mode="fraction",s=lambda)$coef
        Utilde[,i] = coef/ sqrt(sum(coef^2))
    }
    U = svd(Utilde)$u
  }
  U
}
.FunFEM.main <-
function(fd,K,model='AkjBk',init='kmeans',lambda=0,Tinit=c(),maxit=50,eps=1e-8,graph=F){
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  # 
  # New objects
  Lobs = rep(c(-Inf),1,(maxit+1))
  # 	
  # Initialization of T
  if (init=='user'){ T = Tinit}
  else if (init=='kmeans'){
    T = matrix(0,n,K)
    ind = kmeans(Y,K)$cluster
    for (k in 1:K){ T[which(ind==k),k] = 1}
  }
  else if (init=='random'){
    T = t(rmultinom(n,1,c(rep(1/K,K))))
  }
  else if (init=='hclust'){
    T   = matrix(0,n,K)
    ind = cutree(hclust(dist(Y),method='ward'),K)
    for (k in 1:K){ T[which(ind==k),k] = 1}
  }
  
  V         = .fstep(fd,T,lambda)
  prms      = .mstep(fd,V,T,model=model)
  res.estep = .estep(prms,fd,V)
  T         = res.estep$T
  Lobs[1]   = res.estep$loglik
  
  # Main loop
  Linf_new  = Lobs[1]
  for (i in 1:maxit){
    # The three main steps F, M and E
    #cat('.')
    V         = .fstep(fd,T,lambda)
    prms      = .mstep(fd,V,T,model=model)
    res.estep = .estep(prms,fd,V)
    T         = res.estep$T
    Lobs[i+1] = res.estep$loglik
    
    # Stop criterion
    if (i>=2){
      acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new = Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i])
      if (abs(Linf_new - Linf_old) < eps | is.na(Linf_new)) {break}
    }
  }
  
  # Graphic option
  if (graph){
    par(mfrow=c(1,2))
    plot(as.data.frame(as.matrix(Y) %*% V[,1:2]),col=max.col(T),xlab='axis 1',ylab='axis 2',pch=20)
    plot(Lobs[1:i],xlab='iterations',ylab='vraisemblance Espace observe',col=2,pch=20)
  }
  
  # Returning the results
  cls  = max.col(T)
  crit = .criteria(Lobs[(i+1)],T,prms,n);
  res  = list(model=model,K=K,cls=cls,P=T,prms=prms,U=V,aic=crit$aic,bic=crit$bic,icl=crit$icl,
              loglik=Lobs[2:(i+1)],ll=Lobs[i+1],nbprm=crit$nbprm)
}
.mstep <-
function(fd,U,T,model){
  # 12 different submodels: [DkBk] ... [AkjBk]
  # Initialization
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  K = ncol(T)
  d = K-1
  
  mu   = matrix(NA,K,K-1)
  m   = matrix(NA,K,p)
  prop = rep(c(NA),1,K)
  D = array(0,c(K,p,p))
  
  # Projection
  X = Y %*% U
  
  # Estimation
  for (k in 1:K){
    
    nk  = sum(T[,k])
    # Prior Probability
    prop[k] = nk / n
    # Mean in the latent space
    mu[k,]  = colSums((T[,k]*matrix(1,n,d))* X) / nk
    # Observed space
    m[k,]  = colSums(T[,k]*matrix(1,n,p)* Y) / nk
    YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
    Ck  = crossprod(T[,k]*matrix(1,n,p)* YY, YY) / (nk-1) #crossprod(x,y) = t(x) %*% y
    C   = cov(Y)
    
    # Estimation of Delta k amongst 8 submodels
    if (model=='DkBk'){
      D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DkB'){
      D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DBk'){
      D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='DB'){
      D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkjBk'){
      if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkjB'){
      if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkBk'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AkB'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    
    if (model=='AjBk'){
      if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else {
        D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='AjB'){
      if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else{
        D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='ABk'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
      bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
    }
    if (model=='AB'){
      if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
        D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
      bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
      bk[bk<=0] = 1e-3
      D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))	
    }
  }
  prms = list(K=K,p=p,mean=mu,my=m,prop=prop,D=D,model=model)
}
.Random.seed <-
c(403L, 10L, 244503110L, -1932176735L, -39627469L, 1781443746L, 
-247643244L, 1912396311L, 1316757005L, 1637301028L, -67459318L, 
440527477L, 1962383823L, -1407622442L, -1921595776L, 1151813827L, 
738207329L, -1873412400L, -324163538L, 1957404121L, 2099954155L, 
-1124219814L, -1700250084L, -258896097L, -167124107L, 414501388L, 
1034472418L, -1456932675L, 1412982583L, 1618145246L, 1719328024L, 
-1544621957L, -2041321959L, -2018488664L, 1530297430L, -820258607L, 
1942160515L, 748267282L, -1019559260L, 2083425895L, -621805635L, 
662204724L, -173182310L, -1493372699L, 1130738303L, 243301382L, 
1541934992L, -636980685L, 1382754769L, 200342400L, -2043038114L, 
-587797079L, -682935205L, -1307696022L, 2008833292L, -663470257L, 
524696581L, 936108732L, 2088645522L, 1381724301L, 1311065799L, 
659853230L, 1936538120L, 1338867211L, -306235671L, 864014776L, 
-2084082074L, 1473374977L, 1198308755L, -1664330174L, -576730636L, 
57488375L, 1537099757L, 1992695940L, -1150601558L, -686691307L, 
-29238097L, 475468470L, -1669738272L, -1538106333L, -648191167L, 
-837194192L, 49892430L, -1701264199L, -1414086453L, 666742138L, 
-1614745476L, -283978561L, 1958781333L, -516609620L, -2062353534L, 
-1828863651L, -841121705L, 782194494L, 676798328L, -333204389L, 
1050297273L, -1337523704L, -1295931210L, 853446641L, 2030811683L, 
1570419314L, 1077225668L, 550682247L, 2081125981L, -814704556L, 
-414356870L, 141709765L, 1227031327L, 209301158L, -689389520L, 
-402005293L, -998565647L, -1435884512L, -1836064066L, 1367197385L, 
934806651L, 948976202L, -1651507540L, -829408849L, -1271459483L, 
-1329217764L, 950641522L, 1507500013L, 532271143L, -1315551154L, 
1369468072L, 1616675499L, 1060571081L, -1651180584L, 1280952966L, 
1121610593L, 1156804339L, 1647227234L, 1889054676L, -782271785L, 
-401932083L, 1706326116L, -880519606L, -109517899L, -1114538225L, 
-1675166698L, 915293888L, -440622973L, -496395103L, 1091089040L, 
1498820974L, 527871769L, -212311125L, 63000986L, -531798692L, 
617716831L, 1109271605L, 1766809548L, -1320183006L, -663345795L, 
20653687L, 1810349982L, 1062149720L, 925185979L, -1860937383L, 
1075499752L, 300469014L, 2060377233L, -241962685L, -2663342L, 
864624356L, 624052519L, 218802557L, 1208279156L, 1384219226L, 
1793986213L, -761168705L, 972701510L, 117784016L, -321011725L, 
396470033L, -625645760L, -1194389730L, 1901076457L, -2054146789L, 
-256094422L, -1002749876L, -54474737L, 466174789L, 368076028L, 
-816876718L, -799071027L, -1626307449L, -1455816850L, 121196872L, 
-44804789L, 1688680617L, -1464071688L, -1906124506L, -1161967039L, 
1128539091L, -2039877246L, 1628951988L, 1985446199L, 710960173L, 
-900026044L, 2107278954L, -1845079851L, -592298385L, -898794634L, 
479792288L, 1738204515L, 1122395649L, -1967463056L, 1334228238L, 
190363513L, -1975908597L, -1444124614L, -1759905988L, -1481952129L, 
838117333L, -1298478612L, 282880834L, 758166685L, -301093609L, 
22717566L, 280091704L, 1353038619L, -1609830279L, 951105736L, 
991858934L, -1456737999L, -368748189L, 111409420L, 1812041276L, 
-1897413390L, -399075088L, -155216060L, -711552936L, 670065210L, 
-1667577616L, -464297284L, -50673132L, 85232834L, 2061981280L, 
-587455108L, -530078512L, 1861555234L, 1299507368L, -1244275316L, 
-1800257092L, 2132723762L, 401707840L, -271573324L, 572232008L, 
-1090334022L, -1057939696L, 1665147116L, 262256724L, 810469138L, 
-532725344L, 1005840332L, 592617312L, -492090174L, 2014318120L, 
1659546028L, 884761404L, 1194525554L, -1457549520L, -137879516L, 
-77996488L, -777899302L, -263030896L, -31025924L, 2144660884L, 
1942930242L, -926021120L, 829933628L, 1170367056L, 2125193378L, 
-118004920L, -1652917748L, 1459788508L, -982437806L, 2082966656L, 
1011664724L, 1621532648L, 1638658810L, -399251632L, 1096156204L, 
-1332442764L, 139873938L, -1216167328L, 1141977228L, -744173216L, 
-1331152670L, -1354715992L, -1825196340L, -1896121860L, 2120166322L, 
1187709936L, -1169096252L, -225100456L, 461042042L, -456363600L, 
-1546098244L, 211033684L, -1000392254L, -997268704L, -900574276L, 
75837520L, -1971741726L, 986693480L, 1702246220L, 834072316L, 
-756488462L, -532964096L, -1433789324L, 1360596488L, 541894074L, 
1824073168L, -1383061524L, 2040490900L, -776140590L, -547949472L, 
1816964428L, 289755936L, 552487554L, 970052072L, -807494676L, 
-785936836L, 328238706L, -1205703824L, 951561508L, -1628726728L, 
127287706L, 1727655376L, -1828633156L, -28332396L, 1533702786L, 
-1002377792L, -722282692L, 673006224L, -1667848414L, 2073150984L, 
-1953996020L, -1505782884L, 848583698L, -588141696L, -109253868L, 
333591848L, -100220870L, 1596486608L, 1636661868L, -1703654348L, 
854782098L, 504994144L, 2113166668L, -1903381024L, 32283554L, 
497351592L, 1213475980L, -2113182788L, 2101449842L, -1096144912L, 
1172485188L, 1776920664L, -87535558L, -20469648L, -1442443716L, 
740591380L, -1096389950L, 808096480L, 1691993212L, 415706320L, 
-1534164062L, -89522648L, 956420364L, -1151617092L, -93798094L, 
-1408952128L, -1541901772L, -1671091512L, 339022522L, -1733747568L, 
-327829780L, 1220982356L, 117831314L, -1289558240L, 1885247436L, 
-1867110176L, -351428926L, -2072980440L, 1350712236L, -498050372L, 
-325770766L, 1096873520L, 1353226148L, -483798600L, 246000474L, 
-473439088L, 1795185532L, -842961132L, -1907184318L, -1038840064L, 
667614140L, 1815375696L, 477933858L, 1407421000L, -52193268L, 
-634762148L, -1687735854L, 1919088128L, -1226329772L, 25899240L, 
655192314L, -1618610480L, -1751013076L, 999082228L, -1521017454L, 
386001248L, 1174857740L, 300886752L, 1432203874L, 2026687144L, 
-1467970996L, 391534460L, 308714802L, 1677985392L, 452139844L, 
-856369192L, 1541943162L, 840744496L, 1458321340L, 1306830420L, 
-1455843262L, 219451552L, 1643638460L, -873453104L, 2122171234L, 
115221224L, 1857491148L, -1380019716L, -1406043662L, 2068089472L, 
-756066828L, 288294536L, 1237086010L, 1727951696L, 1454552172L, 
-1821111916L, -864234414L, -921073184L, 494607564L, -1850610528L, 
-834529918L, -697442072L, 1195286252L, 610262332L, -1298686862L, 
-1187913138L, 472801948L, 636217357L, -1920355209L, 1951900800L, 
-384142626L, -1205804629L, -1248830803L, 1142505242L, -1831198760L, 
-319902655L, 1762237171L, 1724177684L, -1972066142L, -1029085769L, 
1317992833L, 1603977094L, -1648595516L, -1023030027L, 854943071L, 
760154040L, 402500246L, 1884050595L, -1657180571L, -2098046078L, 
977213680L, -1215004519L, -51603317L, -768851588L, 2105720746L, 
616398815L, 1139993961L, -1499901570L, 315785580L, -727516739L, 
-427021689L, 815558256L, -1129728754L, 881269307L, -382803683L, 
1977127882L, 285560040L, -297013295L, -1167294013L, -1317360956L, 
582920562L, -1480444857L, -591800303L, 161104278L, 713889012L, 
2016239237L, 2036545487L, -711428408L, 1831235430L, -1247373933L, 
-110401803L, -62617902L, -550252256L, -1852899063L, -1494833605L, 
-1211204724L, -784704998L, 1478781583L, -421546535L, 673061486L, 
-168375812L, -1934464147L, 43693335L, 838435040L, -1940359490L, 
-1051446133L, -1977498867L, 603480890L, 2113493880L, 2071472417L, 
-1693711917L, 387182836L, 1017106818L, 1067537431L, 39502561L, 
-1730085850L, 726673700L, 19638741L, -1586526849L, 1570225176L, 
1483993910L, -1887246909L, 76844037L, 480069538L, 1364165520L, 
2074101561L, -628001813L, -9403172L, 1618091274L, -1081261057L, 
797095497L, 309855198L, -1855988852L, -1728136355L, -949548185L, 
-1321664368L, -1230758226L, -1956627621L, 831523645L, 226296234L, 
-477073848L, -1950982543L, -2095366173L, 1840016100L, 2021921298L, 
566232551L, 63672753L, 470120438L, -1472181740L, -1737886427L, 
-1148279249L, 40535016L, 1748128966L, 1070330227L, 1283168725L, 
-543396430L, -1636629632L, -1017071767L, 131128283L, -2051800276L, 
-993338054L, -1368710161L, 556742265L, -1302007538L, -1122914084L, 
-1663170739L, -552565065L, -2130656192L, 139115806L, 791836651L, 
-130376851L, 1417271770L, 2117178264L, 1188448769L, 1288948659L, 
-346436908L, 1667447778L, 1396715127L, -1686063167L, -793549498L, 
-2137164668L, 439114293L, -1072397793L, 885364088L, -793234346L, 
-352311069L, -780247387L, 1686942530L, 1391474992L, -867065511L, 
-1609076021L, -1725506756L, -742263958L, 179298719L, -660301783L, 
-531144002L, -684624980L, 1668161789L, 865963079L, -418759888L, 
-19256242L, -1986215813L, -1804307491L, -1564149238L, -897055177L
)
