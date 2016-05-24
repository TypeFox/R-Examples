.bic <-
function(loglik,prms,n,T){
	K = prms$k;p = prms$p;d = prms$d;model = prms$model;
	tau_bar=1/K*sum(d*(p-(d+1)/2))
	if (model=='AkjBkQkDk'){comp = K*p +  K -1 + K *tau_bar + K*sum(d) + K}
	if (model=='AkjBQkDk'){comp = K*p +  K -1 + K *tau_bar + K*sum(d) + 1}
	if (model=='AkBkQkDk'){comp = K*p +  K -1 + K *tau_bar + K + K }
	if (model=='AkBQkDk'){comp = K*p +  K -1 + K *tau_bar + K + 1}
	if (model=='ABkQkDk'){comp = K*p +  K -1 + K *tau_bar + 1 + K}
	if (model=='ABQkDk'){comp = K*p + K -1 + K *tau_bar + 2}
# 	aic = 2 * loglik - 2* comp
# 	Z = matrix(0,n,K)
# 	for (i in 1:n) Z[i,which.max(T[i,])] = 1
# 	icl = 2 * loglik - comp * log(n) - sum(log(T[Z==1]))
  aic = loglik - comp
	bic = loglik - comp/2 * log(n)
  T[T<1e-6] = 1e-6
  icl = loglik - comp/2 * log(n) - sum(T*log(T)) # ICL criterion
  list(aic=aic,bic=bic,icl=icl)
}
.diago <-
function(v){
	if (length(v)==1){ res = v }
	else { res = diag(v)}
	res
}
.estep <-
function(prms,fd){
	## Initialization
	X = t(fd$coefs)
	model = prms$model; k = prms$k; p = prms$p;
	a = prms$a; b = prms$b; d = prms$d; prop = prms$prop;
	m = prms$m; Q = prms$Q
	A = matrix(NA,nrow(X),k)
	T = matrix(NA,nrow(X),k)
	## Cost function computing
	for(i in 1:k){
		# Projection of test data in the eigenspace Ei
		Qi = as.matrix(Q[i,,1:d[i]])
		Pa = (as.matrix(X - .repmat(m[i,],nrow(X),1)) %*% Qi) %*% t(Qi)
		Pb = Pa + as.matrix(.repmat(m[i,],nrow(X),1) - X)
		
		#Compute cost function A_i(x)
		if (model=='AkjBkQkDk' | model=='AkjBQkDk'){
			ai = a[1:d[i],i];
			A[,i] = t(diag(Pa %*% Qi %*% .diago(1/ai) %*% t(Qi) %*% t(Pa))
				+ (1/b[i] * rowSums(Pb^2)) + sum(log(ai))
				+ (p-d[i]) * log(b[i]) - 2 * log(prop[i]) + p * log(2*pi));
		}
		
		if (model=='AkBkQkDk' | model=='AkBQkDk' | model=='ABkQkDk' | model=='ABQkDk'){
			A[,i] = t(1/a[i] * rowSums(Pa^2) + (1/b[i] * rowSums(Pb^2)) + d[i] * log(a[i])
				+ (p-d[i]) * log(b[i]) - 2 * log(prop[i]) + p * log(2*pi));
		}
	}

	## Posterior probabilities
	for (i in 1:k) {T[,i] = 1 / rowSums(exp(0.5*(t(.repmat(A[,i],k,1))-A)))}
	T
}
.loglikelihood <-
function(prms,fd){
	X = t(fd$coefs)
	## Initialization
	model = prms$model; k = prms$k; p = prms$p;
	a = prms$a; b = prms$b; d = prms$d; prop = prms$prop;
	m = prms$m; Q = prms$Q
	A = matrix(NA,nrow(X),k)
	
	## Cost function computing
	for(i in 1:k){
		# Projection of test data in the eigenspace Ei
		Qi = as.matrix(Q[i,,1:d[i]])
		Pa = (as.matrix(X - .repmat(m[i,],nrow(X),1)) %*% Qi) %*% t(Qi)
		Pb = Pa + as.matrix(.repmat(m[i,],nrow(X),1) - X)
		
		#Compute cost function A_i(x)
		if (model=='AkjBkQkDk' | model=='AkjBQkDk'){
			ai = a[1:d[i],i];
			A[,i] = t(diag(Pa %*% Qi %*% .diago(1/ai) %*% t(Qi) %*% t(Pa))
				+ (1/b[i] * rowSums(Pb^2)) + sum(log(ai))
				+ (p-d[i]) * log(b[i]) - 2 * log(prop[i]) + p * log(2*pi));
		}
		
		if (model=='AkBkQkDk' | model=='AkBQkDk' | model=='ABkQkDk' | model=='ABQkDk'){
			A[,i] = t(1/a[i] * rowSums(Pa^2) + (1/b[i] * rowSums(Pb^2)) + d[i] * log(a[i])
				+ (p-d[i]) * log(b[i]) - 2 * log(prop[i]) + p * log(2*pi));
		}
	}
	A = -1/2 * A
 	l = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
}
.mstep <-
function(fd,T,thd,model){
	X = t(fd$coefs)
	cls = max.col(T)
	## Initialization
	k = ncol(T); N <- nrow(X); p <- ncol(X)
	n = matrix(0,1,k); d = matrix(0,1,k); prop = matrix(0,1,k)
	m = matrix(0,k,p); Tr = matrix(0,1,k); L = matrix(0,k,p)
	V = matrix(0,p,(k*p)); Q = array(NA,c(k,p,p))
	fpcaobj = list()

	## Compute intrinsic dim and eigen-decomposition
	for (i in 1:k){
# 		fdi = fd
		# Class parameters
		n[i] = sum(T[,i])
		prop[i] = n[i] / N
		m[i,] = colSums(t(.repmat(T[,i],p,1)) * X) / n[i]
# 		XX = as.matrix(X - .repmat(m[i,],N,1))
# 		Si = t(T[,i] * XX) %*% XX / n[i]
		# Eigen-decomposition
		dc = .mypca.fd(fd,T[,i],nharm=p)
		L[i,] = dc$val
		L[i,L[i,]<0] = 0
		# barplot(dc$val); Sys.sleep(1)
		Q[i,,] = dc$U
		Tr[i] = sum(L[i,]);
		# Find intrinsic dimensions using the sree-test of Cattell
		dc$val[dc$val<0] = 0
		sc = abs(diff(dc$val))
		d[i] = p-1
		for (j in 1:(p-2)){
			if (prod(sc[(j+1):length(sc)] < thd * max(sc))){
				d[i] = j
				break
			}
		}
		fpcaobj[[i]] = dc
	}
	
	## Computing model parameters
	if (model=='AkjBkQkDk'){a = matrix(0,p,k)} else {a = matrix(0,1,k)}
	b = matrix(0,1,k)

	if (model=='AkjBkQkDk' | model=='AkjBQkDk'){
		for (i in 1:k){
			a[1:d[i],i] = L[i,1:d[i]]
			b[i] = (Tr[i] - sum(L[i,1:d[i]])) / (p-d[i])
		}
		if (model=='AkjBQkDk') b[1:k] = mean(b)
	}

	if (model=='AkBkQkDk' | model=='AkBQkDk' | model=='ABkQkDk' | model=='ABQkDk'){
		for (i in 1:k){
			a[i] = sum(L[i,1:d[i]]) / d[i]
			b[i] = (Tr[i] - sum(L[i,1:d[i]])) / (p-d[i])
		}
		if (model=='AkBQkDk' | model=='ABQkDk') b[1:k] = mean(b)
		if (model=='ABkQkDk' | model=='ABQkDk') a[1:k] = mean(a)
	}

	## Returning model paramters
	prms <- list(model=model,k=k,p=p,a=a,b=b,m=m,prop=prop,d=d,Q=Q,fpcaobj=fpcaobj)
}
.mypca.fd <-
function(fdobj, Ti, nharm = 2, harmfdPar=fdPar(fdobj), centerfns = TRUE){
#  Carry out a functional PCA with regularization
#  Arguments:
#  FDOBJ      ... Functional data object
#  NHARM     ... Number of principal components or harmonics to be kept
#  HARMFDPAR ... Functional parameter object for the harmonics
#  CENTERFNS ... If TRUE, the mean function is first subtracted from each function
#
#  Returns:  An object PCAFD of class "pca.fd" with these named entries:
#  harmonics  ... A functional data object for the harmonics or eigenfunctions
#  values     ... The complete set of eigenvalues
#  scores     ... A matrix of scores on the principal components or harmonics
#  varprop    ... A vector giving the proportion of variance explained
#                 by each eigenfunction
#  meanfd     ... A functional data object giving the mean function

  #  Check FDOBJ
  if (!(inherits(fdobj, "fd"))) stop(
    "Argument FD  not a functional data object.")

  #  compute mean function and center if required
# browser()
meanfd <- mean.fd(fdobj)
# if (centerfns) fdobj <- center.fd(fdobj)

if (centerfns){
	coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj$coefs))) * fdobj$coefs, 1, sum) / sum(Ti)
	fdobj$coefs <- sweep(fdobj$coefs, 1, coefmean)
	meanfd$coefs = as.matrix(data.frame(mean=coefmean))
}

  #  get coefficient matrix and its dimensions
  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]
  coefnames <- dimnames(coef)
  if (nrep < 2) stop("PCA not possible without replications.")

  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type

  #  set up HARMBASIS
  #  currently this is required to be BASISOBJ
  harmbasis <- basisobj

  #  set up LFDOBJ and LAMBDA
  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda

  #  compute CTEMP whose cross product is needed
  ctemp <- coef
  
  #  set up cross product and penalty matrices
#   Cmat <- crossprod(t(ctemp))/nrep
Cmat = (Ti * ctemp) %*% t(ctemp) / nrep

  Jmat <- eval.penalty(basisobj, 0)
  if(lambda > 0) {
    Kmat <- eval.penalty(basisobj, Lfdobj)
    Wmat <- Jmat + lambda * Kmat
  } else {    Wmat <- Jmat  }
  Wmat <- (Wmat + t(Wmat))/2

  #  compute the Choleski factor of Wmat
  Lmat    <- chol(Wmat)
  Lmat.inv <- solve(Lmat)

  #  set up matrix for eigenanalysis
  if(lambda > 0) { Cmat <- t(Lmat.inv) %*% Jmat %*% Cmat %*% Jmat %*% Lmat.inv  }
  else { Cmat <- Lmat %*% Cmat %*% t(Lmat)  }
  
  #  eigenalysis
  Cmat    <- (Cmat + t(Cmat))/2
  result  <- eigen(Cmat)
  eigvalc <- result$values
  eigvecc <- as.matrix(result$vectors[, 1:nharm])
  sumvecc <- apply(eigvecc, 2, sum)
  eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]

  varprop <- eigvalc[1:nharm]/sum(eigvalc)


  harmcoef <- Lmat.inv %*% eigvecc
  U = t(Lmat) %*% eigvecc
  harmscr  <- t(ctemp) %*% U
  
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
    harmnames <- list(coefnames[[1]], harmnames,"values")
    harmfd   <- fd(harmcoef, basisobj, harmnames)

  pcafd  <- list(harmonics=harmfd,values=eigvalc,scores=harmscr,U=U,varprop=varprop,meanfd=meanfd)
  class(pcafd) <- "pca.fd"
  return(pcafd)
}
.Random.seed <-
c(403L, 10L, 168271670L, 10116960L, 276637091L, -1050734911L, 
-1366514000L, -1987247666L, -1897335495L, -856102581L, -1654072582L, 
908942332L, 2095651903L, -914324203L, 458944044L, -1611926782L, 
644272093L, 957729495L, 426132926L, 505184248L, -2002043941L, 
4109881L, -133947512L, -1475879114L, 1138498929L, 1187903395L, 
1507245810L, -1081121212L, 956025095L, 1206096093L, 937833940L, 
-915707654L, -1440346043L, 1433090207L, -286498266L, 961459888L, 
1960074067L, 2079831665L, -43385952L, -1435702466L, -2027457975L, 
2066884347L, 865769162L, 683305004L, -1116000209L, -1484762907L, 
-1911473764L, 785236978L, -828719507L, 1094331303L, -1731127858L, 
1810563880L, -1844934357L, -2092461239L, -311773096L, -107504378L, 
2038659297L, -1943605133L, -722121502L, -1598317484L, -546021033L, 
279445837L, -418114332L, 1514694858L, -1545648843L, -821318769L, 
1712353430L, 1833529408L, 6215939L, 1537717537L, -166636016L, 
673203694L, 1636088217L, 1634465323L, 130514714L, 253961692L, 
-9188385L, 69237685L, 784868684L, -1980771422L, 927079165L, -582736905L, 
-1455203042L, 1378099160L, -923520709L, 491366105L, 42556008L, 
-721802346L, -942780655L, 2118826179L, -2114860846L, 505408868L, 
2033791143L, 410265085L, 1721956340L, -621444902L, -1065725403L, 
671207231L, 673610950L, 1028763984L, -864657037L, -107132015L, 
1277050304L, 1302755230L, -1112290967L, 677415323L, -657089366L, 
1101335244L, 670552975L, -753702971L, -330896516L, 507753426L, 
1567826765L, 329949703L, 825319918L, -204800824L, -1832758325L, 
998864425L, 642153592L, 1235663270L, -548791359L, -230811565L, 
303648514L, -686138828L, 29962167L, -1167904851L, 1311640772L, 
402440426L, 8050517L, -735098641L, -2073985290L, 1315473952L, 
-1778966557L, 1037292673L, -496524560L, -846993010L, -1808251143L, 
-1339190133L, 2038545594L, -170190916L, 395494655L, 2003241045L, 
293763948L, 1485132738L, 1007359005L, -737524585L, 2085687806L, 
695038392L, 1372306843L, 644563961L, -432723128L, -1827854730L, 
-322805839L, -1864957469L, 1427277490L, -896009468L, -808957881L, 
-2018910435L, -158153324L, 125667642L, 379671941L, 1089191519L, 
-1799747354L, -1971156112L, -1528763245L, -904234703L, 1481334752L, 
47768702L, -1712183543L, -357594821L, 184495370L, -2084240148L, 
-1129770641L, 1033590437L, -1803220900L, 352518194L, 614061869L, 
2091922919L, -231692018L, 2099302120L, -1127956245L, 1271644553L, 
-1152753384L, -2029890490L, -1729962847L, -2146360525L, -1451545950L, 
-1501142636L, 2013958167L, 441755661L, 1759288612L, -935892726L, 
1731854453L, -1659581489L, 240082134L, 914140288L, 1016099523L, 
-1188062111L, 752557264L, 1707879982L, -710299687L, 1306578923L, 
12854874L, 533225500L, 229530911L, 1486642549L, -2045178868L, 
1754697186L, 476475069L, 577300791L, 1596632542L, 159097624L, 
1598646907L, 942802969L, 1524154536L, -247143850L, -1256149295L, 
1373683331L, 1191908626L, -695066460L, 464371815L, -1086918211L, 
-148710604L, 1343893658L, -524822811L, 1413529727L, -49333754L, 
-1988893296L, 1666537523L, -543797216L, 1655197884L, -1246410160L, 
-1851033630L, 9100904L, -695529908L, 1640382716L, -2093611022L, 
502129920L, 65753972L, -1016459256L, -85585478L, -1322078768L, 
-604703764L, -542921836L, 874556370L, -1428289696L, 65085772L, 
298086688L, 1428496258L, 428656360L, -2047526676L, -738781636L, 
1534841714L, -1353378192L, -982277084L, -178913480L, 1123055770L, 
-1829953840L, -1866326852L, 686885524L, 1118261890L, 768975040L, 
399529276L, -1484081008L, 1199715874L, 93303304L, -1561543412L, 
2004374172L, 642467346L, -812299904L, 1603148308L, 897044008L, 
-1060310726L, 1768782800L, 138870380L, 984437556L, 272503442L, 
296166752L, 1236213068L, 259740640L, -2045972574L, 1004011432L, 
-1327216756L, -287126852L, -308489358L, 1428941040L, -1208941500L, 
251151704L, 1914847290L, 295490416L, 1567732796L, 743676692L, 
-5904958L, -591176480L, 1463025532L, -974072112L, 870361762L, 
302580776L, 215134476L, -648832068L, 739384626L, -127580224L, 
677353012L, 1100850120L, 25574586L, -1995358320L, -1805319700L, 
-1452477356L, -675804526L, 858237728L, 369863116L, 763057120L, 
758664898L, -1483024856L, -1804005716L, 1776565948L, 1381007858L, 
132198704L, 1970776228L, -201647432L, 1705908058L, -1319219824L, 
-230285956L, 219040020L, 8291650L, -568511744L, 1390453692L, 
1422979920L, -1431823070L, -121239224L, 112172044L, 666989660L, 
1631323346L, 1045427200L, 934266708L, -565167640L, -1112742918L, 
1478312656L, 50990380L, -171621132L, 1300468626L, -1741712288L, 
-1687758580L, 1044632288L, 1943248226L, -985467224L, 493518412L, 
1347310460L, -1151566030L, -1585690000L, -1763995324L, -1842018344L, 
-1339255686L, -1659049424L, 1790867132L, -1317664940L, -245098942L, 
314529696L, 1615258812L, 160944592L, 1608511842L, 1021759464L, 
-1426926644L, 47445756L, -1333734414L, -1101697920L, -1114527500L, 
-1606180728L, 1602441530L, -127956144L, 930173036L, -347640684L, 
-409135790L, -197010720L, 1137188044L, -222788448L, 1119470978L, 
215577832L, -823008020L, 774457404L, -590762894L, -1125618576L, 
-799747548L, 381580856L, -425786854L, 1354050768L, 493792828L, 
-822187884L, 334836866L, 1793971136L, -2124666564L, 1278840464L, 
-712117214L, -1965359096L, 2114716812L, 640482332L, -2108190574L, 
-1839940480L, -1167225836L, 641984296L, 1576767546L, 330930384L, 
548839532L, -429566284L, 1189264402L, 70710240L, 1713849292L, 
1905699296L, 141522978L, 541297960L, -1900199668L, 378053180L, 
-1932908302L, 112307440L, -899793084L, -437887912L, 1536114234L, 
-2005190416L, -1570273604L, -1870626284L, -1190132542L, -1305634208L, 
1657864572L, -1417749296L, -2054418398L, -1668672344L, 103608716L, 
1253134780L, 1881787442L, -355496128L, 1146636980L, -605921976L, 
1058838202L, 660212496L, 204216556L, -43529132L, 84816658L, -1497920096L, 
1237657548L, -865928352L, 1890553026L, -1452412376L, -352625236L, 
1500410684L, 1002651506L, -1812271312L, -1291891676L, -469298632L, 
-289466662L, 499114896L, 2133991676L, 636430740L, 1145760578L, 
1670656512L, -2081131460L, -1728279984L, -798330029L, -2054628748L, 
-17531902L, 1106210455L, 1076287585L, 60552358L, 1420551332L, 
608413013L, 1803390463L, -1256463208L, -1975193162L, 1045547587L, 
-992809083L, 1219767330L, -2028601584L, -1913095495L, -682168469L, 
-123941284L, 2087629706L, -865653889L, 285791945L, 172891998L, 
-266779636L, -1337515043L, -187119897L, 147171088L, -1970178514L, 
-812307237L, -2015919683L, 1354604330L, 1302939592L, -1420788495L, 
2075718243L, -1682401692L, -802464366L, -187827609L, 1017287217L, 
1209712502L, 869573780L, -2109433179L, -1538336849L, -1567382424L, 
810108998L, 2140634611L, 262650197L, -812790990L, 1870114560L, 
-764920855L, 1133889115L, -618440532L, 542592698L, 2117028975L, 
1767852537L, 383904910L, -1390033572L, 1182954957L, 1142470967L, 
1577001920L, 1202330782L, -731119253L, 1715402989L, -1328946854L, 
1062763032L, 1565624193L, -301972173L, -429577388L, -966902430L, 
-1437278729L, 1724995137L, -1979166010L, -437681148L, -1052437579L, 
-1278530401L, -802559752L, 277356758L, -779220637L, -2010946267L, 
1406100162L, -2005626960L, 1954116569L, 792780107L, -1186243652L, 
1217692650L, -1433082849L, -551050583L, 240150078L, 2072484652L, 
-442772099L, 1311202247L, -30004560L, -1123641138L, -1615736325L, 
-1340964259L, 731020682L, 981468584L, -1696412911L, -204191229L, 
-1785324284L, 1959049522L, 1908942727L, -1073769007L, -1928853674L, 
1623217460L, -1271481275L, 1041966223L, -843697656L, -460329434L, 
1067914963L, 171853877L, -1031739502L, -325085344L, -360262967L, 
-1881328133L, 1668620748L, 1407322202L, -1213549233L, -1259553767L, 
-1452342738L, 75382844L, 1290626477L, -470854057L, 1806804896L, 
-70843138L, 1295996875L, -1621466675L, -2145372934L, -1648519368L, 
766018785L, 619915155L, -195740748L, 1525272770L, -976124713L, 
893006113L, 776689894L, 4066788L, 14821653L, -464543425L, 892431576L, 
-2037766922L, -920871165L, -1795729339L, -964167582L, -1194758192L, 
1555215609L, -1830699477L, -493180516L, -2057913910L, -1140917569L, 
1967016969L, 1125494046L, 1840320716L, 1941075101L, -502462169L, 
1081275216L, 1608887790L, 2065259163L, 1076776189L, -736275478L, 
-1940798328L, -1266959823L, 1274686627L, 801509540L, 1990295890L, 
953564327L, -694847247L, -1192218570L, 1673547988L, 1800119397L, 
-1659314577L, -1083292785L)
.repmat <-
function(v,n,p){
	if (p==1){M = cbind(rep(1,n)) %*% v}
	else { cat('!'); M = matrix(rep(v,n),n,(length(v)*p),byrow=T)}
	M
}
