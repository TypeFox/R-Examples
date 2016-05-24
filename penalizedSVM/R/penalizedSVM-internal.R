.calc.mult.inv_Q_mat2 <-
function(U, D_vec, A_vec, mat2=NULL, n.thr=500){

## use always Sherman-Morrison-Woodbury formula
#
# calculate  new.mat = inverse_Q %*% mat2
#   where Q: = U * D * U' + A  - quadratic matrix
# 				mat2 - matrix 
### --> by large Q --> problems!
# D, A diag matrix, input as diag vectors
	
# NOTE : in coments '*' means matrix multiplication! 	
	
	require(MASS)
	
	
	# if mat2 doesn't exist--> use the identity matrix with dim [nrow(U),nrow(U)]
	if (is.null(mat2)) mat2<- diag(1, nrow=nrow(U), ncol=nrow(U))
	
	# if the Q matrix is not large < n.thr invert it directly
	if (nrow(U) < n.thr ){ 
		Q = U %*% diag(D_vec) %*% t(U) + diag(A_vec)
		
		if (rank.condition(Q)$rank == nrow(Q) & nrow(Q)==ncol(Q) ) {
				inv_Q<-solve ( Q )
		} else{
				inv_Q<-pseudoinverse(Q)
		}
		
		# sometimes if rank of Q is much smaller as nrow(Q) --> inv_Q consists only NAs.
		# then use Sherman-Morrison-Woodbury formula
		if (all(is.na(inv_Q) ) == TRUE ) {
			nbw<-.calc.mult.inv_Q_mat2 (U=U, D_vec=D_vec, A_vec=A_vec, mat2=mat2, n.thr=0)
		}else{
			nbw<-inv_Q %*% mat2
		}
#		## test
#		summary(nbw)
#			
#		inv_Q<-ginv(Q)
#		nbw.ginv<-inv_Q %*% mat2
#		summary(nbw.ginv)
#		
#		
#		inv_Q<-pseudoinverse(Q)
#		nbw.p<-inv_Q %*% mat2
#		summary(nbw.p)
#		
#		
#		inv_Q<-solve(Q)
#		nbw.s<-inv_Q %*% mat2
#		summary(nbw.s)
#		
#		
#		###end of test 
		
	}else{
		# apply Sherman-Morrison-Woodbury formula
		# A = Q1
		# U = t(x1)
		# U' = t(U)
		# inv A = A^-1
		# D = diag(D_vec) 
			
		# now 2 cases:
		# 1) A = 0 matrix
		# 2) A is a diagonal matrix with probably some elements = 0
		#    set 0 -> 10^-8 
		
		if (all(A_vec == 0) ) {
			# case 1
			# A is a null matrix
			# inv(A + U*D*U') =  inv( U*D*U') =  inv(U') * inv(D) * inv(U) = t (inv(U)) * inv(D) * inv(U) !!!!!
			# use inv(t(B) ) = t(inv(B) )
			
			A<-diag(A_vec)
			D<-diag(D_vec)
			inv_U<-pseudoinverse(U)  # can not apply solve to a not-square matrix!!
			
			
			#inv_Q <- t(inv_U) %*% diag(1/D_vec)  %*% inv_U  
			# nbw<-inv_Q %*% mat2
			
			nbw <- t(inv_U) %*% ( diag(1/D_vec)  %*% (inv_U %*% mat2 ) )
			
		} else {
			# case 2
			
			# inv(A + U*D*U') =  inv(A)- inv(A)* U * inv( inv(D) + U'*inv(A)*U ) * U' * inv(A)
			#term1 -  term2
						
			# inverse of a diagonal matrix = inv( diag(d1, d2,...) ) = diag(1/d1, 1/d2,...)
			
			# now for nbw
			#  nbw<-inv_Q %*% mat2  
			# -->  nbw<- term1 %*% mat2 - term2 %*% mat2 = nbw.term1 - nbw.term2 
			
			# matvec(M,v) is equivalent to M %*% diag(v) but is faster to execute
			# vecmat(v,M) is equivalent to diag(v) %*% M but is faster to execute. 
			# DO NOT need to storage a HUGE diagonal matrix
			
			# A is a diagonal matrix with some zeros on the diagonal. replace 0 with eps
			# define new A: A1 := A+ eps*I for A elements =0, eps = 10^-8
			
			eps<-10^-8 
			A1_vec<-A_vec
			A1_vec[A1_vec==0]<-  eps
			inv_A1_vec <- 1/ A1_vec
			
			# 1. nbw.term1 inv(A) *  mat2
			nbw.term1 <-  vecmat(inv_A1_vec, mat2 )
			
			# 2. nbw.term2  = inv(A)* U * inv( inv(D) + U'*inv(A)*U ) * U' * inv(A) * mat2
					# tt1 =  inv( inv(D) + U'*inv(A)*U )
					# nbw.term2 =   inv(A)* U * tt1 * U' * inv(A) * mat2
			
			# sometimes problem with  pseudoinverse  Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
			#tt1<-pseudoinverse ( diag(1/D_vec) + matvec(t(U), inv_A1_vec )%*% U )
			## Solution1: apply solve 
			if (exists("tt1")) {rm(tt1)} 
						
			tmp.mat<-diag(1/D_vec) + matvec(t(U), inv_A1_vec )%*% U 
			
			# if the matrix tmp.mat has the full rank --> use 'solve' 
			if (rank.condition(tmp.mat)$rank == nrow(tmp.mat) &  nrow(tmp.mat)==ncol(tmp.mat) ) {
					try(tt1<-solve ( tmp.mat ))
			}						
			
			# if 'solve' fails --> use 'pseudoinverse' / 'ginv'			
			if (!exists("tt1")) try(tt1<-pseudoinverse (tmp.mat ))
			
			if (!exists("tt1")) try(tt1<-ginv ( tmp.mat))
			
			if (!exists("tt1")) stop("Error: can not calculate inverse matix for 'diag(1/D_vec) + matvec(t(U), inv_A1_vec )%*% U'")
			 
			## Solution2: apply Sherman-Morrison-Woodbury formula again   --> doesn't make a lot od sense ,since we are back to the high dim! nrow(U)
			
			# to improve calculation's speed --> use ()
			nbw.term2 <- vecmat(inv_A1_vec, U) %*% ( tt1 %*% ( matvec(t(U), inv_A1_vec) %*% mat2 ) )
			
			nbw<- nbw.term1 - nbw.term2
		} # end of the case 2
	
	} #	end of if (nrow(U) < n.thr )
	return(nbw)
}



`.find.inverse` <-
function(U,D_vec,A_vec, n.thr=500){
### invert a quadratic matrix Q: = U * D * U' + A  
# D, A diag matrix, input as diag vectors
	
	# if the Q matrix is not large < n.thr invert it directly
	if (nrow(U) < n.thr ){
		Q = U %*% diag(D_vec) %*% t(U) + diag(A_vec)
		inv_Q<-pseudoinverse(Q)
	}else{
		# apply Sherman-Morrison-Woodbury formula
		# A = Q1
		# U = t(x1)
		# U' = t(U)
		# inv A = A^-1
		# D = diag(D_vec) 
		
		# now 2 cases:
		# 1) A = 0 matrix
		# 2) A is a diagonal matrix with probably some elements = 0
		#    set 0 -> 10^-8 
		
		if (all(A_vec == 0) ) {
			# case 1
			# A is a null matrix
			# inv(A + U*D*U') =  inv( U*D*U') =  inv(U') * inv(D) * inv(U) = t (inv(U)) * inv(D) * inv(U) !!!!!
			# use inv(t(B) ) = t(inv(B) )
			
			A<-diag(A_vec)
			D<-diag(D_vec)
			inv_U<-pseudoinverse(U) 
			
			# inv_Q <-  matvec(  inv_V, D_vec ) %*% inv_U  
			inv_Q <- t(inv_U) %*% diag(1/D_vec)  %*% inv_U  
			
		} else {
			# case 2
			
			# inv(A + U*D*U') =  inv(A)- inv(A)* U * inv( inv(D) + U'*inv(A)*U ) * U' * inv(A)
			#term1 -  term2
			#
			# inverse of a diagonal matrix = inv( diag(d1, d2,0,...) ) = diag(1/d1, 1/d2,0,...)
			
			# matvec(M,v) is equivalent to M %*% diag(v) but is faster to execute
			# vecmat(v,M) is equivalent to diag(v) %*% M but is faster to execute. 
			# DO NOT need to storage a HUGE diagonal matrix
			
			# A is a diagonal matrix with some zeros on the diagonal.
			# define new A: A1 := A+ eps*I, eps = 10^-8
			
			eps<-10^-8 
			A1_vec<-A_vec+ eps
			inv_A1_vec <- 1/ A1_vec
			
			# tt1 =  inv( inv(D) + U'*inv(A)*U )
			tt1<-pseudoinverse ( diag(1/D_vec) + matvec(t(U), inv_A1_vec )%*% U )
			# term2 = inv(A)* U * inv( inv(D) + U'*inv(A)*U ) * U' * inv(A)
			term2<- vecmat(inv_A1_vec, U) %*% tt1 %*% t(U)
			term22<- matvec( term2, inv_A1_vec  )
			inv_Q<- diag(inv_A1_vec) - term22
		}
	}
	
	return(inv_Q)
}

`.Random.seed` <-
c(403L, 200L, -764817312L, -1954253149L, -1164856895L, 316951984L, 
-1368292146L, -833293767L, 813040715L, 1187268090L, -2126855940L, 
-1528015553L, 980298261L, 2114806060L, 1882116610L, 1578672349L, 
557364183L, 44125886L, 383924472L, -238418213L, 709803321L, -794998648L, 
320966710L, 425321585L, 36900515L, -553884686L, -25831612L, 959010823L, 
1602254301L, -617991468L, -578667014L, -310857403L, -1244703329L, 
1251354406L, 1638377392L, -1576205229L, 1798930289L, -2009385312L, 
480647742L, 301723977L, 449023995L, -766094390L, -2080936660L, 
823025967L, -1250127899L, 1980812956L, -2075108622L, 423910253L, 
771155111L, 1798700750L, -653118424L, -628881365L, 1657948745L, 
256418648L, -985027066L, 543045601L, -1597667981L, 1200106978L, 
-901848748L, 399790167L, -694846387L, 1439182820L, -1797890614L, 
2004739125L, -1052040561L, -1953604714L, 1346692928L, -1458664957L, 
-1661997023L, 1979030288L, 801290478L, -556681063L, -1002934997L, 
2020066330L, 292566748L, -447963425L, -1498097483L, 755871308L, 
1642476194L, 1183955453L, 710644471L, 1172141086L, 89374424L, 
1903690299L, -1150894119L, 1290298728L, -364006250L, -976129519L, 
1852853187L, 1214191570L, 878404196L, 9904551L, 431611645L, -2084222220L, 
-2070094374L, 431649061L, -800929729L, 289335238L, 1999879248L, 
82972787L, 133579921L, 1051784384L, -1626720610L, -186840471L, 
-1930745189L, -1642559574L, -1975389236L, 112725647L, 1539395269L, 
1194370684L, 812387026L, -150790579L, 1057911559L, 1444148974L, 
-1678789688L, 2010784971L, 594836777L, 2048789880L, 1016836774L, 
-1484516159L, 1966499155L, 1891492354L, -787965132L, -143148873L, 
2086600877L, 991124932L, 2137755114L, 149074517L, 1075182575L, 
-1404150794L, -1056707296L, 1446232803L, -1626713215L, 1383187952L, 
77577358L, 1321599993L, -1221795445L, -675660358L, 1652425916L, 
1543550975L, -1927391915L, -1365588372L, -992447806L, 1657189661L, 
-1939240041L, -1057587458L, 1763589304L, -1568428389L, -89817863L, 
-653935032L, 892582262L, 1526571185L, 478319843L, 1849942962L, 
-318022652L, 2048515911L, -1638992355L, -492827500L, -227689414L, 
2103883397L, -1733458081L, -430258202L, 714074736L, -581976685L, 
498338353L, -943036192L, -1083096706L, -436493303L, 350381115L, 
1825966602L, 1617007084L, -2071804305L, 1163438501L, 261427548L, 
-576541390L, 1778136109L, 741001959L, -1050288114L, 510792168L, 
-752307221L, -1221489527L, -1052051432L, 1839433542L, 1820644769L, 
-494515149L, 547471778L, 1326527636L, -2016410857L, 1603964173L, 
-159749084L, 138107914L, -1579327115L, -806004529L, 730381270L, 
-1670564992L, -1350457405L, 1236767073L, -2106697776L, -605495506L, 
747840217L, 1079812331L, -350924454L, 378010396L, 2121272863L, 
-1558694283L, 1236042508L, 2126028514L, -34868291L, 1414833207L, 
-170318114L, -1373736424L, -1785365125L, 1068315929L, 763247016L, 
1735421270L, -780603951L, 2013259139L, -167257070L, 34707364L, 
963653479L, 1226755261L, 1223175220L, -186228838L, 141671397L, 
1610299263L, -9943290L, -1131949936L, 1720035635L, -934658351L, 
-1942333252L, 490668624L, -2105054238L, -428898200L, -1169417140L, 
2043672828L, 1991236594L, -1287554816L, -1465491084L, 1439540232L, 
958742458L, -269708848L, 1258519532L, -667928172L, 1498202578L, 
-894996128L, 1645162316L, 1795598624L, 847837058L, 2059775720L, 
291315948L, -2080381892L, 6183794L, 1397655152L, -1701032924L, 
-1485528264L, 1602804890L, 1659742928L, -296722244L, 1057439380L, 
-831459710L, -1011296064L, 239827772L, 2137918608L, 1367416354L, 
316252168L, 1883915532L, -776174436L, 204430866L, 2137722752L, 
-1154423276L, -854425048L, -132305606L, -1605647920L, -2112996244L, 
-887039692L, 994913938L, -1243098272L, -1735236788L, 202815456L, 
48551842L, -1175685208L, -686060148L, -1035407172L, -1170452110L, 
-2091692816L, -1137897404L, -1936115368L, -1587872198L, -1316873360L, 
-409921476L, -622602476L, -2072953406L, 283644640L, 527560572L, 
1461231312L, -1727696734L, -919624664L, 621921548L, 1635471292L, 
-1314016974L, -1071058496L, 1551661620L, 1676506056L, -190902086L, 
-1169954416L, 1636059628L, -705381804L, 1966581906L, -1076680928L, 
212586956L, -1596439072L, -1044695358L, -107414488L, -1246209364L, 
1604023484L, 1481777650L, -764206800L, 286620836L, -331231560L, 
336997722L, 1109486992L, -1803986052L, -1225031404L, -509521086L, 
-711750912L, 824948668L, -130930864L, 2140107554L, -1088359096L, 
-1147878900L, 2006120540L, -527706926L, 869163520L, -6875820L, 
-1294671384L, 171232250L, -1062801200L, 1292496172L, 1365967604L, 
715108754L, -1607195552L, -2001801972L, -180469024L, 435384674L, 
-451339608L, 1905238604L, 811321724L, 1967488306L, -1518526352L, 
-871579836L, 954586584L, 682294394L, -1582133712L, 347247292L, 
-1826523308L, 131978306L, 1612717472L, -1786639684L, -480374320L, 
-935668382L, -2065563672L, 1113948108L, -1036777732L, -1107895822L, 
1889849472L, 451763444L, -1530306936L, -32207046L, 1889925968L, 
-481501588L, -800122732L, -833492654L, 521066208L, 827923148L, 
267038880L, -329622654L, 1187238632L, 611005164L, -1385329604L, 
821381746L, -409942928L, -794206684L, -169460680L, -1734121446L, 
16538320L, -1960596420L, -1908988780L, -1973527422L, 1411137472L, 
1523675452L, 1622555280L, 292054050L, 1928119304L, -337401716L, 
1703654428L, 1919938194L, -1344985984L, 818688020L, 1404452136L, 
1156542010L, -361673520L, -491232148L, -553680204L, 407293970L, 
176707552L, 1468520396L, 1898151392L, -1725934046L, 1062437160L, 
399300876L, -456251844L, -1875185422L, -1966046480L, 310543172L, 
968323672L, 494426682L, 79298288L, 1279984828L, 1677432340L, 
-1433929534L, -1669867424L, -638118532L, 1105269968L, 1677687842L, 
651395752L, -1644447860L, -135884356L, 1532646450L, 29882176L, 
957757620L, -450262200L, -749700422L, -1304479984L, 614244588L, 
551101524L, 359419666L, 1702607776L, -517517364L, 1515424608L, 
1009127618L, -1186798040L, 1265952172L, 1272987452L, 1714575218L, 
1236331312L, 1263850020L, -1762084808L, 1195336922L, -1089959536L, 
-1701992708L, -1325404268L, 1219185474L, 1537171968L, -1968743876L, 
-162996656L, -301450590L, 347639668L, 924269826L, 390858135L, 
-272760991L, 668560294L, 1164199332L, -1529119659L, -464973569L, 
2024998296L, 1896918710L, -1452880061L, -1404681595L, -1110423774L, 
-590498800L, -319129671L, -1953466261L, -1231727780L, 1958912650L, 
-1681998209L, 888369097L, -1566718882L, -1059238132L, -1443729699L, 
12097511L, 252918800L, -1994895058L, -518089253L, 1685379773L, 
861487146L, 1725626568L, 1467932145L, -1641321629L, -1270158492L, 
437478034L, -1702846617L, 145076529L, -97494410L, -992277612L, 
-1608573019L, 187633839L, 224988520L, -1402753722L, 1466698995L, 
-1132064171L, 1125363250L, -183604224L, -41213719L, -2096736421L, 
722612652L, 1766805946L, -1142526097L, -881880839L, 61040014L, 
-2116085668L, 594295501L, -1493213641L, 530810048L, -586777698L, 
232835179L, -421828627L, 67136090L, -1537340136L, 559921793L, 
65136691L, -34288044L, 631663202L, 168115447L, -1281333439L, 
-1989302842L, 743077636L, 1638318773L, -328229473L, -949661192L, 
1273094102L, 1318517347L, -907148763L, -629274174L, -699658576L, 
-1128040743L, -444223413L, -2116956484L, -2102764310L, 555616543L, 
-153495127L, 119766846L, 1720420908L, -1891387267L, -1459844409L, 
-727518288L, -680581682L, -52226309L, 834522461L, 605406858L, 
-564706136L, 1720254481L, -85765885L, -1169114108L, -626330574L, 
-1392968569L, -1564903215L, 794581590L, 1525039668L, 937158981L, 
2056270735L, 983906056L, -737189594L, -685579821L, -168694987L, 
-1424905582L, 99389024L, -1395386423L, -413767941L, 2013948108L, 
-1484332710L, -2128720817L, -538746599L, -972481746L, 2007115068L, 
-1406360915L, -896402089L, 269163680L, 480393726L, -786545461L, 
-354659635L, -1471746566L, -2101871048L, -1658330143L, 565464211L, 
1369253556L, -1192682558L, 370515927L, -104558559L, -1776453658L, 
-139848988L, 776807957L, 1265638975L, -984764968L, -960433162L, 
168913923L, 1685102405L, 1397709666L, -1245174576L, 818409465L, 
1625183019L, -1631528804L, 137574090L, 1849487295L, 351812361L, 
-203934690L, -405586484L, -766486627L, -291544537L, 988835920L, 
-1443858194L, 629286299L, 1411410941L, 1104459498L, -414187128L, 
1522650929L, 710696867L, -49481820L, 980656722L, -1713374809L, 
1783868913L, 1095342390L, 1049817556L, -1426630811L, -1239169169L, 
-1846400600L, 2055394584L)

`.required` <-
c("corpcor", "statmod", "MASS", "e1071", "mlegp", "tgp", "lhs")

