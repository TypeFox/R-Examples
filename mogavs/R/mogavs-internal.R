.crossover <-
function(numberOfOffsprings, parent1, parent2, crossoverProbability=0.9){
  if(missing(numberOfOffsprings) | numberOfOffsprings==0)
    stop("Arg numberOfOffsprings is missing or zero")
  if(missing(parent1))
    stop("Arg parent1 is missing")
  if(missing(parent2))
    stop("Arg parent2 is missing")
  #Note that lambda should preferably be even. If it isn't then the
  #code produces an extra offspring but returns exactly lambda offsprings
  

  numberOfCrossovers<-ceiling(numberOfOffsprings/2)
  n<-ncol(parent1)
  
  if(!exists("crossoverProbability")){
    crossoverProbability<-1
  }
  
  #initialize offspring
  offspring1<-matrix(0,numberOfCrossovers,n)
  offspring2<-matrix(0,numberOfCrossovers,n)
  
  for(i in 1:numberOfCrossovers){
    a<-1+floor(runif(1)*nrow(parent1))
    b<-1+floor(runif(1)*nrow(parent2))
    crossoverSite<-floor(runif(1)*n)
    if(crossoverSite !=0 && runif(1)<=crossoverProbability){
      offspring1[i,]<-c(parent1[a,1:crossoverSite],parent2[b,(crossoverSite+1):n])
      offspring2[i,]<-c(parent2[b,1:crossoverSite],parent1[a,(crossoverSite+1):n])
    }
    else{
      offspring1[i,]<-parent1[a,]
      offspring2[i,]<-parent2[b,]
    }
  }
  
  offspring<-rbind(offspring1,offspring2)
  
  #If lambda is odd then the last offspring is eliminated
  offspring<-offspring[1:numberOfOffsprings,]
  if(!is.matrix(offspring)){offspring<-matrix(offspring,1,length(offspring))}
  return(offspring)
}
.evaluateOffspring <-
function(offspring,x,y){
  if(missing(offspring))
    stop("Arg offspring is missing")
  if(missing(x))
    stop("Arg x is missing")
  if(missing(y))
    stop("Arg y is missing")
  
  
  
  MSEChildren<-matrix(0,1,nrow(offspring))
  numOfVariablesChildren<-matrix(0,1,nrow(offspring))
  for(i in 1:nrow(offspring)){
    yy<-y
    xx<-cbind(matrix(1,nrow(x),1), x[,which(offspring[i,]==1)])
    beta_coeff_child<-qr.solve(xx,yy)
    
    numOfVariablesChildren[i]=sum(offspring[i,])
    MSEChildren[i]=sum((yy-xx%*%beta_coeff_child)^2)/length(yy)
    
    
  }  
  output<-list(numOfVariablesChildren,MSEChildren)
  return(output)
}
.initializePop <-
function(popSize, x, y){
  
  if(missing(popSize))
    stop("Arg popSize is missing")
  if(missing(x))
    stop("Arg x is missing")
  if(missing(y))
    stop("Arg y is missing")
  
  
  sampleSize<-nrow(x)
  #N is the number of independent variable
  numOfVariables<-ncol(x)
  
  
  obj1<-matrix(0,1,max(numOfVariables,popSize))
  obj2<-matrix(0,1,max(numOfVariables,popSize))
  popMembers<-matrix(0,max(numOfVariables,popSize),numOfVariables)
  
  correlationCoefficient<-cor(y,x)
  sortedCorr<-sort(abs(correlationCoefficient),decreasing=T,index.return=T)
  
  for(i in 1:numOfVariables){
    popMembers[i,sortedCorr$ix[1:i]]=1
  }
  
  for(i in 1:numOfVariables){
    yy<-y
    xx<-cbind(matrix(1,nrow(x),1),x[,which(popMembers[i,]==1)])
    beta_coeff<-qr.solve(xx,yy)
    
    obj1[,i]=i
    obj2[,i]=sum((yy-xx%*%beta_coeff)^2)/length(yy)
  }
  
  if(popSize<numOfVariables){
    popMembers<-popMembers[1:popSize,]
    obj1<-obj1[,1:popSize,drop=FALSE]
    obj2<-obj2[,1:popSize,drop=FALSE]
  }
  if(popSize>numOfVariables){
    popMembers[(numOfVariables+1):popSize,]<-popMembers[1:(popSize-numOfVariables),]
    obj1[,(numOfVariables+1):popSize]<-obj1[1:(popSize-numOfVariables),drop=FALSE]
    obj2[,(numOfVariables+1):popSize]<-obj2[1:(popSize-numOfVariables),drop=FALSE]
  }
  output<-list(popMembers,obj1,obj2)
  return(output)
}
.mutation <-
function(offspring,mutationProbability=1/ncol(offspring)){
  if(missing(offspring))
    stop("Arg offspring is missing")
  
  numOfOffspring<-nrow(offspring)
  numOfVariables<-ncol(offspring)
  
  if(!exists("mutationProbability")){
    mutationProbability<-1/numOfOffspring
  }
  
  randomMatrix<-matrix(runif(numOfOffspring*numOfVariables),numOfOffspring,numOfVariables)
  offspring[which(randomMatrix<=mutationProbability)]<-1-offspring[which(randomMatrix<=mutationProbability)]
  return(offspring)
}
.nonDomination <-
function(obj1Members, obj2Members, members){
  if(missing(obj1Members))
    stop("Arg obj1Members is missing")
  if(missing(obj2Members))
    stop("Arg obj2Members is missing")
  if(missing(members))
    stop("Arg members is missing")
  #allMembers contain [nonDominated;dominatedMembers]
  #The entries with value 1 in the vector named dominatedBoolean will be the dominated solutions
  #The entries with value 0 in the vector named dominatedBoolean will be the non-dominated solutions
  
  noOfMembers<-nrow(members)
  dominatedBoolean<-matrix(0,noOfMembers,1)
  
  for(i in 1:noOfMembers){
    if(dominatedBoolean[i]==0){
      for(j in 1:noOfMembers){
        if(dominatedBoolean[j]==0){
          if((obj2Members[i]<=obj2Members[j]) && (obj1Members[i]<=obj1Members[j]) && ((obj2Members[i]<obj2Members[j]) || (obj1Members[i]<obj1Members[j]))){
            dominatedBoolean[j]<-1
          }
        }
      }
    }
  }
  
  obj1AllMembers<-obj1Members
  obj2AllMembers<-obj2Members
  
  l<-which(dominatedBoolean==0)
  nonDominatedMembers<-matrix(0,noOfMembers-sum(dominatedBoolean),ncol(members))
  if(nrow(nonDominatedMembers)>0){
    for(i in 1:nrow(nonDominatedMembers)){
      nonDominatedMembers[i,]<-members[l[i],]
      obj1AllMembers[i]<-obj1Members[l[i]]
      obj2AllMembers[i]<-obj2Members[l[i]]
    }
  }
  
  l<-which(dominatedBoolean==1)
  dominatedMembers<-matrix(0,sum(dominatedBoolean),ncol(members))
  if(nrow(dominatedMembers)>0){
    for(i in 1:nrow(dominatedMembers)){
      dominatedMembers[i,]<-members[l[i],]
      obj1AllMembers[i+nrow(nonDominatedMembers)]<-obj1Members[l[i]]
      obj2AllMembers[i+nrow(nonDominatedMembers)]<-obj2Members[l[i]]
    }  
  }
  
  if(nrow(dominatedMembers)>0 && nrow(nonDominatedMembers)>0){
    allMembers<-rbind(nonDominatedMembers,dominatedMembers)
  }
  if(nrow(dominatedMembers)==0 && nrow(nonDominatedMembers)>0){
    allMembers<-rbind(nonDominatedMembers)
  }
  if((nrow(dominatedMembers)>0) && (nrow(nonDominatedMembers)==0)){
    allMembers<-rbind(dominatedMembers)
  }
  output<-list(nonDominatedMembers,dominatedMembers,allMembers,obj1AllMembers,obj2AllMembers)
  
  return(output)
}
.Random.seed <-
c(403L, 65L, 24604460L, 2084420359L, -1079130330L, -500643904L, 
1927520024L, 1838257364L, 1381563235L, -1907199619L, -1561983429L, 
1290012267L, 1837777274L, -312183753L, -502523975L, -397865195L, 
1163424760L, 831396104L, -807374424L, -822873510L, -1308131103L, 
729751707L, -246010754L, 2077830194L, -917643598L, 1832113436L, 
1705480766L, 174447329L, 629431613L, -891048624L, 1466640918L, 
-862267919L, 1351270038L, -222710985L, 115866080L, -1664629345L, 
781279902L, -1293560639L, -608745959L, 483497699L, 272240467L, 
-918925548L, 1256733538L, 1024869540L, -70437177L, -1505351443L, 
1464818361L, -254162078L, 741747227L, 2124712367L, 2030387942L, 
250200898L, -1506870209L, 1719684348L, 102243001L, 1742091118L, 
1327702021L, 223665680L, 696484620L, -1599185314L, 1577946923L, 
641449352L, 1029624718L, -1819597843L, -48554678L, 1555442242L, 
-2035622845L, 743429554L, -1439689126L, -802885495L, -301278868L, 
863155810L, 1325591627L, 1307547872L, -1981860299L, -1734273103L, 
-707387251L, 235531675L, 1106852279L, -876427865L, 2142640007L, 
-1699256296L, -994165720L, 1419039395L, 1616988387L, 1023594598L, 
78027264L, -785350884L, -403852674L, 109521421L, 571609447L, 
1835802340L, 1345168490L, -1727041541L, 1010879367L, 1377951706L, 
-2139764846L, 837127417L, 596679826L, -1244529048L, 1990692043L, 
655578531L, 423283527L, -695138861L, 1210278074L, -1088180231L, 
-303999327L, -1133000712L, 1298725637L, 1124134182L, -1477100475L, 
-1903748500L, -1599420936L, 1275697478L, -169585136L, -2121231965L, 
-1266516877L, 1395395497L, -418374911L, -226860782L, -1443017301L, 
-287238066L, 1224447154L, -1280277112L, 1235422252L, 1491197837L, 
940695467L, -908428482L, -1755446852L, -75901246L, -1433413895L, 
859282077L, 689562985L, 1456391202L, -1877894207L, -1083364379L, 
77502789L, 1422411690L, -516041510L, 1726543347L, -1911662005L, 
906276497L, -348279107L, 405397583L, -243484800L, -835972280L, 
-1454435058L, -1130501930L, 624352614L, 202296343L, -1751492068L, 
486888617L, -267770779L, -1946088389L, 828822708L, 1112401913L, 
-39448164L, -1698230122L, 1583912069L, -732817029L, 577126739L, 
-24570346L, -1879563550L, 1432124821L, -1616310711L, -269305329L, 
1799701242L, 972452489L, -369377431L, 404100753L, -472642086L, 
1989562523L, -376995841L, 725164464L, 949063120L, 1415684982L, 
-2053588553L, 479159122L, 1447778159L, -25006718L, -225092060L, 
-1396480004L, -681609781L, 1618155714L, -1511125115L, -1227567721L, 
-398215777L, -1366317527L, -1224204314L, -1611133431L, 1189082879L, 
-1061215083L, 1969062795L, 851971870L, -733104818L, -66260306L, 
1809437159L, 2112337509L, 1586578002L, -690114336L, 1464935413L, 
-1432076882L, -1662499573L, 1671631544L, -1756928524L, -987122280L, 
-729760146L, 1991314821L, 1771298362L, 915625350L, 159282789L, 
645381417L, 886897601L, -491811409L, -815761085L, -778390874L, 
87607196L, 2138851679L, -1286669344L, -1305408961L, 374790485L, 
-660617077L, 217043303L, -603987591L, 813071237L, -1991523890L, 
689465818L, -1372964163L, 1914284731L, 1325557217L, -926701711L, 
2013235359L, 299296934L, -759287225L, -179683592L, 2001348415L, 
474459361L, 2131603321L, -448812891L, 433397913L, -268121685L, 
1302104834L, 1975768128L, 165946213L, 2091654048L, -741640403L, 
257440665L, 1125482938L, -903074029L, -35821014L, -244822167L, 
205059374L, -1645405242L, -1113667614L, 1115203678L, -2105465584L, 
414295577L, 89320122L, 190417050L, -1830930201L, 1564388514L, 
-443843836L, -775501629L, 263352294L, 401574020L, 1745180969L, 
512931490L, -1367926893L, -215178240L, -2110708766L, 542460387L, 
1082706531L, -1622329553L, 1848889841L, 1412047881L, 742515093L, 
1667494708L, 933550953L, 339050729L, -1024727129L, -2043152185L, 
2018119734L, -318507670L, -783605700L, -140563557L, 1296207167L, 
-912462919L, 324150634L, 1802131135L, 900452975L, 1622355720L, 
-1188834617L, -963506112L, -1321115730L, -1145762811L, -630492120L, 
211629946L, 492342374L, -505015889L, 1016125006L, 1792825088L, 
-1133835938L, -755604787L, -385366884L, 360652834L, 831356332L, 
-1326694728L, 103536600L, 551734003L, -10643210L, 1936940589L, 
141042452L, -879192685L, -1800652935L, -1048139939L, 1839921179L, 
925420689L, 992352531L, 1808145204L, 652880444L, 150067042L, 
-63092821L, 1289957814L, -447575865L, -2027928374L, 1980223655L, 
254596599L, 1260946562L, -80342149L, -195193104L, -1144883839L, 
-1038832274L, -1461490331L, 575807779L, -1075842769L, -1040705619L, 
658658078L, -1634921424L, 1695011313L, -1867114304L, -1656100836L, 
-955503903L, 1113641304L, 161375544L, -1569488934L, -1700391122L, 
1675096836L, -1514129750L, -1939558356L, 1832101256L, -1566639570L, 
828877549L, 905031177L, 2038694527L, 971017589L, 1993679749L, 
525235195L, 1591242911L, -742139813L, 1667242462L, -2129205253L, 
1479777227L, -1062101093L, 579517076L, -1497804529L, 513914306L, 
589559182L, 1796404184L, 2116534157L, 999146709L, -1585510506L, 
-1746242590L, -1053058746L, -1478403834L, -1046160032L, 539604899L, 
-1527413063L, 1234123974L, -1200864535L, 1237898679L, 1093133999L, 
-405348144L, 1418868774L, -1542724041L, 54601556L, 267022501L, 
121103006L, 1900771337L, 159706882L, 1040770478L, 1399219263L, 
-1687580760L, 1599145416L, -55193567L, 2147371120L, 141405645L, 
-73980696L, -195091869L, -655835729L, -368193542L, -771638442L, 
-775536523L, -815785368L, 1731925503L, 1227350975L, 404809359L, 
-1536939461L, -1363557275L, -274938436L, -35241078L, 647859754L, 
-1661027045L, 1397028213L, -848288656L, 1832876140L, -1364017109L, 
1239475993L, -275712276L, -1523125334L, -1775666133L, -1814976192L, 
-674102332L, -406167400L, 950076535L, -913812808L, -407993545L, 
-388493354L, -1299040043L, -1383651102L, -1707389744L, 1085820263L, 
-458848596L, -1919521485L, 37692906L, 352614152L, 2039419599L, 
1163660302L, 570773932L, 1030181299L, 183949657L, 134880233L, 
-719686448L, 238845211L, 929691727L, 1378802825L, -1054113412L, 
2106745686L, 1381208050L, 2059605483L, 280130588L, -935256704L, 
-2045696174L, 638944030L, -135184809L, -718499052L, -687869258L, 
645292772L, 954428755L, -409608971L, -385835412L, 337115797L, 
-1730079278L, -572575893L, 372196487L, 1041690277L, -646373625L, 
-235247560L, -807089237L, 2122721678L, 119457174L, 52057432L, 
12883294L, -298850192L, 390459534L, -1459918488L, 757813118L, 
-722768264L, 1594335334L, -250377693L, -1193937759L, -1407760241L, 
1532621024L, -1498853826L, 771955793L, -1249935398L, -913438615L, 
1095700134L, 1377467394L, -860806838L, 1046068215L, 25253168L, 
1996391234L, -191094365L, -553085047L, 2076157561L, 2015267872L, 
150485115L, 2074779540L, 213360481L, -681980327L, 1287851158L, 
1108273184L, 979720574L, -28142959L, 2046179248L, -1367208412L, 
657643331L, 1833678275L, 360418113L, 12268263L, -1655909631L, 
559033692L, 93004187L, 1152366901L, 1847868281L, 1894464301L, 
-471830219L, -1812383453L, 1957783380L, 1058082028L, 342410051L, 
-1851098201L, -1396561790L, -1097054144L, -1681451518L, -28448473L, 
1565187565L, 1925067454L, -1359205319L, 1661020484L, 853222308L, 
-1874475159L, -2094363134L, -1187679972L, 273498564L, 513996393L, 
-1201204110L, 1088442906L, 1255378381L, -938382330L, -1859547650L, 
1392380802L, 1453502868L, -361347465L, 625485664L, -1336038448L, 
-320359036L, 2062479173L, 499985580L, -834645143L, 1628815895L, 
1222133530L, 378722653L, -59150227L, -570894111L, 2104625682L, 
683610662L, -1010670273L, -394068767L, 711719280L, -1927302404L, 
1104930333L, 1595679464L, 1327587183L, -401146485L, -1444330265L, 
1169718869L, -1435735817L, -1522047040L, 1372954985L, 1589257835L, 
-642846342L, 88310125L, 1852321190L, 292880901L, 263402855L, 
1761111080L, 658270484L, 1225301316L, -1993602560L, -939960588L, 
-1622747735L, 1474294783L, -932057258L, 1250870764L, 486358072L, 
1283727760L, 805322789L, -1112365054L, -1548975092L, 1798699182L, 
686444897L, 415465575L, 1465876344L, -938828663L, 9418064L, 244046574L, 
-676713965L, 1672788432L, 1863299107L, -1424375532L, 16095288L, 
-1493206408L, -674571790L, 1786052064L, 1823599036L, -1855083395L, 
-321459838L, 1409039813L, 2004852403L, 96684943L, -791876812L, 
-404818183L, 170571065L, 211102239L, -1934339211L, -1035996233L, 
-414580228L, -2144541755L, 663579081L, 696316637L, 896374719L, 
-764888473L, -1654726514L, -91454727L, -1739388359L)
.removeDuplicates <-
function(vector,obj1Vector,obj2Vector){
  if(missing(vector))
    stop("Arg vector is missing")
  if(missing(obj1Vector))
    stop("Arg obj1Vector is missing")
  if(missing(obj2Vector))
    stop("Arg obj2Vector is missing")
  
  I<-!duplicated(vector,MARGIN=1)
  newVector<-vector[I,]  
  newObj1Vector<-obj1Vector[I]
  newObj2Vector<-obj2Vector[I]
  
  output<-list(newVector,newObj1Vector,newObj2Vector)
  return(output)
}
.sortByComplexity <-
function(members,obj1,obj2){
  if(missing(members))
    stop("Arg members is missing")
  if(missing(obj1))
    stop("Arg obj1 is missing")
  if(missing(obj2))
    stop("Arg obj2 is missing")
  
  
  complexity<-apply(members,1,sum)
  
  sortedComplexity<-sort.int(complexity, index.return=T)
  orderedMembers<-members[sortedComplexity$ix,]
  if(!is.matrix(orderedMembers)){orderedMembers<-matrix(orderedMembers,1,length(orderedMembers))}
  newObj1<-obj1[sortedComplexity$ix]
  newObj2<-obj2[sortedComplexity$ix]
  output<-list(orderedMembers,newObj1,newObj2)
  return(output)
}
.updatePopulationMembers <-
function(nonDominatedSet, dominatedSet, popMembersAugmented, obj1Augmented, obj2Augmented, popSize){
  if(missing(nonDominatedSet))
    stop("Arg nonDominatedSet is missing")
  if(missing(dominatedSet))
    stop("Arg dominatedSet is missing")
  if(missing(popMembersAugmented))
    stop("Arg popMembersAugmented is missing")
  if(missing(obj1Augmented))
    stop("Arg obj1Augmented is missing")
  if(missing(obj2Augmented))
    stop("Arg obj1Augmented is missing")
  if(missing(popSize))
    stop("Arg popSize is missing")
  
  #trueNonDominatedSet variable actually stores the (elite) non dominated members. The variable
  #NonDominatedSet later on takes some dominated members as well to form the next population members.
  
  trueNonDominatedSet<-nonDominatedSet
  
  while(nrow(nonDominatedSet)<popSize){
    temp<-.nonDomination(obj1Augmented[(nrow(nonDominatedSet)+1):length(obj1Augmented)],obj2Augmented[(nrow(nonDominatedSet)+1):length(obj2Augmented)],dominatedSet)
    nonDomTemp<-temp[[1]]
    domTemp<-temp[[2]]

    popMembersAugmented[(nrow(nonDominatedSet)+1):nrow(popMembersAugmented),]<-temp[[3]]
    obj1Augmented[(nrow(nonDominatedSet)+1):length(obj1Augmented)]<-temp[[4]]
    obj2Augmented[(nrow(nonDominatedSet)+1):length(obj2Augmented)]<-temp[[5]]
    rm(temp)

    if(!is.matrix(nonDomTemp)){
      print("nonDomTemp was not a matrix")
      print(nonDomTemp)
    }
    temp<-.sortByComplexity(nonDomTemp,obj1Augmented[(nrow(nonDominatedSet)+1):(nrow(nonDominatedSet)+nrow(nonDomTemp))],obj2Augmented[(nrow(nonDominatedSet)+1):(nrow(nonDominatedSet)+nrow(nonDomTemp))])
    if(!is.matrix(temp[[1]])){
      print("temp[[1]] was not a matrix")
      print(temp[[1]])
    }
    nonDomTemp<-temp[[1]]
    obj1Augmented[(nrow(nonDominatedSet)+1):(nrow(nonDominatedSet)+nrow(nonDomTemp))]<-temp[[2]]
    obj2Augmented[(nrow(nonDominatedSet)+1):(nrow(nonDominatedSet)+nrow(nonDomTemp))]<-temp[[3]]
    rm(temp)

    popMembersAugmented[(nrow(nonDominatedSet)+1):(nrow(nonDominatedSet)+nrow(nonDomTemp)),]<-nonDomTemp
    nonDominatedSet<-rbind(nonDominatedSet,nonDomTemp)
    dominatedSet<-domTemp
  }

  if(nrow(nonDominatedSet)<popSize){ #this will never happen if initial pop doesn't contain duplicates
    randomVector<-ceiling(nrow(nonDominatedSet)*runif(popSize))
    popMembers<-nonDominatedSet[randomVector,]
    
  } else{
    popMembers<-nonDominatedSet[1:popSize,]
  }
  
  if(popSize<=nrow(trueNonDominatedSet)){
    eliteMembers<-popMembers
  } else{
    eliteMembers<-trueNonDominatedSet
  }
  
  numOfVariables<-t(obj1Augmented[1:popSize])
  MSE<-t(obj2Augmented[1:popSize])

  output<-list(popMembers, numOfVariables, MSE, eliteMembers)
  return(output)
}
