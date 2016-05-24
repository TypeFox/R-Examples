removeRegionalBias = function(object, regionalBias, formulaString = value ~1, regCode = "regCode"){
  if ("regionalBias" %in% names(regionalBias)) regionalBias = regionalBias$regionalBias
  depVar = as.character(formulaString[[2]])
  uRegCode = unique(object[[regCode]])
  for (ic in 1:length(uRegCode)){
    rci = uRegCode[ic]
    cBias = regionalBias$wls[uRegCode == rci]
    cVar =  regionalBias$wlsvar[uRegCode == rci]
    cNew = object[[depVar]][object[[regCode]] == rci] - cBias
    object[[depVar]][object[[regCode]] == rci] = cNew
  }
  return(object)
}



    

findRegionalBias = function(object, boundaryLines, formulaString = value~1,
               minKrige = 5, regCode = "regCode", unbias = "default") {
  uRegCode = unique(object[[regCode]])
  gdat = regionalVariograms(object,uRegCode,formulaString,regCode = regCode)
  cDiff = regionalDiff(object,uRegCode,gdat,boundaryLines,formulaString,minKrige=5,regCode = regCode)
  dDiff = eval(call(paste("unBias.",unbias,sep=""),cDiff,uRegCode))
  cBias = dSolve(dDiff)
  dBias = list(regionalBias = data.frame(regCode = uRegCode,cBias),cDiff=cDiff)
  names(dBias$regionalBias)[1] = regCode  
  return(dBias)
}



unBias.default = function(cDiff,uRegCode) {
  D = cDiff$D
  Q = cDiff$Q
  V = cDiff$V
#
  D = rbind(D,0)
  cd = dim(D)[1]
  cd = cd + 1
  D = rbind(D,0)
  D[cd,] = 1
  Q[cd] = 0
  V[cd] = min(V)
  cDiff$D = D
  cDiff$Q = Q
  cDiff$V = V
  return(cDiff)
}



regionalDiff = function(object,uRegCode,gdat,boundaryLines,formulaString,minKrige=5,regCode="regCode") {
# Function to estimate the differences between all neighbouring countries.
# Function gives a kriging estimate on the border of each of the countries,
# and uses the difference between the estimates to define the bias.
ic = 0
nRegCode = length(uRegCode)
res  = data.frame(ic = c(1),i = c(1),j=c(1), c1 = c("A"),c2 = c("A"),a1=c(1),v1=c(1),a2=c(1),v2=c(1),v3=c(0),adiff = c(0),ldim = c(0))
D = matrix(0,ncol = nRegCode,nrow = 1)
Q = array(c(0))
V = array(c(0))
for (i in 1:(nRegCode-1)) {
  rci = as.character(uRegCode[i])
  localDatai = object[object[[regCode]] == rci,]
  ndat = dim(localDatai)[1]
  if (ndat > minKrige) imod = gdat[rci][[2]][[1]]
  for (j in (i+1):nRegCode) {
    rcj = as.character(uRegCode[j])
    boundaries = boundaryLines
    boundVec = (boundaries$c1 == rci & boundaries$c2 == rcj) |
                        (boundaries$c1 == rcj & boundaries$c2 == rci)
 #    cat(paste(i,j,rci,rcj,dim(lbound)[1],"\n"))
    if (sum(boundVec) > 0){
      lbound = boundaries[boundVec,]
#      cat(paste(rci,rcj,dim(lbound)[1],(dim(lbound)[1] > 0),"\n"))
      localDataj = object[object[[regCode]] == rcj,]
      mdat = dim(localDataj)[1]
      if (mdat > minKrige) {
        jmod = gdat[rcj][[2]][[1]]
        bord = paste(rci,rcj,sep="")
#        cat(paste(rci,rcj,dim(lbound)[1],dim(lbound),"\n"))
        xlinNew = SpatialLines(list(Lines(list(Line(coordinates(lbound))),ID = bord)))
        proj4string(xlinNew) = CRS(proj4string(boundaries))
        cat(paste("interpolating border between ",rci,rcj,"\n"))
        stri = krige(formulaString, localDatai,xlinNew,imod)
        strj = krige(formulaString, localDataj,xlinNew,jmod)
      }
      a1 = stri$var1.pred
      v1 = sqrt(stri$var1.var)
      a2 = strj$var1.pred
      v2 = sqrt(strj$var1.var)
      v3 = sqrt(stri$var1.var+strj$var1.var)
      adiff = abs(stri$var1.pred-strj$var1.pred)
      ic = ic+1
      if (ic > 1) D = rbind(D,0)
      D[ic,i] = 1
      D[ic,j] = -1
      Q[ic] = a1-a2
      V[ic] = v3
      res$ic = ic
      res$i = i
      res$j = j
      res$c1 = factor(res$c1,levels = c(levels(res$c1), rci))
      res$c2 = factor(res$c2,levels = c(levels(res$c2), rcj))
      res$c1 = rci
      res$c2 = rcj
      res$a1 = a1
      res$v1 = v1
      res$a2 = a2
      res$v2 = v2
      res$v3 = v3
      res$adiff = adiff
      res$ldim = dim(lbound)[1]
      if (ic == 1) {
        rest = res
      } else {
        rest = rbind(rest,res)
      }
    }
  }
}
#print(acdf)
return(list(D=D,Q=Q,V=V,rest=rest))
}

regionalVariograms = function(object,uRegCode,formulaString,minKrige = 5,regCode){
  nRegCode = length(uRegCode)
  ig = 0
  for (ic in 1:nRegCode) {
    rci = as.character(uRegCode[ic])
    localData = object[object[[regCode]] == rci,]
    ndat = dim(localData)[1]
    if (ndat > minKrige) {
      regionalVariogram = autofitVariogram(formulaString,localData, model = c("Gau","Exp","Sph"))$var_model
      if (ig == 0) {
        gdat = gstat(id = rci,formula = formulaString,model = regionalVariogram,data = localData)
        ig = 1
      } else {
        gdat = gstat(gdat,id = rci,formula = formulaString,model = regionalVariogram,data = localData)
      }
    }
  }
  return(gdat)
}
