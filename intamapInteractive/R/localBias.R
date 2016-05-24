


removeLocalBias = function(object,localBias,gid = "group", formulaString=value~1,regCode = "regCode") {
#  if (is.list(localBias)) localBias = localBias$bias
  depVar=as.character(formulaString[[2]])
  for (i in 1:length(localBias)) {
    lBias = localBias[[i]]$bias
    rci = names(localBias)[i]
    for (j in 1:dim(lBias)[1]){
      if (rci == "single") {
        lNew = object[[depVar]][object[[gid]] == lBias$groups[j]] - lBias$wls[j]
        object[[depVar]][object[[gid]] == lBias$groups[j]] = lNew
      } else {
        lNew = object[[depVar]][object[[regCode]] == rci & 
             object[[gid]] == lBias$groups[j]] - lBias$wls[j]
        object[[depVar]][object[[regCode]] == rci & 
             object[[gid]] == lBias$groups[j]] = lNew
      } 
    }
  }
  return(object)
}



findLocalBias = function(object,gid = "group",formulaString=value~1,regCode="regCode",...) {
  if (regCode %in% names(object)) {
    uRegCode = unique(object[[regCode]])
    nRegCode = length(uRegCode)
  } else nRegCode = 1
  
  igid = which(names(object) == gid)
  localBias = list()
  
  for (ic in 1:nRegCode) {
    if (nRegCode > 1) {
      regCodei = as.character(uRegCode[ic])
      localData = object[object[[regCode]] == regCodei, ]
    } else {
      regCodei = "single"
      localData = object
    }
#      localData@data[,igid] = factor(localData@data[,igid])
      ngroups = length(unique(localData@data[,igid]))
    if (ngroups > 1) {
      localBias[[regCodei]] = localNetworkBias(localData,igid,formulaString,...)
    }
  }
  return(localBias)
}

localNetworkBias = function(localData,igid,formulaString,minKrige = 3,...)  {
#  groups = sort(unique(localData@data[,igid]))
  depVar = formulaString[[2]]
  groups = as.character(sort(unique(localData@data[,igid])))
  nGroups = length(groups)                                     
  D = matrix(0,ncol = nGroups,nrow = 1)
  varModel = vgm(0,"Nug",0,100)
  V = array()
  for (i in 1:nGroups) {
    ig = groups[i]
    groupData = localData[localData@data[,igid] == ig,]
    class(groupData) = class(localData)[1]
    ndat = dim(groupData)[1]
#    cat(paste(i,ig,"\n"))

    if (ndat > minKrige) varModel = autofitVariogram(formulaString,groupData,model="Sph")$var_model
    print(ndat)
    print(varModel)
    if (i == 1) {
      gdat = gstat(NULL,id = as.character(ig),formula = formulaString,model = varModel,data = groupData)
    } else {
      gdat = gstat(gdat,id = ig,formula = formulaString,model = varModel,data = groupData)
    }
  }
  im = 0
  Q = c(1)

  for (i in 1:(nGroups-1)) {
    ig = groups[i]
    ndat = sum(localData@data[,igid] ==ig)
    groupDatai = gdat[[1]][ig][[1]]$data
    if (ndat > minKrige) imod = gdat[[2]][ig][[1]]
    for (j in (i+1):nGroups) {
      jg = groups[j]
      mdat = sum(localData@data[,igid]==jg)
      groupDataj = gdat[[1]][jg][[1]]$data
      if (mdat > minKrige) jmod = gdat[[2]][jg][[1]]
      cArea = commonArea(groupDatai,groupDataj)
      bi = cArea[[1]]
      bj = cArea[[2]]

      if (bj > 0.5 & ndat > minKrige) {
        krigj <- krig(formulaString, groupDatai,groupDataj,imod)
        im = im + 1
        D[im,i] = 1
        D[im,j] = -1
        Q[im] = krigj[[1]]
        V[im] = krigj[[2]]

        kresj = data.frame(icol = i, jcol = j,krigj[[3]])
        if (im == 1 & FALSE) {
          write.csv(round(kresj),file = "sl_varios.txt",append=TRUE)
        } else if (FALSE) {
          write.csv(round(kresj,2),file = "sl_varios.txt",col.names=FALSE,append=TRUE)
        }
        D = rbind(D,matrix(0,ncol=nGroups,nrow=1))

      }
      if (bi > 0.5 & mdat > minKrige) {
        krigi <- krig(formulaString, groupDataj,groupDatai,jmod)
        im = im + 1
        D[im,j] = 1
        D[im,i] = -1
        Q[im] = krigi[[1]]
        V[im] = krigi[[2]]
        kresi = data.frame(icol = j, jcol = i,krigi[[3]])
        D = rbind(D,matrix(0,ncol=nGroups,nrow=1))
      }
    }
  }                                                                                                               
#  
  im = im+1
  D[im,] = 1
  Q[im] = 0
  V[im] = min(V)
  rDiff = list(D=D,Q=Q,V=V)
  locBias = dSolve(rDiff)  
  locBias = list(bias = data.frame(groups,ols = locBias$ols,ols.std = locBias$ols.std,
                       wls = locBias$wls,wls.std = locBias$wls.std),
                       D = D, V = V, Q = Q)
  return(locBias)

}





krig = function(form, groupDatai, groupDataj, imod) {
  depVar = as.character(form[[2]])
  jlen = dim(groupDataj)[1]
  strj = krige(form, groupDatai, groupDataj, imod)
  strj$obs = groupDataj[[depVar]]
  weight = 1/strj$var1.var
  aest = sum(weight*strj$var1.pred - weight*groupDataj[[depVar]])/sum(weight)
  # prediction variance - avar - is wrong
  errvar = 1/(jlen-2)*sum((groupDataj[[depVar]] - strj$var1.pred+aest)^2)
#  avar = (1/jlen + ((mean(strj$var1.pred))^2)/(jlen-1)/var(strj$var1.pred))*errvar
#  avar = var(strj$var1.pred-groupDataj[[depVar]])/jlen
  avar = ((weight%*%((strj$var1.pred-strj$obs)^2))/sum(weight)-aest^2)/jlen
  avar = (weight%*%((strj$var1.pred-strj$obs-aest)^2))/sum(weight)/jlen
#  print(asdf)
  return(list(aest,avar,strj)) 
}








