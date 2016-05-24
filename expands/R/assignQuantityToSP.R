assignQuantityToSP<-function (cbs, dm, colName = "PM_cnv", keepAmbigSeg = FALSE) 
{
  print("Assigning copy number to SPs...")
  SPs = sort(unique(c(dm[, "SP_cnv"],dm[,"SP"])))
  spSize = unique(round(SPs * 1000)/1000)
  x = colnames(cbs)
  x[(length(x) + 1):(length(x) + length(SPs))] = paste("SP_", 
                                                       as.character(spSize), sep = "")
  ploidy = matrix(cbind(cbs, matrix(NaN, nrow(cbs), length(SPs))), 
               nrow = nrow(cbs), ncol = length(x), dimnames = list(1:nrow(cbs), 
                                                                   x))
  toD = c()
  allAssigned=c();
  for (k in 1:nrow(ploidy)) {
    if (mod(k, 100) == 0) {
      print(paste("Finding overlaps for CBS segment", k, 
                  "out of ", nrow(ploidy), "..."))
    }
    idx = which(dm[, "chr"] == ploidy[k, "chr"] & dm[, "startpos"] >= 
                  ploidy[k, "startpos"] & dm[, "startpos"] <= ploidy[k, "endpos"])
    if (length(idx) == 0) {
      next
    }
    for (j in 1:length(idx)) {
      dmx = dm[idx[j], ]
      ##First assign SP with abberant copy number
      if (is.na(dmx["SP_cnv"]) || is.na(dmx[colName])) {
        next
      }
      sp = paste("SP_", as.character(round(dmx["SP_cnv"] * 
                                             1000)/1000), sep = "")
      thisAssigned=c(k,dmx["SP_cnv"] ,ploidy[k,"chr"],ploidy[k,"startpos"],ploidy[k,"endpos"],dmx[colName]);
      allAssigned = rbind(allAssigned, thisAssigned)
      if (is.na(ploidy[k, sp]) || ploidy[k, sp] == as.double(dmx[colName])) {
        ploidy[k, sp] = as.double(dmx[colName])
      } else {
        toD = rbind(toD, thisAssigned)
      }
      
      if (is.na(dmx["SP"]) || dmx["SP"]==dmx["SP_cnv"]) {
        next
      }
      ##Next assign SP with copy number, if co-occurence assumption violation applies 
      if (is.na(dmx["PM"]) ){
        next
      }
      sp = paste("SP_", as.character(round(dmx["SP"] *
                                             1000)/1000), sep = "")
      thisAssigned=c(k,dmx["SP"] ,ploidy[k,"chr"],ploidy[k,"startpos"],ploidy[k,"endpos"],dmx["PM"]);
      allAssigned = rbind(allAssigned, thisAssigned)
      if (is.na(ploidy[k, sp]) || ploidy[k, sp] == dmx["PM"]) {
        ploidy[k, sp] = as.double(dmx["PM"])
      } else {
        toD = rbind(toD, thisAssigned)
      }
      
    }
  }
  
  ##Either remove ambiguous segments or calculate their median ploidy based on SNV ploidy 
  printErr=FALSE;
  if (!is.null(nrow(toD))) {
    toD=matrix(toD,nrow=nrow(toD),ncol=ncol(toD),dimnames=list(
      rows=c(1:nrow(toD)),cols=c("Idx","SP","chr","startpos","endpos",colName)))
    colnames(allAssigned)=colnames(toD)  
    uD=unique(toD);
    for (i in 1:nrow(uD)){
      ii=which(as.numeric(uD[i,"SP"])==as.numeric(allAssigned[,"SP"]) & 
                 as.numeric(uD[i,"chr"])==as.numeric(allAssigned[,"chr"]) &
                 as.numeric(uD[i,"startpos"])==as.numeric(allAssigned[,"startpos"]) &
                 as.numeric(uD[i,"endpos"])==as.numeric(allAssigned[,"endpos"]));
      sp = paste("SP_", as.character(round(as.numeric(uD[i,"SP"] )* 1000)/1000), sep = "")
      tmpI=as.numeric(uD[i,"Idx"]);
      if (!keepAmbigSeg){
        ploidy[tmpI,sp]=NA;
        printErr=TRUE;    
      }else{
        ploidy[tmpI,sp]=round(median(as.numeric(allAssigned[ii,colName])));
      }
    }
  }else if (is.null(nrow(toD)) && length(toD) > 0) {
    if (!keepAmbigSeg){
      ploidy[as.numeric(toD[1]),]=NA;
      printErr=TRUE;    
    }
  }
  
  if(printErr){
    print(paste("Ambiguous SP specific ploidies found for ",length(toD)," segment-SP pairs."));
    print("Ploidies not assigned for these segments in corresponding SPs.")
  } 
  

  
  dm1=try(.assignPloidyToSPwithSNV(dm,ploidy),silent=FALSE)
  if(class(dm1)!="try-error"){
    dm=dm1;
  }
  out=list("dm"=dm, "ploidy"=ploidy);
  
  print("... Done.")
  if (keepAmbigSeg){
    print("Warning: parameter <keepAmbigSeg> set to TRUE. Output includes segment-assignements where subpopulation specific ploidy is ambiguous.Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements.")
  }
  
  return(out)
}

.assignPloidyToSPwithSNV <-function(dm,ploidy){
  ##PM may either be 2 or may be equal to the copy number of SP_cnv, if and only if this SP is a descendant of SP_cnv
  ii=which(!is.na(dm[,"SP"]) & (is.na(dm[,"SP_cnv"]) | dm[,"SP"]!=dm[,"SP_cnv"])); ## PM_cnv refers to ploidy of SP_cnv, not to ploidy of SP. Latter may be 2
  spIdx=grep("SP",colnames(ploidy))
  if (!isempty(ii)){
    for (k in ii){
      idx = which(ploidy[, "chr"] == dm[k, "chr"] & ploidy[, "startpos"] <= 
                    dm[k, "startpos"] & ploidy[, "endpos"] >= dm[k, "startpos"])
      if(!isempty(idx)){
        sp = paste("SP_", as.character(round(as.numeric(dm[k,"SP"] )* 1000)/1000), sep = "")
        idx=idx[which.min(ploidy[idx,"endpos"]-ploidy[idx,"startpos"])];##Smallest overlapping segment
        dm[k,"PM"]=max(dm[k,"PM_B"],ploidy[idx,sp]);
      }
    }
  }
  return(dm);
}
