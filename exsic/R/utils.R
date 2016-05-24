
update.pb <- function(pb, i, txt=""){
  if(!is.null(pb)){
    if(txt!="")   cat(paste("\n   Processing table: ",txt,"\n",sep=""))
    setTxtProgressBar(pb, i)
  }
}

#' Strips a text of any final dots.
#' 
#' This is a small helper function to prepare the final exsic table.
#' 
#' Intended for eliminating any duplicated dots in the final document due to final dots usually in the 
#' 'locnotes' (location notes) column.
#' 
#' @param txt a text string
#' @return a text string without final dot
#' @author Reinhard Simon
#' @aliases strip.final.dot
#' @export
#' @family helper
strip.final.dot <- function(txt=""){
  res = ""
  #if(missing("txt")) return("")
  if(is.null(txt)) return("")
  if(is.na(txt) ) return("")
  if(!is.character(txt) ) return("")
  
  if(is.character(txt)){
    txt = str_trim(txt)
    if(txt=="" | txt==".") return("")
    len = str_length(txt)
    lst = str_sub(txt,len,len+1)
    if(lst==".") {
      res = str_sub(txt,1,len-1)  
    } else {
      res = txt
    }
    
  }
  return(res)
  
}

add.dots <-function(s){
  x = toupper(s)
  x = str_trim(s)
  x = str_split(x,"")[[1]]
  x = x[-c(1)]
  x = paste(x,collapse=".")
  x = paste(x,".",sep="")
  str_trim(x)
}


#' Constructs a final citation of collector names.
#' 
#' Uses two of the BRAHMS fields (collector, addcoll). Expects names in both fields to conform to the
#' pattern: "Lastname, A" where A are the list of initials without dots; if no comma present the function
#' assummes all words to be part of the lastname like in: van de Bergh. 
#' 
#' The 'collector' field must have only one name; the 'addcoll' field may have several names separated by ; or
#' it may be empty. Like in the case of 'collector' a missing comma indicates that all words till the next ; 
#' are part of the last name.
#' 
#' Two collector names will be separated by '&'; more than two names will result in a citation with just
#' the (principal) collectors name followed by 'et al.'.
#' 
#' @note The function does not yet handle any other special cases like including an 'obligatory' second
#' author as permitted in some specimen citations when the author of the treatment is also a secondary
#' collector.
#' 
#' @note In any case, using this function or editing the final citation manually, the final citation should
#' be in a column called 'collcite'.
#' 
#' @param collector a name
#' @param addcoll a name or list of names; may be empty
#' @param initials use 'none' or 'before' to indicate if initials should be used for the citation
#' @param dots boolean; should dots be used to separate the initials?
#' @return a text
#' @author Reinhard Simon
#' @aliases coll.cite
#' @export
#' @family helper
coll.cite <-function(collector, 
                     addcoll="", 
                    initials=c("none","before"),  
                    dots = TRUE){
  res=""
  coll = collector
  coll.last = strsplit(coll,",")[[1]][1]
  if("before" == initials[1]){
    coll.init = strsplit(coll,",")[[1]][2]
    if(!is.na(coll.init)){
      res = coll.init
      if(dots){
        res = add.dots(res)
      }
      res = paste(res," ",sep="")
    }
  }
  res = paste(res,coll.last,sep="")
  addc.cnt = addcoll
  if(str_detect(addcoll,";")) {
    addc.cnt = str_split(addcoll,";")[[1]]
  } else {
    
  }
  if(length(addc.cnt)>0){
    n = length(addc.cnt)
    if(n>1){
      res=paste(res," et al.",sep="")
    }
    if(n==1){
      if(str_detect(addc.cnt[1],",")){
        add.ln = str_split(addc.cnt[1],",")[[1]][1]
        
        if("before" == initials[1]){
          add.init =  str_split(addc.cnt[1],",")[[1]][2]
          if(!is.na(add.init)){
            res.add = add.init
            if(dots){
              res.add = add.dots(res.add)
            }
          }
          add.ln = paste(res.add, add.ln)
          add.ln = str_trim(add.ln)
        }
        
        res = paste(res," & ", add.ln, sep="")
      } else {
        addcoll = str_trim(addcoll)
        if(addcoll!=""){
          res = paste(res," & ", addcoll, sep="")  
        }
        
      }
    }
  } else {
    #res = paste(res," & ", addcoll, sep="")
  }
  res = str_replace(res," & NA","")
  res = str_trim(res)
  res
}

# 
# grp = as.data.frame(cbind(
#   number = c(7,9,10,11,15,16,17,18,21,22),
#   species= c("3","1", "2", "2", "3", "3", "3","5","3","3")
#   ), stringsAsFactors=FALSE)


group.specimens <- function(astr, group.join = "-") {
  #print(astr)
  #print(group.join)
  ss = astr
  try({
  s = str_sub(ss,1, str_length(ss)-1)
  s1=str_trim(str_split(s,",")[[1]])
  # separate the collector number from the species and make a table
  s2=t(as.data.frame(sapply(str_trim(s1),str_split," ")))
  s2 = as.data.frame(s2, stringsAsFactors=FALSE)
  row.names(s2) = 1:nrow(s2)
  s2[,1] = as.character(s2[,1])
  s2[,2] = as.character(s2[,2])
  s2[s2[,1]=="s.n.",1]=""
  
    s2[,1] = as.integer(s2[,1])  
  
  
  names(s2) = c("number", "species")
  #print(s2)
  grp = s2
  n = nrow(grp) - 1
  if(n==0) return(astr)
  res=""
  grp$number = as.integer(grp$number)
  grp = cbind(grp, spr = rep(", ", nrow(grp)))
  grp$spr = as.character(grp$spr)
  grp$species = as.character(grp$species)
  grp = cbind(grp, grp = rep(1, nrow(grp)))
  
  
  # 1st summarize ids by species: 
  # 1b order within species & separate by comma
  grp = grp[with(grp, order(species, number)), ]
  
  
  # 2nd condense seq
  # for each species check which numbers are consecutive and then 
  
  sps = unique(grp$species)
  for(s in 1:length(sps)){
    sgrp = grp[grp$species == sps[s], ]
    gc = 0
    if(nrow(sgrp)>2){
      
      for(g in 2:nrow(sgrp)) {
        if(is.na(sgrp$number[g])) {
          sgrp$grp[g] = g
        } else {
          if((sgrp$number[g-1]+1) != sgrp$number[g]) {
            gc = gc+1
            sgrp$grp[g] = g
            
          } else {
            #sgrp$spr[g] = "-"
            sgrp$spr[g-1] = "-"
          }
          
        }
      }
      grp[grp$species == sps[s], ] = sgrp
    }
  }
  # keep only the first and last of each group bigger than 2!
  for(s in 1:nrow(grp)){
    mrk = grp$spr[s]
    if(mrk == "-"){
      if(s < (nrow(grp)-1)){
        mrk_1 = grp$spr[s+1]
        if(mrk_1 == "-") {
          grp$spr[s+1] = "x"
        }
      }
      
    }
  }
  grp$spr[nrow(grp)] = "."
  grp = grp[grp$spr!="x",]
  
  # Now build the string
  # before: handle NA -> s.n.
  grp$number = as.character(grp$number)
  if(any(is.na(grp$number))) {
    grp[is.na(grp$number),"number"] = "s.n."
  }
  
  res = ""
  for(s in 1:length(sps)){
    sgrp = grp[grp$species == sps[s], ]
    sgrp[nrow(sgrp), 3] = " "
    for(ss in 1:nrow(sgrp)){
      res = paste(res, sgrp$number[ss], sgrp$spr[ss], sep="")  
    }
    res = paste(res, sps[s], sep="" )
    if(s<length(sps)){
      res = paste(res,", ", sep="")
    } else {
      res = paste(res,".", sep="")
    }
  }
  astr = res
  
#   lwr = grp$number[1]
#   upr = grp$number[1]
#   for(i in 1:(nrow(grp)-1)){
#     this.species = grp$species[i]
#     next.species = grp$species[i + 1]
#     this.number  =  grp$number[i]
#     next.number  =  grp$number[i + 1]
#     if(!is.na(next.number) & !is.na(this.number)){
#       if((next.species == this.species & (next.number == (this.number+1))) ){
#         upr = next.number
#       } else {
#         if(lwr == upr){
#           res = paste(res,lwr," ",this.species,", ", sep="")
#         } else {
#           res = paste(res,lwr,group.join,upr," ",this.species,", ", sep="")  
#         }
#         
#         lwr = next.number
#         upr = lwr
#       }
#       
#     }
#     #print(paste(i, this.number, this.species, next.number, next.species, lwr, upr))    
#   }
#   if(lwr == upr){
#     res = paste(res,lwr," ",next.species,".", sep="")
#   } else {
#     res = paste(res, lwr, group.join, upr, " ", next.species,".", sep="")  
#   }
#     
#     astr = res
   },silent=TRUE)
   
  return(astr)
}
