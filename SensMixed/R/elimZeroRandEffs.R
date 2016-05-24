
#### check if the term has zero variance 
checkIsZeroVar <- function(model, rand.term)
{  
  vcr <- VarCorr(model)
  ind.term <- findIndTerm(vcr, getGrSlrand(rand.term))
  stddev <- attr(vcr[[ind.term]], "stddev")
  if(abs(stddev) < 1e-07)
    return(TRUE)  
  return(FALSE)
}




### get the random terms
getRandTerms <- function(fmodel)
{
  terms.fm <- attr(terms(fmodel),"term.labels")
  ind.rand.terms <- which(unlist(lapply(terms.fm, 
                                        function(x) 
                                          substring.location(x, "|")$first))!=0)
  return(unlist(lapply(terms.fm[ind.rand.terms], 
                       function(x) paste("(",x,")",sep=""))))
}

### get names of variables of the slope and group part of random term
getGrSlrand <- function(rand.term)
{
  # find the names of variables for slope part sl and for group part gr
  rand.term1 <- substring2(rand.term,2,nchar(rand.term)-1)
  splGrSl <- unlist(strsplit(rand.term1, "|", fixed = TRUE))
  gr <- substring2(splGrSl[2],2,nchar(splGrSl[2]))
  sl <- unlist(strsplit(splGrSl[1], "+", fixed = TRUE))
  for(i in 1:length(sl))
  {
    if(i==1)
      sl[i] <- substring2(sl[i],1,nchar(sl[i])-1)
    else
    {
      sl[i] <- substring2(sl[i],2,nchar(sl[i])-1)
    }
  }
  # change 1 to (Intercept) or eliminate 0 in slope part
  sl[which(sl=="1")] <- "(Intercept)"
  if(length(which(sl == "0"))!=0)
    sl <- sl[-which(sl == "0")]
  return(list(gr = gr, sl = sl))
}


### Find Index of the random term in VarCorr matrix
findIndTerm <- function(vcr, GrSl)
{
  indsGr <- which(names(vcr) == GrSl$gr)
  for(indGr in indsGr)
  {
    if(length(which((names(attr(vcr[[indGr]], "stddev")) == GrSl$sl) 
                    == FALSE)) == 0)
      return(indGr)
  }
}

### check if there are no random terms in the model
checkPresRandTerms <- function(mf.final)
{
  sub.loc.rand <- substring.location(paste(mf.final)[3], "|")
  if(sub.loc.rand$first==0 && sub.loc.rand$last==0)
    return(FALSE)
  return(TRUE)
}



#### eliminate components with zero variance
elimZeroVar <- function(model)
{
  fmodel <- formula(model)    
  rand.terms <- getRandTerms(fmodel)
  fm <- paste(fmodel)
  for(rand.term in rand.terms)
  {
    if(checkIsZeroVar(model, rand.term))
    {
      fm[3] <- paste(fm[3], "-", rand.term)
        warning(paste("\n Random term",rand.term, 
                      "was eliminated because of standard deviation being equal to 0 \n", 
                      sep=" "), call. = FALSE, immediate. = TRUE)
    }
  }
  mf.final <-  as.formula(paste(fm[2], fm[1], fm[3], sep=""))
  mf.final <- update.formula(mf.final,mf.final)
  is.present.rand <- checkPresRandTerms(mf.final)
  if(!is.present.rand)
    return(lm(mf.final, data = model.frame(model)))
  else
    return(lmer(mf.final, data = model.frame(model)))

}
