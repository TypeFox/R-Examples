
options(contrasts = c("contr.sum", "contr.poly"))

is.balanced <- function(design,slack=1) {
  Lvln=sapply(1:NCOL(design), function(i) length(levels(design[,i])))
  Lcounts=sapply(1:NCOL(design),function(i)
    min(unlist(apply(design,2,table)[i])))
  check.for.missing <- apply(as.matrix(apply(design,2,table)),2,length)
  
  minimums = round(NROW(design)/Lvln)-slack
  if( all(Lcounts >= minimums ))
    if(all(check.for.missing >= Lvln))
      ret.value=TRUE
    else
      ret.value=FALSE
  else
    ret.value=FALSE
  ret.value
}


is.singular <- function(design) {
  offset.column=rep(1,NROW(design))
  aug.matrix= cbind(offset.column,data.matrix(design))
  Info = t(aug.matrix) %*% aug.matrix
  if( abs(det(Info)) < 1e-20)
    return.val = TRUE
  else
    return.val = FALSE
  return.val
}


is.close.to.orthogonal <- function(design,tol=0.2) {

  x=suppressWarnings(cor(data.matrix(design)))
  x[is.na(x)] <- 1.0
  y=max(c(max(abs(x)[upper.tri(x)]),max(abs(x)[lower.tri(x)])))
  if (y < tol)
    return.value = TRUE
  else
    return.value = FALSE
  return.value
}



mc.good.designs = function(orig.set,
  cards=NULL,slack=1,tol=.2,no.replace=TRUE,
  size = 100,max.trials=1000000
  ) {
  set.seed(0)
  samps=NULL
  designs=NULL
  out=NULL
  n = max.trials
  if(is.null(cards)) {
    cards=nrow(orig.set)-3
  }
  if((no.replace==FALSE) && (10*size > choose(NROW(orig.set),cards)) ){
    warning("the number of designs required is more than one tenth\n  of the designs available \n  With no replace=FALSE there is a high chance of repeats")
  }
    
  if(no.replace){
    kk=NROW(orig.set)
    all.samps = combn(kk,cards)
    n = dim(all.samps)[2]
    samp.order=sample(n)
    if (size > n)
       warning("no.replace=TRUE can produce at most ", n," cards, however ",
                size, "were requested")
  }
  cindex=1
  for(k in 1:n) {
    if( (k %% 5000) == 0 ) print(paste(k,cindex-1))
    if(no.replace){
      samp=all.samps[,samp.order[k]]
    }
    else {
      samp= sample(1:NROW(orig.set),cards)
    }
    design=orig.set[samp,]
    if(!is.singular(design)&&
       is.close.to.orthogonal(design,tol) &&
       is.balanced(design,slack)){
        samps[[cindex]]=samp
        designs[[cindex]]=design
        cindex = cindex+1
        #print(paste(cindex-1,size))
        if(cindex>size) break
      }    
  }

  if (is.null(designs)){
    warning("No valid designs found after looking at ",
            max.trials," designs")
    NULL
  }
  else {

    if((cindex -1) < size)
      warning("Only ",cindex-1," designs found after looking at ",
              n," designs")
    out=NULL
    out$samps=samps
    out$designs=designs
    out
  }
}


M.Conjoint=function(despack,data,type="linear"){
  despack=mc.despack.linear.conjoint(despack,data)
  despack=mc.despack.linear.utils(despack)
  despack=mc.importances(despack)
  despack$utils=mc.mean.over.design.utils(despack$all.utils)
  despack$imps=mc.mean.over.design.imps(despack$all.imps)
  
  
  cat("Average values are output \nThey may or may not be meaningful\n\n")
  
  output.utils( mean.utils(despack$utils))
  cat("\nImportances\n")
  
  output.imps(rowMeans(despack$imps))
  invisible(despack)
}


mc.get.one.design = function(all.cards,cards,slack=1,tol=.2,max.tries=1000000) {
  set.seed(0)
  designs=NULL
  tries=0
  while(TRUE) {
    tries=tries+1
    if ((tries%%10000) == 0) print(paste(tries,"designs checked"))
    design=all.cards[sample(1:NROW(all.cards),cards),]
    if(!is.singular(design)&&
       is.close.to.orthogonal(design,tol) &&
       is.balanced(design,slack)){

      break
    }


    if(tries > max.tries){
      print("maximum attempts exceeded")
      design=NULL
      break
    }

    
  }
  design
}
do.linear.conjoint = function(design,data){
        options(contrasts = c("contr.sum", "contr.poly"))
        out=NULL
        tmp <- paste("design$",names(design),sep="",collapse="+")
        form = as.formula(paste("data","~",tmp,sep=""))
        fit=lm(form)
        out$coeffs = coefficients(fit)
        out
      }



mc.despack.linear.conjoint= function(despack,data=NULL) {

   if(is.null(data)) {
     data=despack$data
   }

   
   new.despack=despack
   new.despack$designs=NULL
   new.despack$coeffs=NULL
   index=0

   for(i in 1:length(despack$designs)){

     design = despack$designs[[i]]
     if(is.null(despack$samps))
       samp=as.numeric(row.names(design))
     else
       samp=despack$samps[i]
     tdata= as.matrix(data)[unlist(samp),]
     tdata=apply(apply(as.matrix(tdata),2,order),2,order)
     temp= do.linear.conjoint(design,tdata)
     if(!any(is.na(temp$coeffs))) {
       #only use the design if no coeffs are NA
       index=index+1
       new.despack$designs[[index]]=design
       new.despack$coeffs[[index]]=temp$coeffs
     }
   }
   new.despack
 }
       
     

mc.despack.linear.utils= function(despack) {

despack.all.utils = NULL

  for(i in 1:length(despack$designs)){
    design=despack$designs[[i]]
    coeffs=despack$coeffs[[i]]
    tutils=linear.utils(design,as.matrix(coeffs))
    despack$all.utils[[i]]=tutils
  }
  despack
}

linear.utils = function(design,coeffs){
  Lvln=sapply(1:NCOL(design), function(i) length(levels(design[,i])))
  utils.names=sapply(1:NCOL(design), function(i) levels(design[,i]))
  blank.utils=sapply(1:length(Lvln), function(i) rep(0,Lvln[i]))
  for (i in 1:length(Lvln)) {
    names(blank.utils[[i]])=utils.names[[i]]
  }

  blank.utils=c(list(0),blank.utils)
  names(blank.utils[[1]])= c("intercept")
  names(blank.utils) = c("",colnames(design))
  all.utils=NULL

  for(i in 1:NCOL(as.matrix(coeffs))) {
    index=2
    blank.utils[[1]][1] = as.matrix(coeffs)[1,i]
    for(j in 1:length(Lvln)) {
      total=0

      for(k in 1:(Lvln[j]-1)){
        blank.utils[[j+1]][k] = as.matrix(coeffs)[index,i]
        total=total + blank.utils[[j+1]][k]
        index=index+1
      }
      blank.utils[[j+1]][Lvln[j]] = -1*total
    }
    index=0
    
    all.utils[[i]]=blank.utils
  }
  
  all.utils
}

mc.importances= function(despack){

  despack$all.imps=NULL
  
  for(k in 1:length(despack$all.utils)){
    utils= despack$all.utils[[k]]
    imps = matrix(rep(0,(length(utils[[1]])-1)*length(utils)))
    dim(imps) = c(length(utils[[1]])-1,length(utils))
    
    for(i in 1:length(utils)){
      
      for(j in 1:length(utils[[1]])-1){
        imps[j,i]=max(utils[[i]][[j+1]]) - min(utils[[i]][[j+1]])
      }
    }
    
  
    imps=sweep(imps,2,colSums(imps),"/")
    rownames(imps) = colnames(despack$design[[1]])
    despack$all.imps[[k]]=imps
  }
  despack
}
  
    
    
    
    
output.utils= function(utils){


cat("Utilities\n")

temp=round(as.matrix(unlist(utils)),3)
colnames(temp)=""
print(temp)

}



mc.add.to.design.fast = function(all.possible, old.design, cards.to.add = 3,slack=1, tol=0.2) {
  if (!design.good(old.design,slack,tol)){
    stop("old.design not usable")
  }
  extra.lines = df.diff(all.possible,old.design)
  num.extra.lines=NROW(extra.lines)
  
  search.order= sample(num.extra.lines)
  
  
  
  lines.found = 0
  new.lines = NULL
  
  for(i in 1:num.extra.lines){
    next.line.found = FALSE
    for(j in 1:NROW(old.design)){
      
      test.design = rbind(old.design[-j,],extra.lines[search.order[i],])
      
      if(design.good(test.design, slack, tol)){
        next.line.found = TRUE
        lines.found= lines.found +1
        new.lines = rbind(new.lines, extra.lines[search.order[i],])
        break
      }
      
      
      
    }
    if(lines.found == cards.to.add) break
    if (next.line.found) next
  }
  if (lines.found < cards.to.add){
    warning("Not enough lines found")
  }
  else{
    new.design = rbind(old.design,new.lines)
    out=NULL
    out$base.design = old.design
    out$added = new.lines
    out$design = new.design[sample(nrow(new.design)),]
    row.names(out$design)=1:nrow(out$design)
    out
    
  }
      
      
    
}

design.good = function(design,slack=1, tol=0.2) {
  if(!is.singular(design)&&
       is.close.to.orthogonal(design,tol) &&
       is.balanced(design,slack)){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
    
}

vector.of.rows <- function(A){
  # turn A into a vector of colon separated rows
  # A can be a data.frame or matrix
  # used to perform set operations
  # on data frames
  f <- function(...) paste(..., sep=":")
  a <- do.call("f", as.data.frame(A))
  a
}
  
df.diff <- function(A,B) {
  a=vector.of.rows(A)
  b=vector.of.rows(B)
  A=A[match(setdiff(a,b),a),]
  A
}


mc.add.to.design= function(all.possible, old.design, 
                              cards.to.add = 3,slack=1, tol=0.2,max.trials=100,max.good.designs=100) {
  
  out=NULL
  
  if (!design.good(old.design,slack,tol)){
    stop("old.design not usable")
  }
  extra.lines = df.diff(all.possible,old.design)
  num.extra.lines=NROW(extra.lines)
  
  all.samps = combn(num.extra.lines,cards.to.add)
  samp.order=sample(NCOL(all.samps))
  
  test.design=NULL
  best.design=NULL
  best.length=0
  
  
  for(i in 1:NCOL(all.samps)){
    
      
      test.design = rbind(old.design,extra.lines[all.samps[,samp.order[i]],])
      
      if(design.good(test.design, slack, tol)){
        
        y=suppressWarnings(mc.good.designs(test.design, NROW(old.design),size=max.good.designs))
        
        new.length=length(y$designs)
        if(new.length>best.length){
          best.design=test.design
          best.length=new.length
        }
        
      }
      
      if(i > max.trials) break
      if(best.length == max.good.designs) break  #no test.design will do better
      
   }
  if (is.null(best.design)){
    print("No design found")
    NULL
  }
  else{
   print(paste("Number good designs",best.length))
   
   out$base.design = old.design
   out$added = df.diff(best.design,old.design)
   out$design = best.design[sample(nrow(best.design)),]
   row.names(out$design)=1:nrow(out$design)
   out
   
  }
  
}

mc.get.initial.design= function(full.design,cards=NULL,extra.cards=3, slack=1, tol = .2, max.trials=100) {
    out = NULL
    if(is.null(cards)) {
      
      min.cards=1+ sum(sapply(1:NCOL(full.design), function(i) length(levels(full.design[,i])))-1)
      cards=min.cards + 3
      
      
    }


    base.design = mc.get.one.design(full.design,cards,slack, tol, 10000)
    out = mc.add.to.design(full.design, base.design,extra.cards,slack,tol, max.trials)

    out
}
  

mean.utils.utility.i=function(utils,i){
  x=lapply(1:length(utils),function(k) utils[[k]][[i]])
  colMeans(do.call(rbind,x))
}

mean.utils.case.i=function(all.utils,i){
  z=lapply(1:length(all.utils), function(j) all.utils[[j]][[i]])
  lapply(1:length(z[[1]]), function(i) mean.utils.utility.i(z,i))
}

mc.mean.over.design.utils = function(all.utils){
  lapply(1:length(all.utils[[1]]), function(i) mean.utils.case.i(all.utils,i))
}

mean.utils = function(utils){
  if(typeof(utils[[1]]) == "list") {
    lapply(1:length(utils[[1]]), function(i) mean.utils.utility.i(utils,i))
  }
  else {
    utils
  }
}

mc.mean.over.design.imps= function(all.imps){
  Reduce('+',all.imps)/length(all.imps)
}


output.imps=function(imps){
  x=as.matrix(imps)
  colnames(x)=rep("",ncol(x))
  print(x)
}

