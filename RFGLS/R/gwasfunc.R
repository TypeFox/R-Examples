#Function for use in doing single-SNP analyses as in GWAS.
#Adapted by Rob Kirkpatrick (May 2013) from syntax by Xiang Li and Saonli Basu
######################################################################################
bare_fgls <- function(Y,X,vmat,gkmat=NULL){
  if(is.null(gkmat)){
    list.vmat<-listbdsmatrix(vmat,diag=T,id=F)
    vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)
    vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
    vmat.Inv<-forceSymmetric(vmat.Inv)
    gkmat<-as(chol(vmat.Inv),"sparseMatrix")
  }
  
  newz <-gkmat%*%X
  newy <- gkmat%*%Y
  lfit <- lm(newy[,1]~0+as.matrix(newz))
  ls <- summary(lfit)
  fit <- list(
    ctable = ls$coefficients,    
    df.residual = lfit$df.residual
  )
  return(fit)
}

# gwasfunc <- function(snp, pheno, id, x.covar.full, vmat, tlist, sizelist, med){
#   out <- NULL
#   if(is.null(x.covar.full)){test2.dat <- as.data.frame(na.omit(cbind(snp,pheno,id)))}
#   else{
#     test2.dat <- as.data.frame(na.omit(cbind(snp,pheno,id,x.covar.full)))
#     #covarcol <- 4+ncol(x.covar.full)-1
#     x.covar <- as.matrix(test2.dat[,-c(1:3)])
#   }
#   freqtabl <- table(test2.dat[,1]) #Check if SNP is monomorphic.
#   if(length(freqtabl)==1){
#     out <- c("snp",rep(NA,5))
#     return(out)
#   }
#   else{
#     if(nrow(test2.dat)!=length(id)){ #Check if anyone is missing the current SNP; if so, cut them out of res cov matrix.
#       vmat0 <- vmat[which(vmat@Dimnames[[1]] %in% test2.dat[,3]),which(vmat@Dimnames[[2]] %in% test2.dat[,3])]
#     }
#     else{vmat0 <- vmat}
#   }
#   id <- test2.dat[,3]
#   if(is.null(x.covar.full)){
#     lme.out <- try(fgls(test2.dat[,2]~test2.dat[,1],#data=test2.dat,
#                         vmat=vmat0,tlist=tlist,sizelist=sizelist,med=med))
#   }
#   else{
#     lme.out <- try(fgls(test2.dat[,2]~test2.dat[,1]+x.covar,
#                         #data=test2.dat,
#                         vmat=vmat0,tlist=tlist, sizelist=sizelist, med=med))
#   }
#   if(class(lme.out)!="try-error"){
#     tmp <- c(lme.out$ctable[2,1],lme.out$ctable[2,2],lme.out$ctable[2,3],lme.out$df.residual,lme.out$ctable[2,4])
#     out <- c("snp",tmp)
#   }
#   else{out <- c("snp",rep(NA,5))} #If fgls() fails for some reason.
#   return(out)
# }
