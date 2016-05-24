########## R-function: summary.asp ##########

# For summarising results of asp2(). Based on summary.spm() of package SemiPar


summary.asp <- function(object,test1=FALSE,test2=FALSE,signif=0.05,...)
{
   # Create `presence' indicators

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)
   krige.present <- !is.null(object$info$krige)
   random.present <- !is.null(object$info$random)

   single.pen.present <- ((length(object$info$pen$name)==1)
                          &(!lin.present)&(!krige.present)
                          &(!random.present))

   nonlin.table=nonlin.table2=NULL

# scbTest
  if (pen.present & (test1|test2)){
    if (test1) {
      testobj1=scbM(object,level=1-signif)
      for (j in 1:length(testobj1$crit)){
        testobj1$tstat[[j]]= max(abs((testobj1$fitted[[j]]))*(testobj1$Stdev[[j]] )^(-1) )
        testobj1$pval[[j]] <- .C("stailp",crit=as.numeric(testobj1$tstat[[j]]), k0 = as.numeric(c(testobj1$k0[[j]],1)), d = as.integer(1), m = as.integer(2), rdf = as.numeric(testobj1$df), x = numeric(1), k = as.integer(1),PACKAGE="AdaptFitOS")$x
      }
    }
    if (test2){
      if (object$info$pen$basis!="os") {warning("Specification test (test2) only supported if B-splines were used for fitting the model."); test2=F}
      else testobj2=scbTest(object,level=1-signif)
    }
  }


   # First deal with constant only model.

   if ((!lin.present)&(!pen.present)&(!krige.present)&(!random.present))
   {
      num.obs <- length(object$info$y)

      fixed.coef <- object$fit$coef$fixed

      se.val <- sqrt(diag(object$aux$cov.mat))[1]

      const.table <- data.frame(signif(fixed.coef[1],5),
            	  	       signif(se.val,5),
	   		       row.names="intercept")

      const.table <- cbind(const.table,
                           signif(const.table[,1]/const.table[,2],5),
                           round((1-pt(abs(const.table[,1]/const.table[,2]),
                                   num.obs))*2,4))

      dimnames(const.table)[[2]] <- c("coef","se","ratio","p-value")

      cat("\n\n")
      cat("Summary for linear components:\n\n")
      print(signif(const.table,4))
      cat("\n\n")
   }

   # Next, deal with scatterplot smoothing models.

   if((!lin.present)&(!krige.present)&(single.pen.present))
   {
      dfs <- unlist(object$aux$df)+1
      spars <- unlist(object$info$pen$spar)

      num.knots <- unlist(lapply(object$info$pen$knots,length))

      df.omit.inds <- 1

      if (random.present)
          df.omit.inds <- c(df.omit.inds,length(dfs))
     
#      num.pen <- ncol(as.matrix(object$info$pen$x))
#      names.vec <- NULL
#      for (i in 1:num.pen)
#         names.vec <- c(names.vec,paste("f(",object$info$pen$name[i],")",
#                        sep=""))
#
#      nonlin.table <- data.frame(signif(dfs[-df.omit.inds],4),
#		       spars,num.knots,row.names=names.vec)
#      dimnames(nonlin.table)[[2]] <- c("df","spar","knots")

      names.vec <-pval.vec1<-tstat.vec1<- crit.vec1<-pval.vec2<-tstat.vec2<- crit.vec2<-basis.vec<-deg.vec<-pen.vec<- adap.vec<- NULL
#      for (i in 1:num.pen) {
        i=1
         names.vec <- c(names.vec,paste("f(",object$info$pen$name[i],")", sep=""))
         basis.vec= c(basis.vec,object$info$pen$basis)
         deg.vec = c(deg.vec,object$info$pen$degree[[i]][1])
         pen.vec = c(pen.vec,object$info$pen$degree[[i]][2])
         adap.vec = c(adap.vec,object$info$pen$adap[[i]])
         if (test1) {
           pval.vec1= c(pval.vec1,min(1,testobj1$pval[[i]]))
           tstat.vec1= c(tstat.vec1,testobj1$tstat[[i]])
           crit.vec1= c(crit.vec1,testobj1$crit[[i]])
         }
         if (test2) {
           pval.vec2= c(pval.vec2,min(1,testobj2$pval[[i]]))
           tstat.vec2= c(tstat.vec2,testobj2$tstat[[i]])
           crit.vec2= c(crit.vec2,testobj2$crit[[i]])
         }
#      }
      if (test1) nonlin.table <- data.frame(basis.vec,deg.vec ,pen.vec,adap.vec,signif(dfs[-df.omit.inds],3),
                            signif(spars,3),num.knots,"  |",round(tstat.vec1,3),round(crit.vec1,3),round(pval.vec1,3),row.names=names.vec)
      else   nonlin.table <- data.frame(basis.vec,deg.vec,pen.vec,adap.vec ,signif(dfs[-df.omit.inds],3),
                            signif(spars,3),num.knots,row.names=names.vec)

     if (test1) dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots","|","tstat",paste("crit(",signif,"%)",sep=""),"pval")
     else dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots")

      cat("\n\n")
      cat("Summary for non-linear components:\n\n")
#      print(signif(nonlin.table,4))
      print(nonlin.table)
      cat("\n")
      cat("Note this includes 1 df for the intercept.\n\n")
      if (test2){
        cat("\nTest for a polynomial of degree...:\n\n")
        nonlin.table2 <- data.frame(pen.vec-1,adap.vec,round(tstat.vec2,3),round(crit.vec2,3),round(pval.vec2,3),row.names=names.vec)
        dimnames(nonlin.table2)[[2]] <- c("degree","adap","tstat",paste("crit(",signif,"%)",sep=""),"pval")
        print(nonlin.table2)
      } else nonlin.table2=NULL
      cat("\n")
   }

   if (lin.present) num.lin <- ncol(as.matrix(object$info$lin$x)) 
   if (!lin.present) num.lin <- 0

   if (pen.present) num.pen <- ncol(as.matrix(object$info$pen$x)) 
   if (!pen.present) num.pen <- 0

   if (krige.present) num.krige <- 1
   if (!krige.present) num.krige <- 0
   
   if (random.present) num.random <- 1
   if (!random.present) num.random <- 0

   if ((pen.present)|(krige.present)|(random.present))
   {
      # Determine df vectors for each type
   
      dfs <- object$aux$df
   
      if (pen.present)
         df.pen <- dfs[(1+num.lin+1):(1+num.lin+num.pen)]
   
      if (krige.present)
         df.krige <- dfs[(1+num.lin+num.pen+1):
                         (1+num.lin+num.pen+num.krige)]
   
      if (random.present)
         df.random <- dfs[(1+num.lin+num.pen+num.krige+1):
                             (1+num.lin+num.pen+num.krige+num.random)]
   }

   # Create linear components table
   lin.table=NULL
   if (lin.present)
   {
      num.obs <- nrow(object$info$lin$x)

      fixed.coefs <- object$fit$coef$fixed

      se.vals <- sqrt(diag(object$aux$cov.mat))[1:(num.lin+1)]

      lin.table <- data.frame(signif(fixed.coefs[1:(num.lin+1)],5),
            	  	       signif(se.vals,5),
	   		       row.names=c("intercept",object$info$lin$name))

      lin.table <- cbind(lin.table,signif(lin.table[,1]/lin.table[,2],5),
                           round((1-pt(abs(lin.table[,1]/lin.table[,2]),
                                   num.obs))*2,4))

      dimnames(lin.table)[[2]] <- c("coef","se","ratio","p-value")

      cat("\n\n")
      cat("Summary for linear components:\n\n")
      print(signif(lin.table,4))
      cat("\n\n")
   }

   # Create table for penalised and kriging terms

   if (pen.present&(!single.pen.present)){
      nonlin.table <- NULL
      spars <- unlist(object$info$pen$spar)
      num.knots <- unlist(lapply(object$info$pen$knots,length))

      num.pen <- ncol(as.matrix(object$info$pen$x))
      names.vec <-pval.vec1<-tstat.vec1<- crit.vec1<-pval.vec2<-tstat.vec2<- crit.vec2<-basis.vec<-deg.vec<-pen.vec<- adap.vec<- NULL
      for (i in 1:num.pen) {
         names.vec <- c(names.vec,paste("f(",object$info$pen$name[i],")", sep=""))
         basis.vec= c(basis.vec,object$info$pen$basis)
         deg.vec = c(deg.vec,object$info$pen$degree[[i]][1])
         pen.vec = c(pen.vec,object$info$pen$degree[[i]][2])
         adap.vec = c(adap.vec,object$info$pen$adap[[i]])
          if (test1) {
           pval.vec1= c(pval.vec1,min(1,testobj1$pval[[i]]))
           tstat.vec1= c(tstat.vec1,testobj1$tstat[[i]])
           crit.vec1= c(crit.vec1,testobj1$crit[[i]])
         }
         if (test2) {
           pval.vec2= c(pval.vec2,min(1,testobj2$pval[[i]]))
           tstat.vec2= c(tstat.vec2,testobj2$tstat[[i]])
           crit.vec2= c(crit.vec2,testobj2$crit[[i]])
         }
      }
#      if (test) nonlin.table <- rbind(nonlin.table,
#                            data.frame(basis.vec,deg.vec ,pen.vec,adap.vec,signif(df.pen,3),
#                            signif(spars,3),num.knots,"  |",round(tstat.vec,3),round(crit.vec,3),round(pval.vec,3),row.names=names.vec))
#      else   nonlin.table <- rbind(nonlin.table,
#                            data.frame(basis.vec,deg.vec,pen.vec,adap.vec ,signif(df.pen,3),
#                            signif(spars,3),num.knots,row.names=names.vec))

      if (test1) nonlin.table <- data.frame(basis.vec,deg.vec ,pen.vec,adap.vec,signif(df.pen,3),
                            signif(spars,3),num.knots,"  |",round(tstat.vec1,3),round(crit.vec1,3),round(pval.vec1,3),row.names=names.vec)
      else   nonlin.table <- data.frame(basis.vec,deg.vec,pen.vec,adap.vec ,signif(df.pen,3),
                            signif(spars,3),num.knots,row.names=names.vec)


   }

   if (krige.present&(!single.pen.present))
   {
      spar.krige <- object$info$krige$spar
 
      num.knots.krige <- nrow(object$info$krige$knots)

      name.vec <-  paste("f(",object$info$krige$name[1],",",
			         object$info$krige$name[2],")",sep="")

      nonlin.table <-  rbind(nonlin.table,
                                c(signif(df.krige,4),
                                  spar.krige,num.knots.krige))

      nonlin.table <- as.data.frame(nonlin.table)

      dimnames(nonlin.table)[[1]][nrow(nonlin.table)] <- name.vec
   }

   if (((pen.present)|(krige.present))&(!single.pen.present))
   {
#      if (test) dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots","|","tstat","crit","pval")
#      else dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots")

      if (test1) dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots","|","tstat",paste("crit(",signif,"%)",sep=""),"pval")
      else dimnames(nonlin.table)[[2]] <- c("basis","deg","pen","adap","df","spar","knots")

      cat("\n\n")
      cat("Summary for non-linear components:\n\n")
     # print(signif(nonlin.table,4))
      print(nonlin.table)
      if (test2){
        cat("\nTest for a polynomial of degree...:\n\n")
        nonlin.table2 <- data.frame(pen.vec-1,adap.vec,round(tstat.vec2,3),round(crit.vec2,3),round(pval.vec2,3),row.names=names.vec)
        dimnames(nonlin.table2)[[2]] <- c("degree","adap","tstat",paste("crit(",signif,"%)",sep=""),"pval")
        print(nonlin.table2)
      }  else nonlin.table2=NULL
      cat("\n\n")
   }

   if (random.present)
   {
      random.table <- as.data.frame(signif(df.random,4))

      dimnames(random.table)[[1]] <- "random intercept"
 
      dimnames(random.table)[[2]] <- "df"

      cat("Summary for random intercept component:\n\n")
      print(signif(random.table,4))
      cat("\n\n")
   }
   invisible(list(lin=lin.table,pen=nonlin.table,test2=nonlin.table2))
}

######### End of summary.spm ##########










