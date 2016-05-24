########## R-function: summary.spm ##########

# For summarising results of spm().

# Last changed: 21 JAN 2005 by MPW

summary.spm <- function(object,...)
{
   # Create `presence' indicators

   lin.present <- !is.null(object$info$lin)
   pen.present <- !is.null(object$info$pen)
   krige.present <- !is.null(object$info$krige)
   random.present <- !is.null(object$info$random)

   single.pen.present <- ((length(object$info$pen$name)==1)
                          &(!lin.present)&(!krige.present)
                          &(!random.present))

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
     
      num.pen <- ncol(as.matrix(object$info$pen$x))      
      names.vec <- NULL
      for (i in 1:num.pen)
         names.vec <- c(names.vec,paste("f(",object$info$pen$name[i],")",
                        sep="")) 

      nonlin.table <- data.frame(signif(dfs[-df.omit.inds],4),
		       spars,num.knots,row.names=names.vec)

      dimnames(nonlin.table)[[2]] <- c("df","spar","knots")

      cat("\n\n")
      cat("Summary for non-linear components:\n\n")
      print(signif(nonlin.table,4))
      cat("\n")
      cat("Note this includes 1 df for the intercept.\n\n")
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

   if (lin.present)
   {
      num.obs <- length(object$info$lin$x)

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

   nonlin.table <- NULL
   if (pen.present&(!single.pen.present))
   {
      spars <- unlist(object$info$pen$spar)
      num.knots <- unlist(lapply(object$info$pen$knots,length))

      num.pen <- ncol(as.matrix(object$info$pen$x))      
      names.vec <- NULL
      for (i in 1:num.pen)
         names.vec <- c(names.vec,paste("f(",object$info$pen$name[i],")",
                          sep="")) 

       nonlin.table <- rbind(nonlin.table,
                              data.frame(signif(df.pen,4),
                              spars,num.knots,row.names=names.vec))
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
      dimnames(nonlin.table)[[2]] <- c("df","spar","knots")

      cat("\n\n")
      cat("Summary for non-linear components:\n\n")
      print(signif(nonlin.table,4))
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
   invisible() 
}

######### End of summary.spm ##########










