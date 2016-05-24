.logit4Pred <-
function(lm.out, nm, mydata, my.formula, brief, res.rows,
         n.vars, n.pred, n.obs, n.keep, digits.d, pre, line,
         new.data, pred, pred.all, 
         numeric.all, in.data.frame, X1.new, 
         X2.new, X3.new, X4.new, X5.new, X6.new,
         pdf.file, pdf.width, pdf.height, ...) {

  pred.sort <- TRUE  # data must be sorted to find cases close to fitted=0.5

# ----------
# prediction
# ----------
     
  cat( "\n\n", "  FORECASTS", "\n\n")

  cat("Data, Fitted Values, Standard Errors\n")
  cat("   [sorted by fitted value]\n")
  if (n.keep > 50 && pred.all == FALSE && !new.data) 
    cat("   [to save space only some intervals printed, pred.all=TRUE to see all]\n")
  .dash(68)
  
  if (!new.data) {
    p.int <- data.frame(predict(lm.out, type="response", se.fit=TRUE))

    # classification
    prd <- integer(length=nrow(p.int))
    for (i in 1:nrow(p.int))
      if (p.int$fit[i]<0.5) prd[i] <- 0 else prd[i] <- 1

    if (all(prd == 0)  ||  all(prd==1)) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "All predicted values are ", prd[1], ".\n",
        "Something is wrong here.\n\n", sep="")
    }

    out <- cbind(lm.out$model[c(nm[seq(2,n.vars)],nm[1])],prd,p.int$fit,p.int$se.fit)
  }

  else {
    Xnew.val <- list(X1.new)
    if (n.vars > 2) for (i in 2:(n.pred)) {
      pp <- eval(parse(text=paste("X", toString(i),".new",sep="")))
      Xnew.val <- c(Xnew.val, list(pp))
    }
    Xnew <- expand.grid(Xnew.val)
    for (i in 1:(n.pred)) names(Xnew)[i] <- nm[i+1]

    p.int <- data.frame(predict(lm.out, type="response", se.fit=TRUE, newdata=Xnew))
    prd <- integer(length=nrow(p.int))
    for (i in 1:nrow(p.int))
      if (p.int$fit[i]<0.5) prd[i] <- 0 else prd[i] <- 1

    Ynew <- character(length = nrow(Xnew))
    Ynew <- ""
    out <- cbind(Xnew, Ynew, prd, p.int$fit,p.int$se.fit)
  }

  names(out)[n.vars+1] <- "predict"
  names(out)[n.vars+2] <- "fitted"
  names(out)[n.vars+3] <- "std.err"
  out <- data.frame(out)
  if (pred.sort) {
    o <- order(out[,n.vars+2])  # fitted value
    out <- out[o,]
  }


  if (n.keep < 25  || pred.all == TRUE || new.data)
    print(out, digits=digits.d)
  else {
    print(out[1:4,], digits=digits.d)
    cat("\n... for the rows of data where fitted is close to 0.5 ...\n\n")
    i.mid <- which.min(abs(0.5-sort(p.int$fit)))  # requires out is sorted by fit
    print(out[(i.mid-2):(i.mid+2),], digits=digits.d)
    cat("\n... for the last 4 rows of sorted data ...\n\n")
    print(out[(n.keep-3):n.keep,], digits=digits.d)
  }
  .dash(68)

  # classification table
  if (!new.data) {
    mytable <- table(out[,n.vars], out[,n.vars+1])
    hit0 <- mytable[1,1]; hit1 <- mytable[2,2]
    mis0 <- mytable[1,2]; mis1 <- mytable[2,1]
    per0 <- hit0 / (hit0 + mis0); per1 <- hit1 / (hit1 + mis1)
    hitT <- hit0 + hit1
    tot0 <- mytable[1,1] + mytable[1,2]
    tot1 <- mytable[2,1] + mytable[2,2]
    totG <- tot0 + tot1
    perT <- hitT / totG
    per0G <- tot0 / totG
    per1G <- tot1 / totG
    ln <- nchar(nm[1])
    cat("\n\nClassification Table for", nm[1], "\n\n")
    if (is.factor(lm.out$model[,nm[1]]))
      cat(" 0: ", levels(lm.out$model[,nm[1]])[1], "\n",
          " 1: ", levels(lm.out$model[,nm[1]])[2], "\n", sep="")
    cat(.fmtc(" ",ln+8), "Baseline         Predicted", "\n")
    .dash(51)
    cat(.fmtc(" ",ln+7), "Total  %Tot        0      1  %Correct", "\n")
    cat(.fmtc(" ",ln), "  0  ", .fmti(tot0,6), .fmt(100*per0G,1,5), " ",
        .fmti(hit0,6), .fmti(mis0,6), "   ",
        .fmt(100*per0,1), "\n")
    cat(.fmtc(nm[1],ln), "  1  ", .fmti(tot1,6), .fmt(100*per1G,1,5), " ",
        .fmti(mis1,6), .fmti(hit1,6), "   ",
        .fmt(100*per1,1), "\n")
    .dash(51)
    cat(.fmtc(" ",ln), "Total", .fmti(totG,6), .fmtc(" ",25),
        .fmt(100*perT,1), "\n")
    cat("\n")
  }


  # graphics
  if (pred && n.pred==1 && !is.factor(lm.out$model[,nm[2]]) && is.null(X1.new)) {

    .opendev(pdf.file, pdf.width, pdf.height)

    x.values <- lm.out$model[,nm[2]]
    if (!is.factor(lm.out$model[,nm[1]])) {
      y.values <- lm.out$model[,nm[1]] 
      y.label <- nm[1]
    }
    else {
      y.values <- as.numeric(lm.out$model[,nm[1]])
      min.y <- min(y.values)
      y.label <- paste(nm[1], " ",
         "(0=", levels(lm.out$model[,nm[1]])[1], ",",
         " 1=", levels(lm.out$model[,nm[1]])[2], ")", sep="")
      for (i in 1:length(y.values))
        if (y.values[i] == min.y) y.values[i] <- 0 else y.values[i] <- 1
    }

    # plot
    plot(x.values,y.values, type="n", axes=FALSE, ann=FALSE, ylim=c(-.10,1.10), ...)

    usr <- par("usr")
    col.bg <- getOption("col.bg")
    rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border="black")

    col.grid <- getOption("col.grid")
    abline(v=axTicks(1), col=col.grid, lwd=.5)
    abline(h=axTicks(2), col=col.grid, lwd=.5)

    .axes(NULL, NULL, axTicks(1), axTicks(2),
          par("usr")[1], par("usr")[3], cex.axis=.8, col.axis="gray30", ...)

    main.lab <- "Logistic Fit and Scatterplot"
    sub.lab <- NULL
    .axlabs(nm[2], y.label, main.lab, sub.lab, max.lbl.y=3,
            cex.lab=getOption("lab.size"), cex.main=1.0, ...) 

    col.fill <- getOption("col.fill.pt")
    col.stroke <- getOption("col.stroke.pt")
    points(x.values,y.values, pch=21, col=col.stroke, bg=col.fill, cex=0.8)
    lines(x.values, p.int$fit, col=col.stroke, lwd=2)

    if (!is.null(pdf.file)) {
      dev.off()
      .showfile(pdf.file, "fitted values and scatter plot")
    }
  }

  else {  # scatterplot matrix for multiple regression
    if (numeric.all && in.data.frame) {

      .opendev(pdf.file, pdf.width, pdf.height)

      panel2.smooth <- function (x, y, pch=par("pch"), cex=.9,
        col.pt=getOption("col.stroke.pt"), col.smooth=getOption("col.stroke.bar"),
        span=2/3, iter=3, ...) 
      {
          points(x, y, pch=pch, col=col.pt, cex=cex)
          ok <- is.finite(x) & is.finite(y)
          if (any(ok)) 
            lines(lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
      }
      pairs(lm.out$model[c(nm)], panel=panel2.smooth)

      if (!is.null(pdf.file)) {
        dev.off()
        .showfile(pdf.file, "scatter plot matrix")
      }
    }
    else {
      cat("\n\n>>> No scatterplot matrix reported because not all variables are ")
      if (!in.data.frame) cat("in the data frame.\n")
      if (!numeric.all) cat("numeric.\n")
    }
  }


}
