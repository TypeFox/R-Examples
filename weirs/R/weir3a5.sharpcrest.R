"weir3a5.sharpcrest" <-
function(h, ht=NULL, b=NULL, B=NULL, P=NULL, L=NULL,
         r=0, A=NULL, alpha=1,
         slopeus="vertical",
         kc=NULL, kt=NULL, C=NULL,
         contractratio=NULL,
         extended=TRUE,
         header="", resetkts=TRUE,
         flowdigits=2, coedigits=3,
         verbose=FALSE, eps=0.001, maxit=20) {


  if(slopeus != "vertical") {
    if(as.logical(length(grep(":",slopeus)))) {
      slopeus <- as.numeric(unlist(strsplit(slopeus,":")));
      slopeus <- slopeus[1]/slopeus[2];
    } else {
      stop("slopeus does not contain a colon (:)");
    }
  } else {
    slopeus <- 0;
  }

  if(slopeus < 0) stop("upstream slope can not be negative");

  if(is.null(h)) {
     stop("head is NULL");
  }
  if(is.null(ht)) {
    ht <- rep(0, length(h));
  } else if(length(ht) == 1) {
    ht <- rep(ht, length(h));
  }
  if(length(h) != length(ht)) {
     stop("length of head vector is not equal to length of tail water head vector");
  }
  if(is.null(L)) {
     stop("weir length (along flow) is NULL");
  }
  if(is.null(B) | B <= 0) {
     stop("channel width is NULL or <= 0");
  }
  if(is.null(P) | P <= 0) {
     stop("weir height is NULL or <= 0");
  }
  if(! is.null(contractratio) & length(contractratio) != length(h)) {
     stop("user provided contraction ratio vector is not equal to length of head vector");
  }
  if(! is.null(A) & length(A) != length(h)) {
     stop("user provided approach area vector is not equal to length of head vector");
  }
  if(B < b) {
     stop("user provided channel width and weir width are incompatible");
  }
  if(length(L) == 1) L <- rep(L, length(h));
  if(length(L) != length(h)) {
     stop("length of weir-length vector is not equal to length of head vector");
  }
  if(! is.null(r) & (r < 0 | r > b)) {
     stop("implausible radius of curvature on abutment");
  }

  if(length(alpha) == 1) alpha <- rep(alpha, length(h));
  if(length(alpha) != length(h)) {
     stop("alpha vector is not equal to length of head vector");
  }

  if(! is.null(kc)) {
    if(length(kc) == 1) kc <- rep(kc, length(h));
    if(length(kc) != length(h)) {
       stop("contraction coefficient is not equal to length of head vector");
    }
  }

  if(! is.null(kt)) {
    if(length(kt) == 1) kt <- rep(kt, length(h));
    if(length(kt) != length(h)) {
      stop("submergence coefficient is not equal to length of head vector");
    }
  }

  if(! is.null(C)) {
    if(length(C) == 1) C <- rep(C, length(h));
    if(length(C) != length(h)) {
      stop("discharge coefficient is not equal to length of head vector");
    }
  }

  Qs <- Qo <- Qerr <- vector(mode="numeric", length=length(h));
  messages <- Cs <- kcs <- Vels <- bBs <- htHs <- Qs;
  kts <- rep(NA, length(h));
  g <- 32.2; g2 <- 2*g;

  for(i in 1:length(h)) {
    it <- 0;
    messages[i] <- "ok"; # assumed for initialization
  	hh <- h[i];

  	r.over.b  <- r/b;
    h.over.L  <- hh/L[i];
    h.over.P  <- hh/P;

  	if(h.over.P > 5) { # TWRI3(A5),p.4
  	  Qo[i]   <- NA;
      Qs[i]   <- NA;
      Qerr[i] <- NA;
      Vels[i] <- NA;
      Cs[i]   <- NA;
      kcs[i]  <- NA;
      kts[i]  <- NA;
      messages[i] <- "too much head";
  	  next;
  	}
  	if(hh == 0) {
  	  Qo[i]   <- 0.00;
      Qs[i]   <- 0.00;
      Qerr[i] <- 0.00;
      Vels[i] <- 0.00;
      Cs[i]   <- NA;
      kcs[i]  <- NA;
      kts[i]  <- NA;
      messages[i] <- "head zero";
  	  next;
  	}

    # Now that h.over.L and h.over.P are available---determine weir flow type
    if(slopeus == 0) {
      fig6 <- get("0.0000", .weir.nomographs$fig6);
      h.over.L.critical <- approx(fig6$h.over.P, fig6$h.over.L, h.over.P, rule=2)$y;
    } else if(slopeus == 1) {
      fig6 <- get("1.0000", .weir.nomographs$fig6);
      h.over.L.critical <- approx(fig6$h.over.P, fig6$h.over.L, h.over.P, rule=2)$y;
    } else {
      slopes6 <- as.numeric(ls(.weir.nomographs$fig6));
      tmp <- slopes6[slopes6 <= slopeus];
      if(length(tmp) != 0) slopes6.min <- max(tmp);
      tmp <- slopes6[slopes6 >= slopeus];
      if(length(tmp) != 0) slopes6.max <- min(tmp);
      #if(verbose) message("slopes6.min=",slopes6.min);
      #if(verbose) message("slopes6.max=",slopes6.max);
      if(slopeus < 1) {
        if(slopes6.min == slopes6.max) {
          fig6   <- get(sprintf("%.4f",slopes6.min), .weir.nomographs$fig6);
          h.over.L.critical <- approx(fig6$h.over.P, fig6$h.over.L, h.over.P, rule=2)$y;
        } else {
          fig6   <- get(sprintf("%.4f",slopes6.min), .weir.nomographs$fig6);
          h.over.L.min <- approx(fig6$h.over.P, fig6$h.over.L, h.over.P, rule=2)$y;
          fig6   <- get(sprintf("%.4f",slopes6.max), .weir.nomographs$fig6);
          h.over.L.max <- approx(fig6$h.over.P, fig6$h.over.L, h.over.P, rule=2)$y;
          h.over.L.critical <- approx(c(slopes6.min,slopes6.max),
                                      c(h.over.L.min, h.over.L.max), slopeus)$y;
        }
      } else { #  (slopeus > 1)
        h.over.L.critical <- 2.4; # TWRI3A5p.8 (a stretch in the assumption);
        # h/L = 2.4 is large and even though
      }
    }

    if(h.over.L < h.over.L.critical) {
       Qo[i]   <- NA;
       Qs[i]   <- NA;
       Qerr[i] <- NA;
       Vels[i] <- NA;
       Cs[i]   <- NA;
       kcs[i]  <- NA;
       kts[i]  <- NA;
       messages[i] <- "broad-crested";
       next;
    }

    if(is.null(kc)) {
      # Now that we know sharp-crested flow---determine kc
      if(is.null(contractratio)) {
         b.over.B <- b/B;
      } else {
         b.over.B <- contractratio[i];
      }
      bBs[i] <- b.over.B;

      b.over.B.ratios <- as.numeric(ls(.weir.nomographs$fig3));

      tmp <- b.over.B.ratios[b.over.B.ratios < b.over.B];
      if(length(tmp) != 0) b.over.B.ratios.min <- max(tmp);
      tmp <- b.over.B.ratios[b.over.B.ratios > b.over.B];
      if(length(tmp) != 0) b.over.B.ratios.max <- min(tmp);

      if(b.over.B < 0.90 & b.over.B >= 0.20) {
        if(b.over.B == 0.20) b.over.B.ratios.min <- 0.20;
        fig3   <- get(sprintf("%.2f",b.over.B.ratios.min),
                                     .weir.nomographs$fig3);
        kc.min <- approx(fig3$h.over.P, fig3$kc, h.over.P, rule=2)$y;
        fig3   <- get(sprintf("%.2f",b.over.B.ratios.max),
                                     .weir.nomographs$fig3);
        kc.max <- approx(fig3$h.over.P, fig3$kc, h.over.P, rule=2)$y;
        the.kc     <- approx(c(b.over.B.ratios.min,b.over.B.ratios.max),
                         c(kc.min, kc.max), b.over.B)$y;
      } else if(b.over.B <= 0.20) {
        Qo[i]   <- NA;
        Qs[i]   <- NA;
        Qerr[i] <- NA;
        Vels[i] <- NA;
        Cs[i]   <- NA;
        kcs[i]  <- NA;
        kts[i]  <- NA;
        messages[i] <- "too much contraction";
        next;
      } else {
        fig3   <- get(sprintf("%.2f",b.over.B.ratios.min),
                                     .weir.nomographs$fig3);
        kc.min <- approx(fig3$h.over.P, fig3$kc, h.over.P, rule=2)$y;
        the.kc     <- approx(c(b.over.B.ratios.min, 1),
                             c(kc.min, 1), b.over.B)$y;
      }

      if(r.over.b > 0.12) {
        the.kc <- 1.00; # TWRI3A5p.5
      } else if(r.over.b > 0) {
        the.kc <- approx(c(0,0.12), c(the.kc,1), r.over.b, rule=2)$y;
      }
    } else {
      the.kc <- kc[i];
    }


    if(is.null(C)) {
      # Now ready to determine C
      if(slopeus == 0) {
        fig2 <- get("0.0000", .weir.nomographs$fig2);
        the.C <- approx(fig2$h.over.P, fig2$C, h.over.P, rule=2)$y;
      } else if(slopeus == 1) {
        fig2 <- get("1.0000", .weir.nomographs$fig2);
        the.C <- approx(fig2$h.over.P, fig2$C, h.over.P, rule=2)$y;
      } else {
        slopes2 <- as.numeric(ls(.weir.nomographs$fig2));
        #if(verbose) message("Cval: slopes2=",slopes2);
        tmp <- slopes2[slopes2 <= slopeus];
        if(length(tmp) != 0) slopes2.min <- max(tmp);
        tmp <- slopes2[slopes2 >= slopeus];
        if(length(tmp) != 0) slopes2.max <- min(tmp);
        #if(verbose) message("Cval: slopes2.min=",slopes2.min);
        #if(verbose) message("Cval: slopes2.max=",slopes2.max);
        if(slopeus < 2) {
          if(slopes2.min == slopes2.max) {
            fig2  <- get(sprintf("%.4f",slopes2.min), .weir.nomographs$fig2);
            the.C <- approx(fig2$h.over.P, fig2$C, h.over.P, rule=2)$y;
          } else {
            fig2  <- get(sprintf("%.4f",slopes2.min), .weir.nomographs$fig2);
            C.min <- approx(fig2$h.over.P, fig2$C, h.over.P, rule=2)$y;
            #if(verbose) message("Cval: C.min=",C.min);
            fig2  <- get(sprintf("%.4f",slopes2.max), .weir.nomographs$fig2);
            C.max <- approx(fig2$h.over.P, fig2$C, h.over.P, rule=2)$y;
            #if(verbose) message("Cval: C.max=",C.max);
            the.C     <- approx(c(slopes2.min, slopes2.max),
                            c(C.min, C.max), slopeus)$y;
          }
        } else { #  (slopeus > 1)
          Qo[i]   <- NA;
          Qs[i]   <- NA;
          Qerr[i] <- NA;
          Vels[i] <- NA;
          Cs[i]   <- NA;
          kcs[i]  <- NA;
          kts[i]  <- NA;
          messages[i] <- "slopeus too shallow to determine C";
          next;
        }
      }
    } else {
      the.C <- C[i];
    }

    Cs[i]  <- the.C;
    kcs[i] <- the.kc;
    bigOval <- the.kc*the.C*b;

    Qint  <- bigOval*the.C*hh^1.5;
    Qo[i] <- Qold <- Qint;
    the.alpha <- alpha[i];
    the.A     <- ifelse(is.null(A), (hh+P)*B, A[i]);
    "afunc" <- function(Q) {
       vel     <- Q/the.A;
       vhead   <- vel^2/g2;
       Vels[i] <- vhead;
       H       <- hh + the.alpha * vhead;
       Qtmp    <- bigOval*H^1.5;
       return(Qtmp);
    }
    while(1) {
      it <- it + 1;
      Qnew <- afunc(Qold);
      if(Qnew == Inf) {
        Qo[i]  <- NA;
        Qs[i]  <- NA;
        kts[i] <- NA;
        messages[i] <- "nonconvergence";
        break;
      }
      err <- abs(Qnew - Qold);
      Qerr[i] <- err;
      if(it > maxit || err < eps) break;
      Qold <- Qnew;
    }
    Qs[i]     <- Qnew;
    vhead     <- (Qnew/the.A)^2/g2;
    H         <- hh + vhead;
    Vels[i]   <- vhead;
    ht.over.H <- ht[i]/H;
    htHs[i]   <- ht.over.H;
    H.over.P  <- H/P;
    #message("ht ", ht[i]);
    #message("ht.over.H ", ht.over.H);
    #message("H.overP ", H.over.P);
    
    if(is.na(kts[i]) & messages[i] != "nonconvergence") {
      if(ht.over.H > 0.95) {
        kts[i] <- NA;
        messages[i] <- "too much submergence to est. kt";
        next;
      } else if(ht.over.H == 0) {
        kts[i] <- 1;
        next;
      } else {
        if(H.over.P < 0.20) {
          kts[i] <- NA;
          messages[i] <- "H.over.P too small to est. kt";
          next;
        } else if(H.over.P > 2) {
          kts[i] <- NA;
          messages[i] <- "H.over.P too large to est. kt";
          next;
        } else {
          if(H.over.P == 0.20) {
            fig4 <- get("0.20", .weir.nomographs$fig4);
            the.kt <- approx(fig4$ht.over.H, fig4$kt, ht.over.H, rule=2)$y;
          } else if(slopeus == 2) {
            fig4 <- get("2.00", .weir.nomographs$fig4);
            the.kt <- approx(fig4$ht.over.H, fig4$kt, ht.over.H, rule=2)$y;
          } else {
            slopes4 <- as.numeric(ls(.weir.nomographs$fig4));
            #if(verbose) message("ktval: slopes4=",slopes4);
            tmp <- slopes4[slopes4 <= H.over.P];
            if(length(tmp) != 0) slopes4.min <- max(tmp);
            tmp <- slopes4[slopes4 >= H.over.P];
            if(length(tmp) != 0) slopes4.max <- min(tmp);
            #if(verbose) message("ktval: slopes4.min=",slopes4.min);
            #if(verbose) message("ktval: slopes4.max=",slopes4.max);
            fig4  <- get(sprintf("%.2f",slopes4.min), .weir.nomographs$fig4);
            kt.min <- approx(fig4$ht.over.H, fig4$kt, ht.over.H, rule=2)$y;
            #if(verbose) message("ktval: kt.min=",kt.min);
            fig4  <- get(sprintf("%.2f",slopes4.max), .weir.nomographs$fig4);
            kt.max <- approx(fig4$ht.over.H, fig4$kt, ht.over.H, rule=2)$y;
            #if(verbose) message("ktval: kt.max=",kt.max);
            the.kt <- approx(c(slopes4.min, slopes4.max),
                             c(kt.min, kt.max), H.over.P)$y;
            #if(verbose) message("kt: ",the.kt);
          }
          kts[i] <- the.kt;
        }
      }
    }
  }

  if(resetkts) kts[kts > 1] <- 1;
  
  if(extended) {
    z <- data.frame(head=h,
                    flow=round(Qs*kts,    digits=flowdigits),
                    delta=c(NA,diff(Qs*kts)),
                    flowfree=round(Qs,    digits=flowdigits),
                    flowo=round(Qo,       digits=flowdigits),
                    error=Qerr,
                    velheadfree=round(Vels, digits=flowdigits),
                    Hfree=round(h + Vels, digits=flowdigits),
                    ht=ht,
                    L=L,
                    b.over.B=bBs,
                    h.over.L=h/L,
                    h.over.P=h/P,
                    ht.over.H=round(htHs, digits=flowdigits),
                    H.over.P=round((h + Vels)/P, digits=flowdigits),
                    C=round(Cs,    digits=coedigits),
                    kc=round(kcs,  digits=coedigits),
                    kt=round(kts,  digits=coedigits),
                    source="weir3a5.sharpcrest",
                    message=messages);
  } else {
    z <- data.frame(head=h,
                    flow=round(Qs*kts,  digits=flowdigits),
                    flowfree=round(Qs,  digits=flowdigits),
                    flowo=round(Qo,     digits=flowdigits),
                    velheadfree=round(Vels, digits=flowdigits),
                    C=round(Cs,    digits=coedigits),
                    kc=round(kcs,  digits=coedigits),
                    kt=round(kts,  digits=coedigits),
                    source="weir3a5.sharpcrest",
                    message=messages);
  }
  att <- attributes(z);
  att$header <- header;
  attributes(z) <- att;
  return(z);
}
