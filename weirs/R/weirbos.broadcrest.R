"weirbos.broadcrest" <-
function(h, ht=NULL, b=NULL, B=NULL, P=NULL, L=NULL,
         R=0, r=0, A=NULL, alpha=1,
         slopeus="vertical", slopeds="vertical",
         kc=NULL, kr=NULL, ks=NULL, C=NULL,
         contractratio=NULL,
         hhptest=TRUE, extended=TRUE,
         header="",
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

  if(slopeds != "vertical") {
    if(as.logical(length(grep(":",slopeds)))) {
      slopeds <- as.numeric(unlist(strsplit(slopeds,":")));
      slopeds <- slopeds[1]/slopeds[2];
    } else {
      stop("slopeds does not contain a colon (:)");
    }
  } else {
    slopeds <- 0;
  }

  if(slopeus < 0) stop("upstream slope can not be negative");
  if(slopeds < 0) stop("downstream slope can not be negative");


  if(is.null(h)) {
     stop("head is NULL");
  }
  if(is.null(ht)) ht <- rep(0, length(h));
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
  if(! is.null(R) & (R < 0 | R > min(L))) {
     stop("implausible radius of curvature on crest");
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

  if(! is.null(kr)) {
    if(length(kr) == 1) kr <- rep(kr, length(h));
    if(length(kr) != length(h)) {
      stop("rounding coefficient is not equal to length of head vector");
    }
  }

  if(! is.null(ks)) {
    if(length(ks) == 1) ks <- rep(ks, length(h));
    if(length(ks) != length(h)) {
      stop("downstream slope coefficient is not equal to length of head vector");
    }
  }

  if(! is.null(C)) {
    if(length(C) == 1) C <- rep(C, length(h));
    if(length(C) != length(h)) {
      stop("discharge coefficient is not equal to length of head vector");
    }
  }

  Qs <- Qo <- Qcv <- Qerr <- vector(mode="numeric", length=length(h));
  messages <- Cs <- Cvs <- kcs <- krs <- kss <- Vels <- bBs <- Qs;
  g <- 32.2; g2 <- 2*g;
  BosTerms <- (2/3)*sqrt(2/3*g)*b;


  Rhkr <- get("broadcrest.roundingtable", .weir.nomographs);

  BosC  <- get("Bos1989C",  .weir.nomographs);
  BosCv <- get("Bos1989Cv", .weir.nomographs);

  for(i in 1:length(h)) {
    it <- 0;
    messages[i] <- "ok"; # assumed for initialization
  	hh <- h[i];
  	R.over.h <- R/hh;
  	r.over.b <- r/b;

    Astar <- hh*b;

  	if(ht[i]/hh >= 0.66) { # Sturm (2001,p.55) attributed to Bos(1988).
  	  Qo[i]   <- NA;
      Qcv[i]  <- NA;
      Qs[i]   <- NA;
      Qerr[i] <- NA;
      Vels[i] <- NA;
      Cs[i]   <- NA;
      Cvs[i]  <- NA;
      kcs[i]  <- NA;
      krs[i]  <- NA;
      kss[i]  <- NA;
      messages[i] <- "submerged";
  	  next;
  	}
  	if(hh == 0) {
  	  Qo[i]   <- 0.00;
      Qcv[i]  <- 0.00;
      Qs[i]   <- 0.00;
      Qerr[i] <- 0.00;
      Vels[i] <- 0.00;
      Cs[i]   <- NA;
      Cvs[i]  <- NA;
      kcs[i]  <- NA;
      krs[i]  <- NA;
      kss[i]  <- NA;
      messages[i] <- "head zero";
  	  next;
  	}
    if(hhptest & hh/(hh+P) > 0.35) {
  	  Qo[i]   <- NA;
      Qcv[i]  <- NA;
      Qs[i]   <- NA;
      Qerr[i] <- NA;
      Vels[i] <- NA;
      Cs[i]   <- NA;
      Cvs[i]  <- NA;
      kcs[i]  <- NA;
      krs[i]  <- NA;
      kss[i]  <- NA;
      messages[i] <- "C is indeterminent";
  	  next;
    }
    if(is.null(kr)) {
      the.kr <- approx(Rhkr$R.over.h, Rhkr$kr, R.over.h, rule=2)$y
    } else {
      the.kr <- kr[i];
    }

    h.over.L <- hh/L[i];
    h.over.P <- hh/P;

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

    if(h.over.L >= h.over.L.critical) {
       Qo[i]   <- NA;
       Qcv[i]  <- NA;
       Qs[i]   <- NA;
       Qerr[i] <- NA;
       Vels[i] <- NA;
       Cs[i]   <- NA;
       Cvs[i]  <- NA;
       kcs[i]  <- NA;
       krs[i]  <- NA;
       kss[i]  <- NA;
       messages[i] <- "sharp-crested";
       next;
    }

    if(is.null(kc)) {
      # Now that we know broad-crested flow---determine kc
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
        Qcv[i]  <- NA;
        Qs[i]   <- NA;
        Qerr[i] <- NA;
        Vels[i] <- NA;
        Cs[i]   <- NA;
        Cvs[i]  <- NA;
        kcs[i]  <- NA;
        krs[i]  <- NA;
        kss[i]  <- NA;
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
        the.kc <- 1.00; # TWRI3A5p.10
      } else if(r.over.b > 0) {
        the.kc <- approx(c(0,0.12), c(the.kc,1), r.over.b, rule=2)$y;
      }
    } else {
      the.kc <- kc[i];
    }


    if(is.null(C)) {
      if(h.over.L < 0.08) {
        Qo[i]   <- NA;
        Qcv[i]  <- NA;
        Qs[i]   <- NA;
        Qerr[i] <- NA;
        Vels[i] <- NA;
        Cs[i]   <- NA;
        Cvs[i]  <- NA;
        kcs[i]  <- NA;
        krs[i]  <- NA;
        kss[i]  <- NA;
        messages[i] <- "h.over.L too small to determine C";
        next;
      } else if(h.over.L > 1.5) {
        Qo[i]   <- NA;
        Qcv[i]  <- NA;
        Qs[i]   <- NA;
        Qerr[i] <- NA;
        Vels[i] <- NA;
        Cs[i]   <- NA;
        Cvs[i]  <- NA;
        kcs[i]  <- NA;
        krs[i]  <- NA;
        kss[i]  <- NA;
        messages[i] <- "h.over.L too small to determine C";
        next;
      } else {
        the.C <- approx(BosC$h.over.L, BosC$C, h.over.L, rule=2)$y;
      }
    } else {
      the.C <- C[i];
    }

    the.A   <- ifelse(is.null(A), (hh+P)*B, A[i]);
    CAstar.over.A <- the.C*Astar/the.A;

    the.Cv <- approx(BosCv$CAstar.over.A, BosCv$Cv,
                     CAstar.over.A, rule=2)$y;
    #message("the.Cv",the.Cv);
    #message("the.A", the.A);
    #message("Astar", Astar);


    if(is.null(ks)) {
      if(slopeds > 1) {
        if(slopeds < 2) {
          the.ks <- 1;
        } else if(slopeds > 5) {
          Qo[i]   <- NA;
          Qcv[i]  <- NA;
          Qs[i]   <- NA;
          Qerr[i] <- NA;
          Vels[i] <- NA;
          Cs[i]   <- the.C;
          Cvs[i]  <- the.Cv;
          kcs[i]  <- the.kc;
          krs[i]  <- the.kr;
          kss[i]  <- NA;
          messages[i] <- "slopeds too shallow to determine ks";
          next;
        } else {
          tmpe      <- get("broadcrest.downstreamtable", .weir.nomographs);
          slopes    <- as.numeric(get("slopes",   tmpe));
          h.over.Ls <- get("h.over.L", tmpe);
          tmp       <- slopes[slopes <= slopeds];
          if(length(tmp) != 0) slopes.min <- max(tmp);
          tmp       <- slopes[slopes >= slopeds];
          if(length(tmp) != 0) slopes.max <- min(tmp);
          if(slopes.min == slopes.max) {
            tmp    <- get(sprintf("%.4f", slopes.min), tmpe);
            the.ks <- approx(h.over.Ls, tmp, h.over.P, rule=2)$y;
          } else {
            tmp    <- get(sprintf("%.4f", slopes.min), tmpe);
            ks.min <- approx(h.over.Ls, tmp, h.over.L, rule=2)$y;
            tmp    <- get(sprintf("%.4f", slopes.max), tmpe);
            ks.max <- approx(h.over.Ls, tmp, h.over.L, rule=2)$y;
            the.ks <- approx(c(slopes.min, slopes.max),
                             c(ks.min, ks.max), slopeds)$y;
          }
        }
      } else {
        the.ks <- 1;
      }
    } else {
      the.ks <- ks[i];
    }

    Cs[i]  <- the.C;
    kcs[i] <- the.kc;
    krs[i] <- the.kr;
    kss[i] <- the.ks;
    Cvs[i] <- the.Cv;

    bigOval <- the.kc*the.kr*the.ks*BosTerms;
    Qint    <- bigOval*the.C*hh^1.5;
    Qbycv   <- Qint*the.Cv;
    Qo[i]   <- Qint;
    Qcv[i]  <- Qbycv;
    Qold    <- Qbycv;
    the.alpha <- alpha[i];
    "afunc" <- function(Q) {
       vel     <- Q/the.A;
       vhead   <- vel^2/g2;
       Vels[i] <- vhead;
       H <- hh + alpha[i]*vhead;
       Qtmp <- bigOval*the.C*1*H^1.5; # 1 is the Cv=1
       return(Qtmp);
    }
    while(1) {
      it <- it + 1;
      Qnew <- afunc(Qold);
      if(Qnew == Inf) {
        Qo[i]   <- NA;
        Qs[i]   <- NA;
        Qcv[i]  <- NA;
        messages[i] <- "nonconvergence";
        break;
      }
      err <- abs(Qnew - Qold);
      Qerr[i] <- err;
      Vels[i] <- ifelse(is.null(A), Qnew/((hh+P)*B), Qnew/A[i]);
      if(it > maxit || err < eps) break;
      Qold <- Qnew;
    }
    Qs[i] <- Qnew;
  }
  if(extended) {
    z <- data.frame(head=h,
                    flow=round(Qs,      digits=flowdigits),
                    delta=c(NA,diff(Qs)),
                    flowo=round(Qo,     digits=flowdigits),
                    flowcv=round(Qcv,   digits=flowdigits),
                    error=Qerr,
                    velhead=round(Vels, digits=flowdigits),
                    H=round(h+Vels,     digits=flowdigits),
                    ht=ht,
                    L=L,
                    b.over.B=bBs,
                    h.over.L=h/L,
                    h.over.P=h/P,
                    C=round(Cs,    digits=coedigits),
                    Cv=round(Cvs,  digits=coedigits),
                    kc=round(kcs,  digits=coedigits),
                    kr=round(krs,  digits=coedigits),
                    ks=round(kss,  digits=coedigits),
                    source="weirbos.broadcrest",
                    message=messages);
  } else {
    z <- data.frame(head=h,
                    flow=round(Qs,      digits=flowdigits),
                    flowo=round(Qo,     digits=flowdigits),
                    flowcv=round(Qcv,   digits=flowdigits),
                    velhead=round(Vels, digits=flowdigits),
                    C=round(Cs,    digits=coedigits),
                    Cv=round(Cvs,  digits=coedigits),
                    kc=round(kcs,  digits=coedigits),
                    kr=round(krs,  digits=coedigits),
                    ks=round(kss,  digits=coedigits),
                    source="weirbos.broadcrest",
                    message=messages);
  }
  att <- attributes(z);
  att$header <- header;
  attributes(z) <- att;
  return(z);
}
