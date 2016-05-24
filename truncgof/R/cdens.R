"cdens" <-
function(distn, H)
{
    if (!is.function(try(get(distn), silent = TRUE)))
       stop("'distn' must be a character of a distribution function")
    pdens <- distn
    ddens <- paste("d", substring(distn, 2), sep = "")
    nd <- formals(get(ddens))
    np <- formals(get(pdens))
    mn <- match(names(nd), names(np))
    args <- c(np[!is.na(mn)], nd[which("log" == names(nd))])
    argd <- paste(names(args), collapse = ", ")
    argp <- paste(head(names(args),length(args)-1), collapse = ", ")
    f <- function(){}
    formals(f) <- c(alist(x=), args, list(H = H))
    if(substring(distn, 2)=="gamma"){
      body(f) <- parse(text = "
        {i <- which(x <= H);
         logarg <- \"log\" %in% names(as.list(match.call()));
         if(!missing(rate)){
          den <- 1 - pgamma(H, shape, rate)
          if (logarg && log)
           res <- dgamma(x, shape, rate, log=log) - log(den)
          else res <- dgamma(x, shape, rate, log=log)/den
          res[i] <- 0
         }else{
          den <- 1 - pgamma(H, shape, scale)
          if (logarg && log)
          res <- dgamma(x, shape, scale, log=log) - log(den)
         else res <- dgamma(x, shape, scale, log=log)/den
         res[i] <- 0      
         };
         return(res)}")
    }else{
      body(f) <- parse(text = paste(
                         "{i <- which(x <= H);
        logarg <- \"log\" %in% names(as.list(match.call()));
        den <- 1 -", pdens, "(H, ", argp, ");
        if (logarg && log) res <-", ddens, "(x, ", argd, ") - log(den)
        else res <-", ddens, "(x, ", argd, ")/den;
        res[i] <- 0; return(res)}"))
    }
    return(f)
}
