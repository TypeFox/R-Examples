`[[<-.setupSNP` <-
function(x,i,j,info,value){
    out<-NextMethod("[[")

    lab<- attr(x, "label.SNPs")
    val<- names(out)[!(names(out) %in% lab) & !(names(out) %in% names(x))]
    if(!is.null(val)){ # new var
       inf<- attr(x, "gen.info")
       if(length(val)>1) warning("More than 1 column added. May break consistency of colSNPs", call.=FALSE)
       if (is.snp(value)){
          lab<-c(lab,val[1])
          if(!is.null(inf)) {
             if (missing(info)){
                 info<-rep(NA,ncol(inf))
                 warning("info was filled with NA")
             }
             inf<-rbind(inf,info)
          }
       }
    }
    if (!is.null(dim(out))){
      k<- match(lab, names(out)) # nuevas columnas con snps
      k<-k[!is.na(k)]
      ik<- match(names(out), lab) # nuevas columnas con snps
      ik<-ik[!is.na(ik)]
      
      attr(out, "colSNPs")        <- sort(k)
      attr(out, "label.SNPs")     <- lab[ik]
      attr(out, "gen.info")       <- inf[ik,]
    }
    out }

