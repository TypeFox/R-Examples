BootsChapt<-function(x, stat1, stat2 = NULL, B, replace = FALSE, alternative = c("one.change", "epidemic"), 
        adj.Wn = FALSE, tol = 1.0e-7, maxit = 50,trace = FALSE,... )
{
    if(B>=10) cat("It will take a while to do the bootstrap test(s).\n Please be patient.\n")
    alternative <- match.arg(alternative)
    nc<-0
    ifelse(alternative == "epidemic", nc<-2, nc<-1)
    if(is.vector(x)){
        n<-length(x)
        z<-sort(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        n<-nrow(x)
        z<-sort(x[,1])
    }
    if(nc == 1){
        if(stat1 < 0 && (is.null(stat2) || stat2<0)) stop("no valid input value for Tn")
        if(stat1 >= 0) Tn<-stat1
        else{ 
            if(!is.null(stat2) && stat2>=0){ Tn<-stat2
                warning("input stat2 is taken as Tn")}
        }
        Freq<-0
    }
    else{ 
        if(stat1 < 0 || (is.null(stat2) || stat2<0)) stop("no valid input value for Vn or Wn")
        Vn<-stat1; Wn<-stat2
        Freq.Vn<-0; Freq.Wn<-0}
    for(b in 1:B){
        y<-x[sample(1:n, n,replace = replace ),]
        Temp<-SemiparChangePoint(y, alternative = alternative, adj.Wn = adj.Wn, tol = tol, maxit = maxit,trace = F)
        if(nc == 1) Freq<-Freq+(Temp$Sn>Tn)
        else{
            Freq.Vn<-Freq.Vn+(Temp$Vn>Vn); 
            Freq.Wn<-Freq.Wn+(Temp$Wn>Wn)}
#        cat("Iteration ",b,"k-hat=",Temp$k.hat,"P-value=",Freq/b,"\n")
    }
    if(nc==1){ p.boots<-Freq/B; names(p.boots)<-"bootstrap p-value of Tn"; return(p.boots)}
    else{ Freq.Vn<-Freq.Vn/B; Freq.Wn<-Freq.Wn/B
        names(Freq.Vn)<-"bootstrap p-value of Vn"
        names(Freq.Wn)<-"bootstrap p-value of Wn"
        list(p.boots.Vn = Freq.Vn, p.boots.Wn = Freq.Wn)}
}
