agreement <-
function(x,y,error,target,CCC_a=0.95,TDI_a,alpha=0.05,CP_a=0.9,H_label="x",V_label="y",min=NA,max=NA,by=NA,dec=2){
if (tolower(error)=="prop"){
    if ((min(y) & min(x))<=0){
        stop("For proportional error, all x and y should be greater than 0");
    }
    y <- log(y);
    x <- log(x);
}

if (tolower(error)!="prop"){
    if (is.na(min)) min <- min(min(y),min(x));
    if (is.na(max)) max <- max(max(y),max(x));
    if (is.na(by)) by <- round((max-min)/5,dec);
}else{
    if (is.na(min)) min <- exp(min(min(y),min(x)));
    if (is.na(max)) max <- exp(max(max(y),max(x)));
    if (sum(is.na(by))==1) by <- round(seq(min,max,length=7),dec);
}

d <- y-x;
m1 <- mean(y);
m2 <- mean(x);
yyy <- mean(d);
v1 <- var(y);
v2 <- var(x);
yy <- var(d);
xx <- sum(y^2);
xxx <- sum(x^2);
e2 <- sum(d^2);
n <- length(d);

if (length(y)!=length(x)) stop("x and y should have the same length!");

mu_d <- m1-m2; 
d2 <- mu_d^2;
e2 <- e2/n; 
sd2 <- (e2-d2)*n/(n-3); 
s12 <- v1*(n-1)/n; 
s22 <- v2*(n-1)/n; 
U <- mu_d/sqrt(sqrt(s12*s22)); 
V <- sqrt(s12/s22); 
Ca <- 2/(V+1/V+U**2); 
rc <- 1-e2/(d2+s12+s22); 
r <- (rc/Ca); 
d2_s2 <- (d2/sd2); 
b1 <- r*V; 
b0 <- m1-b1*m2; 
se2 <- s12*(1-r**2)*n/(n-3); 
e2 <- e2*n/(n-1);
kp <- TDI_a;
if (error=="prop"){
    kp <- log(TDI_a/100+1);
}

if (tolower(target)=="fixed"){
    di <- b0+(b1-1)*x; 
    d2i <- di^2;
    cpi <- pchisq(kp^2/se2,1,d2i/se2);
    se <- sqrt(se2); 
    kpm <- (kp+di)/se; 
    kmm <- (kp-di)/se;
    C0 <- dnorm(-kpm)-dnorm(kmm); 
    C1 <- C0*x;
    C2 <- -kpm*dnorm(-kpm)-kmm*dnorm(kmm);
    cp <- mean(cpi);
    c0 <- mean(C0);
    c1 <- mean(C1);
    c2 <- mean(C2);
    s_p <- c0^2+(c0*m2-c1)^2/s22+c2^2/2;
    if (cp<1){
        s_T <- sqrt((s_p/((n-3)*cp^2*(1-cp)^2)));
    }else{s_t <- NA;}
    if (cp<1){
        T <- log(cp/(1-cp));
    }else{T <- NA;}
    T_l <- T-qnorm(1-alpha)*s_T;
    cp_l <- exp(T_l)/(1+exp(T_l));
}

#For CCC, accuracy and precision;
L_ca <- log(Ca/(1-Ca));
z_rc <- log((1+rc)/(1-rc))/2; 
z_r <- log((1+r)/(1-r))/2; 

#For MSD & approx TDI & CP;
W <- log(e2);
k <- qnorm(1-(1-CP_a)/2)*sqrt(e2);
cp_a <- pchisq(kp^2/e2,1); 

if (tolower(target)=="fixed"){
    SZ_rc <- (1-r^2)*rc^2/((n-2)*(1-rc^2)^2*r^2)*(V*U^2*rc^2+(1-rc*r*V)^2+V^2*rc^2*(1-r^2)/2); 
    SZ_r <- (1-r^2/2)/(n-3);
    SL_ca <- (U^2*V*Ca^2*(1-r^2)+(1-V*Ca)^2*(1-r^4)/2)/((n-2)*(1-Ca)^2);
    if (s22==0){
        s_W <- 2*(1-d2^2/e2^2)/(n-2);
    }else{s_W <- 2*(1-(d2+s22*(1-b1)^2)^2/e2**2)/(n-2);}   
}else{
    SZ_rc <- ((1-r^2)*rc^2/(1-rc^2)/r^2+2*rc^3*(1-rc)*U**2/r/(1-rc^2)^2-rc^4*U^4/r^2/(1-rc^2)^2/2)/(n-2);
    SZ_r <- 1/(n-3);
    SL_ca <- ((Ca^2*U^2*(V+1/V-2*r)+Ca^2*(V^2+1/V^2+2*r^2)/2+(1+r**2)*(Ca*U^2-1))/((n-2)*(1-Ca)^2));
    s_W <- 2*(1-d2^2/e2^2)/(n-2);
}

SZ_rc <- sqrt(SZ_rc); 
SZ_r <- sqrt(SZ_r); 
SL_ca <- sqrt(SL_ca); 
s_W <- sqrt(s_W);
z_rc_l <- z_rc-qnorm(1-alpha)*SZ_rc; 
z_r_l <- z_r-qnorm(1-alpha)*SZ_r;
l_ca_l <- L_ca-qnorm(1-alpha)*SL_ca;
rc_l <- tanh(z_rc_l); 
r_l <- tanh(z_r_l); 
ca_l <- exp(l_ca_l)/(1+exp(l_ca_l));

W_u <- W+qnorm(1-alpha)*s_W;
e2_u <- exp(W_u);
k_u <- qnorm(1-(1-CP_a)/2)*sqrt(e2_u);
cp_a_l <- pchisq(kp^2/e2_u,1); 

if (tolower(error)=="prop"){
    k <- 100*(exp(k)-1); 
    k_u <- 100*(exp(k_u)-1);
}

if (tolower(target)!="fixed"){
    cp <- pchisq(kp^2/sd2,1,d2_s2); 
    sd <- sqrt(sd2); 
    kpm <- (kp+mu_d)/sd; 
    kmm <- (kp-mu_d)/sd;
    if (cp<1){ 
        s_T <- ((dnorm(-kpm)-dnorm(kmm))^2+(kmm*dnorm(kmm)+kpm*dnorm(-kpm))^2/2)/((n-3)*cp^2*(1-cp)^2); 
    }else {s_t <- NA;}  
    s_T <- sqrt(s_T);
    if (cp<1){
        T <- log(cp/(1-cp));
    }else{T <- NA;}
    T_l <- T-qnorm(1-alpha)*s_T;
    cp_l <- exp(T_l)/(1+exp(T_l));
}
conf <- (1-alpha)*100;

Data <- list(x=x,y=y,conf=conf,error=error,target=target,xlab=H_label,ylab=V_label,min=min,max=max,by=by,dec=dec)
Estimate=list(CCC=rc, Precision=r, Accuracy=Ca, TDI=k, CP=cp, RBS=d2_s2)
class(Estimate) <- "agreement"
Conf_Limit <- list(CCC=rc_l, Precision=r_l, Accuracy=ca_l, TDI=k_u, CP=cp_l, RBS=NA)
class(Conf_Limit) <- "agreement"
Allowance <- list(CCC=CCC_a, Precision=NA, Accuracy=NA, TDI=TDI_a, CP=CP_a, RBS=NA)
class(Allowance) <- "agreement"
Report <- list(Data=Data,Estimate=Estimate,Conf_Limit=Conf_Limit,Allowance=Allowance)
class(Report) <- "report"
return(Report)
}

