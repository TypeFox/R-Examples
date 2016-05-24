"ss.nonadh" <-
function(mu0 = NULL, mu1 = NULL, delta=NULL, sigma0.sq=1, rho0=0, rho1=0, ss.ratio=1, var.ratio=1, deltaB=0, sig.level = 0.05, power = .80, 
    alternative = c("two.sided", "one.sided"),refinement=NULL,error.fisher=10^-6) 
{
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    if (tside==2 && deltaB !=0){stop("two-sided test should have deltaB=0") }
    if (!is.null(refinement)){
        rchoices<-c("Normal","Bernoulli","Fisher.exact")
        if (refinement=="normal") refinement<-"Normal"
        else if (refinement=="bernoulli") refinement<-"Bernoulli"
        else if (refinement=="fisher.exact") refinement<-"Fisher.exact"
        refinement<-char.expand(refinement,rchoices,
              stop("Refinement choices are:\n NULL, 'Normal', 'Bernoulli', and 'Fisher.exact'\n(upper case first letter is OK)"))
        }
    NOTE <- "n0 is number in *control* group\n\tn1 is number in *treatment* group"
    METHOD <- "Two-sample difference in means sample size with nonadherence"
        if (is.null(mu0) && is.null(mu1) && is.null(delta)){
        stop("either give 'mu0' and 'mu1' or give 'delta' ")
        }
        else if (is.null(delta)){ delta<- mu0-mu1 }
        else if (!is.null(mu0) && !is.null(mu1) && !is.null(delta)){
            if ( (mu0-mu1) != delta ){ stop(" mu0-mu1 should equal delta, or either delta should be null or mu0 and mu1 should be null") }
        }
        delta.star<- (1-rho0-rho1)*delta
        delta.star.minus.deltaB<- delta.star - deltaB
        if (delta.star> deltaB){ 
            if (tside==1) alternative<-paste(alternative,", delta > ",deltaB,sep="")   
             } 
        else if (delta.star< deltaB) {
            delta.star.minus.deltaB<- abs( delta.star.minus.deltaB) 
            if (tside==1) alternative<-paste(alternative,", delta < ",deltaB,sep="")   
             }
        else{ stop("delta.star - deltaB cannot equal 0") }
        sigma1.sq<- sigma0.sq* var.ratio
        V0<- (1-rho0)*sigma0.sq + rho0*sigma1.sq + rho0*(1-rho0)*delta^2
        V1<- (1-rho1)*sigma1.sq + rho1*sigma0.sq + rho1*(1-rho1)*delta^2

        Za<- qnorm(1-sig.level/tside)
        Zb<- qnorm(power)
 
    if (is.null(refinement) ){
        N0<- ((V0 + V1/ss.ratio)*(Za+Zb)^2) / delta.star.minus.deltaB^2
       }
    else if (refinement=="Normal" ){
        N0.start<- floor( ((V0 + V1/ss.ratio)*(Za+Zb)^2) / ((1-rho0-rho1)*delta-deltaB )^2)
        rootfunc1<-function(n0,alpha=sig.level/tside,v0=V0,v1=V1,r=ss.ratio,nominal.power=power,Delta.star.minus.deltaB=delta.star.minus.deltaB){
            ### calculate the power by t distribution minus nominal power
            n0<-ceiling(n0)
            n1<-ceiling(n0*r)
            tau.sq<- (v0/n0 + v1/n1)
            nu<- ((tau.sq)^2)/( v0^2/((n0-1)*n0^2) + v1^2/((n1-1)*n1^2 ) ) 
            nu[nu<1]<-1
            ncp<- Delta.star.minus.deltaB/ sqrt(tau.sq) 
            POWER<-1-pt( qt(1-alpha,nu), nu,  ncp )
            #temp<-c(n0,n1,tau.sq,nu,POWER)
            #names(temp)<-c("n0","n1","tau.sq","nu","POWER")
            #print(temp)
            nominal.power-POWER
          } 
        if (rootfunc1(N0.start)<0) N0.start<-1
        N0<-uniroot.integer(rootfunc1, lower =N0.start, upper =Inf,step.power=1)$root
       }
    else if (refinement=="Bernoulli"  || refinement=="Fisher.exact" ){
       if (is.null(mu0)) stop("For binary data, need a value for mu0")
       if (is.null(mu1)) mu1<- mu0-delta
       if (!is.null(sigma0.sq)){
	     if (sigma0.sq !=1) warning("sigma0.sq input value not used\n because refinement='Bernoulli' or 'Fisher.exact' ")
           sigma0.sq<-NULL
     }
       mu0.star<- mu0 - rho0*(mu0-mu1)
       mu1.star<- mu1 + rho1*(mu0-mu1)
       V0<- mu0.star*(1-mu0.star) 
       V1<- mu1.star*(1-mu1.star)
       N0<- ((V0 + V1/ss.ratio)*(Za+Zb)^2) / delta.star.minus.deltaB^2
       if (refinement=="Bernoulli" && tside==2){
          NOTE<-paste(NOTE,"\n No continuity correction.\n Consider Using refinement=Fisher.exact,\n using 'Esc' if computation time is too long.")
       }
       if (refinement=="Fisher.exact"){
           if (tside==1) stop("Refinement 'Fisher.exact' only defined when alternative='two-sided'")
           METHOD<-"Sample Size with Fisher's Exact and Nonadherence"
           N0.start<- max(1,ceiling(N0)-1)
           print("Intermediate Calculations:POWER is exact power of Fisher's exact test")
           rootfunc2<-function(n0,p0=mu0.star,p1=mu1.star,alpha=sig.level,nominal.power=power,r=ss.ratio,eps=error.fisher){
		    prob.reject<-0
                n0<-ceiling(n0)
                n1<-ceiling(n0*r)
		    ilow<-qbinom(eps/4,n0,p0)
		    ihigh<-qbinom(1-eps/4,n0,p0)
	          jlow<-qbinom(eps/4,n1,p1)
		    jhigh<-qbinom(1-eps/4,n1,p1)
		    for (i in ilow:ihigh){
 				for (j in jlow:jhigh){
				    if (fisher.test(matrix(c(n0-i,i,n1-j,j),2,2))$p.value<=alpha){ 
                               prob.reject<-prob.reject+ dbinom(i,n0,p0)*dbinom(j,n1,p1) }
				}
			}
                temp<-c(n0,n1,prob.reject)
                names(temp)<-c("n0","n1","POWER")
                print(temp)
                nominal.power - prob.reject
            }
           if (rootfunc2(N0.start)<0) N0.start<-1
          N0<-uniroot.integer(rootfunc2, lower =N0.start, upper =Inf,step.power=3)$root
        } ### end of Fisher.exact if statement
       } ### end of Bernoulli or Fisher.exact if statement
    else { stop("refinement should equal either NULL, 'normal', 'Bernoulli' or 'Fisher.exact'") }

     output<-list(n0 = ceiling(N0), n1=ceiling(ss.ratio*N0), 
        mu0=mu0,
        mu1=mu1,
        delta = delta,
        sigma0.sq = sigma0.sq, var.ratio=var.ratio, 
        #ss.ratio=ss.ratio, 
        rho0=rho0, rho1=rho1, sig.level = sig.level, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD,refinement=refinement)

    structure(output,class = "power.htest")


}

