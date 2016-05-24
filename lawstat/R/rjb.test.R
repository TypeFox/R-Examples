`rjb.test` <-
function (x, option=c("RJB", "JB"), crit.values=c("chisq.approximation", "empirical"), N=0)
{

option<-match.arg(option) 
crit.values=match.arg(crit.values)


   if (NCOL(x) > 1){stop("x is not a vector or univariate time series")}

   if (any(is.na(x))){stop("NAs in x")}

   if ((crit.values=="empirical")&(N==0)) 
   {stop("number of Monte Carlo simulations N should be provided for the empirical critical values")}

   DNAME <- deparse(substitute(x))


    ### Calculate the first 4 central moments ###
    n <- length(x)
    m1 <- sum(x)/n
    m2 <- sum((x - m1)^2)/n
    m3 <- sum((x - m1)^3)/n
    m4 <- sum((x - m1)^4)/n


    ### User can choose the Standard Jarque Bera Test or Robust Jarque Bera Test ###
    ### Robust Jarque Bera Test is default ###
    if(option=="JB")
    {
    b1 <- (m3/m2^(3/2))^2;
    b2 <- (m4/m2^2);
    METHOD <- "Standard Jarque Bera Test"
    statistic<-n * b1/6 + n * (b2 - 3)^2/24
    }

    else
    {
    option="RJB";
    J<-sqrt(pi/2)*mean(abs(x-median(x)));
    J2<-J^2;
    b1 <- (m3/(J2)^(3/2))^2;
    b2 <- (m4/(J2)^2);
    vk<-64/n;
    METHOD <- "Robust Jarque Bera Test"
    vs<-6/n
    ek<-3
    statistic <- b1/vs + (b2 - ek)^2/vk
    }

 
if(crit.values=="empirical")
  {
     if(option=="JB")
     {
       #### computes empirical critical values for the JB statistic####

   jb<-double(N)

   for (k in 1:N)
      {
       e <- rnorm(length(x), mean=0, sd = sqrt(1))  
       m1 <- sum(e)/n
       m2 <- sum((e - m1)^2)/n
       m3 <- sum((e - m1)^3)/n
       m4 <- sum((e - m1)^4)/n
       b1 <- (m3/m2^(3/2))^2;
       b2 <- (m4/m2^2);
       vk<-24/n;
       vs<-6/n
       ek<-3
       jb[k]<-b1/vs + (b2 - ek)^2/vk
      }
   y<-sort(jb)
            if (statistic>=max(y)) {p.value=0}
             else if (statistic<=min(y)) {p.value=1}
             else
             {
  
   bn<-which(y==min(y[I(y>=statistic)]))
   an<-which(y==max(y[I(y<statistic)]))
   a<-max(y[I(y<statistic)])
   b<-min(y[I(y>=statistic)])
   pa<-(an - 1) / (N - 1)
   pb<-(bn - 1) / (N - 1)
   alpha<-(statistic-a)/(b-a)  
   p.value=1-alpha*pb-(1-alpha)*pa
             }
      }


     else
        {
      #### computes empirical critical values for the RJB statistic####

   rjb<-double(N)

   for (k in 1:N)
      {
       e <- rnorm(length(x), mean=0, sd = sqrt(1))  
       J<-sqrt(pi/2)*mean(abs(e-median(e)));
       J2<-J^2;
       m1 <- sum(e)/n
       m2 <- sum((e - m1)^2)/n
       m3 <- sum((e - m1)^3)/n
       m4 <- sum((e - m1)^4)/n
       b1 <- (m3/(J2)^(3/2))^2;
       b2 <- (m4/(J2)^2);
       vk<-64/n;
       vs<-6/n
       ek<-3
       rjb[k]<-b1/vs + (b2 - ek)^2/vk
      }
   y<-sort(rjb)
             if (statistic>=max(y)) {p.value=0}
             else if (statistic<=min(y)) {p.value=1}
             else
             {
   bn<-which(y==min(y[I(y>=statistic)]))
   an<-which(y==max(y[I(y<statistic)]))
   a<-max(y[I(y<statistic)])
   b<-min(y[I(y>=statistic)])
   pa<-(an - 1) / (N - 1)
   pb<-(bn - 1) / (N - 1)
   alpha<-(statistic-a)/(b-a)  
   p.value=1-alpha*pb-(1-alpha)*pa
              }
          }

}

 

else {p.value <- 1 - pchisq(statistic, df = 2)}
if(option=="JB"){METHOD <- "Jarque Bera Test"}
else {METHOD <- "Robust Jarque Bera Test"}

    ### Display Output ###
    STATISTIC=statistic
    names(STATISTIC) <- "X-squared"
    PARAMETER <- 2
    names(PARAMETER) <- "df"
    
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
        p.value = p.value, method = METHOD, data.name = DNAME),
        class = "htest")
}

