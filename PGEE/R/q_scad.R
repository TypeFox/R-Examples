q_scad <-
function(theta,lambda,a=3.7)
       {
        #length of parameter
        p<-length(theta)
        #get the absolute value
        theta<-abs(theta)
        #create vector of zeros
        b1<-rep(0,p)
        #if theta is greater then lambda set it to 1
        b1[theta>lambda]<-1
        #create an another vector of zeros
        b2<-rep(0,p)
        #if theta is less than a*lambda, set it to 1.
        b2[theta<(lambda*a)]<-1
        lambda*(1-b1)+((lambda*a)-theta)*b2/(a-1)*b1
}
