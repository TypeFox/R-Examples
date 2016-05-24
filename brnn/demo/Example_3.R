#Example 3
#2 Inputs and 1 output
#the data used in Paciorek and
#Schervish (2004). The data is from a two input one output function with Gaussian noise
#with mean zero and standard deviation 0.25

data(twoinput)

#Formula interface
out=brnn(y~x1+x2,data=twoinput,neurons=10)

#With the default S3 method
#out=brnn(y=as.vector(twoinput$y),x=as.matrix(cbind(twoinput$x1,twoinput$x2)),neurons=10)
   
f=function(x1,x2) predict(out,cbind(x1,x2))
x1=seq(min(twoinput$x1),max(twoinput$x1),length.out=50)
x2=seq(min(twoinput$x2),max(twoinput$x2),length.out=50)
z=outer(x1,x2,f) # calculating the density values

transformation_matrix=persp(x1, x2, z, 
                            main="Fitted model", 
                            sub=expression(y==italic(g)~(bold(x))+e),
                            col="lightgreen",theta=30, phi=20,r=50, d=0.1,
                            expand=0.5,ltheta=90, lphi=180, 
                            shade=0.75, ticktype="detailed",nticks=5) 
points(trans3d(twoinput$x1,twoinput$x2, f(twoinput$x1,twoinput$x2), 
               transformation_matrix), col = "red") 
