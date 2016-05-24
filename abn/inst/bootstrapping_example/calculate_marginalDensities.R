##################################################################################################################
## Estimate all posterior estimates from a given DAG
##################################################################################################################
# Step 1. Define the DAG
##################################################################################################################
library(abn);
# define a DAG and call it mydag
mydag<-matrix(c(
  # D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 Year Loc.x Loc.y
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D1
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D2
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D3
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   1,    0,    # D4
  0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0,   0,    0,    # D5
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D6
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,   0,    0,    # D7
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0,   0,    0,    # D8
  0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # D9
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0,   0,    0,    # D10
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,   0,    0,    # Year
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0,    # Loc.x
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,   0,    0     # Loc.y
),byrow=TRUE,ncol=13); 
colnames(mydag)<-rownames(mydag)<-names(pigs.1par);#set names

## setup distribution list for each node
mydists.pigs <- list( D1 = "binomial",
                      D2 = "binomial",
                      D3 = "binomial",
                      D4 = "binomial",
                      D5 = "binomial",
                      D6 = "binomial", 
                      D7 = "binomial",
                      D8 = "binomial",
                      D9 = "binomial",
                      D10 = "binomial",
                      Year = "gaussian",
                      Loc.x = "gaussian",
                      Loc.y = "gaussian"
);

## use fitabn to print mlik and check if mydag is the correct model 
print(fitabn(dag.m=mydag,data.df=pigs.1par,data.dists=mydists.pigs)$mlik);
## if correct mlik.pigs <- -79118.38
##################################################################################################################
# Step 2. Compute marginals for each and every parameter in the model
#################################################################################################
## N.B.: with compute.fixed=TRUE this tries to estimate the posterior densities. Note that the
## algorithm used to locate the x values is very crude and it may be better to provide these
## manually in some case. Using n.grid=1000 uses a fine grid of 1000 points. 
## this works well here for the "sensible" variables but terrible for the "bad" variables
##################################################################################################################
marg<-fitabn(dag.m=mydag,data.df=pigs.1par,data.dists=mydists.pigs,compute.fixed=TRUE,n.grid=1000);#, marginal.node=1,marginal.param=1,variate.vec=seq(-150,5,len=1000),verbose=TRUE);

### Reproduce fitabn in order to get a grid of discrete values d and f(d):
marnew <- marg$marginals[[1]]

for(i in 2: length(marg$marginals)){
  marnew <- c(marnew, marg$marginals[[i]])
}

### Now create the DATA for going to JAGS!!!
print(names( marnew));
##  we want to bind all the marginals the same nodes into a matrix
m <- marnew;

D1.p <- cbind( m[[ "D1|(Intercept)"]], m[["D1|D2"]]);
D2.p <- cbind( m[[ "D2|(Intercept)"]], m[["D2|D3"]]);
D3.p <- cbind( m[[ "D3|(Intercept)"]], m[["D3|D4"]]);
D4.p <- cbind( m[[ "D4|(Intercept)"]], m[["D4|Loc.x"]]);
D5.p <- cbind( m[[ "D5|(Intercept)"]], m[["D5|D6"]]);
D6.p <- cbind( m[[ "D6|(Intercept)"]], m[["D6|D4"]]);
D7.p <- cbind( m[[ "D7|(Intercept)"]], m[["D7|Year"]]);
D8.p <- cbind( m[[ "D8|(Intercept)"]], m[["D8|D10"]]);
D9.p <- cbind( m[[ "D9|(Intercept)"]], m[["D9|D2"]]);
D10.p <- cbind( m[[ "D10|(Intercept)"]], m[["D10|D9"]]);
Year.p <- cbind( m[[ "Year|(Intercept)"]]); 
prec.Year.p <- cbind( m[[ "Year|precision" ]]);
Loc.x.p <- cbind( m[[ "Loc.x|(Intercept)"]], m[["Loc.x|D7"]]);
prec.Loc.x.p <- cbind( m[[ "Loc.x|precision"]]);
Loc.y.p <- cbind( m[[ "Loc.y|(Intercept)"]], m[["Loc.y|D7"]]);
prec.Loc.y.p <- cbind( m[[ "Loc.y|precision"]]);

###############################################################################################
# Step 3. dump all the densities into a format which can be read be JAGS
# dump this into a file called post_params.R
###############################################################################################
## now dump it all to a format which JAGS then import

dump(c("D1.p","D2.p","D3.p","D4.p","D5.p","D6.p","D7.p","D8.p","D9.p","D10.p",
       "Year.p","prec.Year.p","Loc.x.p","prec.Loc.x.p","Loc.y.p","prec.Loc.y.p"),
     file="post_params.R");

proc.time()