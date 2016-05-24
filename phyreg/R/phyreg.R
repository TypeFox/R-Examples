phyreg <-
function(control, test, data, subset, phydata, taxmatrix, heightsdata, rho=-1, lorho=0.3, hirho=0.6, errrho=0.02, minrho=0.0001, tolerance=1e-6, oppf=5, opdf=0, parmx=0, parmxz=0, opfunccall=0, addDF=0, linputs=FALSE, sinputs=FALSE, means=FALSE, lmshortx=FALSE, lmshortxz=FALSE, lmlongx=FALSE, lmlongxz=FALSE, hinput=FALSE, paper=FALSE, dfwarning=TRUE, oprho=FALSE, plot=FALSE, reset=FALSE) {

# ======= First define inphyreg and other functions, then get to body

 
inphyreg<-function(cont, intest, dataframe, insubset, phyvar, taxmat, inheights, rho=-1, lorho=0.1, hirho=0.9, errrho=0.08, minrho=0.01, 
tolerance=1e-6, addDF=0, oppf=5, opdf=0, oprho=0, parmx=0, parmxz=0, opfunccall=0, linputs=FALSE, sinputs=FALSE, means=FALSE, lmshortx=FALSE, lmshortxz=FALSE, lmlongx=FALSE, lmlongxz=FALSE, hinput=FALSE, dfwarning=TRUE, paper=FALSE, plot=FALSE) {


#====================================== definition of opsi 

##         NO -- NOW ONLY ONE defintion of opsi, at top level
 

## =============== This is ginv, taken from R 3..0.2 on 26 January 2014, to avoid having to load MASS

#ginv<-  function (X, tol = sqrt(.Machine$double.eps)) 
#{
#    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
#        stop("'X' must be a numeric or complex matrix")
#    if (!is.matrix(X)) 
#        X <- as.matrix(X)
#    Xsvd <- svd(X)
#    if (is.complex(X)) 
#        Xsvd$u <- Conj(Xsvd$u)
#    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
#    if (all(Positive)) 
#        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
#    else if (!any(Positive)) 
#        array(0, dim(X)[2L:1L])
#    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
#        t(Xsvd$u[, Positive, drop = FALSE]))
#}

VCMginv <- function(VCM, tolerance = sqrt(.Machine$double.eps))  ## retyped and shortened from the MASS ginv function
{VCMsvd<-svd(VCM)                                                                                    ## is VCM the right name?
 pos<-(VCMsvd$d > max( tolerance*VCMsvd$d[1L],0))
 if (all(pos)) VCMsvd$v%*% (1/VCMsvd$d * t(VCMsvd$u))
  else if (!any(pos)) array(0,c(dim(VCM)[1L],dim(VCM)[1L]))
    else VCMsvd$v[, pos, drop=FALSE] %*% ((1/VCMsvd$d[pos]) * t(VCMsvd$u[, pos,drop=FALSE]))
}


## ================================ Functions now defined inside inphyreg, namely findmax and getlik
##================================  and other functions internal to them

findmax<-function(xlo, xhi, xmin, xerr) {

xtriple<-c(1,3,79); ftriple<-c(1,3,79);
xtriple[1]<-xlo; xtriple[2]<-sqrt(xlo*xhi); xtriple[3]<-xhi
ftriple[1]<-getlik(xtriple[1]);  ftriple[2]<-getlik(xtriple[2]);  ftriple[3]<-getlik(xtriple[3]);
if(any(is.na(ftriple))) {cat("\nSome ftriple is NA\n"); browser()}

  # surround makes sure the maximum is in between the first and third elements of xtriple and ftriple
   surround<-function() {
     while ((ftriple[1]>ftriple[2]) && (xtriple[3]>xmin)){  
      newxlo<-xtriple[1]**2/xtriple[2];
      newflo<-getlik(newxlo);
      xtriple[3]<<-xtriple[2]; xtriple[2]<<-xtriple[1]; xtriple[1]<<-newxlo;
      ftriple[3]<<-ftriple[2]; ftriple[2]<<-ftriple[1]; ftriple[1]<<-newflo;
      }
    while (ftriple[3]>ftriple[2]) {  
      newxhi<-xtriple[3]**2/xtriple[2];
      newfhi<-getlik(newxhi);
      xtriple[1]<<-xtriple[2]; xtriple[2]<<-xtriple[3]; xtriple[3]<<-newxhi;
      ftriple[1]<<-ftriple[2]; ftriple[2]<<-ftriple[3]; ftriple[3]<<-newfhi;
     }
    }
   
   # try checks if the left or right half contains a higher value, and if so focusses in on that half
   try<- function(xtriple, ftriple, xmid, fmid, lopos, hipos) {        ## Notice I should change this name, as its an R function I use elsewhere in the code!
    if (fmid>ftriple[2]) {
      xtriplo<-xtriple[lopos]; xtriphi<-xtriple[hipos];
      xtriple[1]<<-xtriplo; xtriple[2]<<-xmid; xtriple[3]<<-xtriphi;
      ftriplo<-ftriple[lopos]; ftriphi<-ftriple[hipos];
      ftriple[1]<<-ftriplo; ftriple[2]<<-fmid; ftriple[3]<<-ftriphi;
      return(1);
      }
      else {return (0)}
    }
  
  # centre focusses in when the old centre is the new centre i.e. when the left and right halves don't show a higher value
   centre<- function(xtriple,ftriple,lxmid,lfmid,rxmid,rfmid) {
    xtriple[1]<<-lxmid; xtriple[3]<<-rxmid;
    ftriple[1]<<-lfmid; ftriple[3]<<-rfmid;
   }


 # begin the actual work by "surrounding"
 
  bizarretry <- surround() # pure calls are illegal. Clearly I need to do something else here... # I've <- to <<-'ed
 
 # The logic of the following block is to start with the existing triple of points. Decide which half to fill in by which end has highest f.
 # If that half has higher, then go that way (i.e. set the new triple to the left hand side, and you're done). If that half has lower, then
 # fill in the other half, and if it has highest f, then go that way. Otherwise, you're left with pulling both ends in to make a new triple
 # centred in the same place but narrower.
 
 if (xtriple[3]<=xmin) {
    xout<-xtriple[3]; fout<-ftriple[3]; edge<-1 } else   # The [3] is justified as the user chooses xmin and this is below it...
  {while (xtriple[3]-xtriple[1]>xerr) {                 #  level 1 DO */
    if (ftriple[1]>ftriple[3]) {                     #  level 2 DO */
      lxmid<-sqrt(xtriple[1]*xtriple[2]);
      lfmid<-getlik(lxmid);
      if (! try(xtriple,ftriple,lxmid,lfmid,1,2)) { # level 3 DO */
        rxmid<-sqrt(xtriple[2]*xtriple[3]);
        rfmid<-getlik(rxmid);
        if (! try(xtriple, ftriple,rxmid,rfmid,2,3)) {
           centre(xtriple,ftriple,lxmid,lfmid,rxmid,rfmid)}
        }                                           # level 3 END */
     }    else    # level 2 END */  BEFORE the else
          {                                          # level 2 DO  */
      rxmid<-sqrt(xtriple[2]*xtriple[3]);
      rfmid<-getlik(rxmid);
      if (! try(xtriple,ftriple,rxmid,rfmid,2,3)) { # level 3 DO  */
        lxmid<-sqrt(xtriple[1]*xtriple[2]);
        lfmid<-getlik(lxmid);
        if (! try(xtriple, ftriple,lxmid,lfmid,1,2)) {
          centre(xtriple,ftriple,lxmid,lfmid,rxmid,rfmid) }
        }                                             # level 3 END  */
      }                                               # level 2 END  */
  }                                                   # level 1 END  */
  xout<-xtriple[2];  
  fout<-ftriple[2];   
  edge<-0;  
  }
  
  return(list(rho=xout,lik=fout,edge=edge))

  }  # end of findmax 

# getlik calculates the likelihood concentrated for rho for a given rho, and is used to find the maximum likelihood
#  estimate for rho.

getlik<-function(rho) {    

paperv<- 1 - hmatmarca^rho  # This doesn't affect the value of paperv in the outer environment -- just within this call of getlik
   
 for (ii in 1:nontopnodes) {
  if (spn[ii]) {
  desct<-getdesc(ii,specsetactive,txp)
  papercvec[ii]<-getsigm(desct,paperv)-(1-heights[txp[ii]]**rho) 
  } }
  
  ## TEMPORARY because in a heavily-used section: indeed, I've abandoned even the following line, so persuaded am I all is well here...
  ## if (!all(papercvec>=0))  { print("some of papercvec is negative"); print(papercvec);print(specsetactive);print(txp);browser()}
  #  mybool<-try(all(papercvec>=0))
  #  if (class(mybool)=="try-error" || is.na(mybool) || !mybool) {cat(paste0(" getlik says papercvec has problems, possibly negative values")); browser()}

#  /* note that the heights are assumed scaled so the root is at 1 */
#  /* note also we need to use on[] to get the original name of the
#     node, as txp delivers the number in the reduced phylogeny that
#     includes only the included nodes (though all species) */

lik<-0;

sgliml2<-0;
for (ii in 1:nontopnodes) {
 if (spn[ii])  sgliml2<-sgliml2+log(papercvec[ii])
 }

gliml3<- c(rep(0,nnodes-nspec));
for (ii in 1:nontopnodes) {
 if (spn[ii])  gliml3[txp[ii]-nspec]<-gliml3[txp[ii]-nspec]+1/(papercvec[ii]**2)
 }
 
rhoterm<-(-sgliml2-sum(log(gliml3))/2)/2;

#  /* Note we assume in what follows that we are controlling for something -- and that
#     the constant will be there, so this is OK */
invv<-VCMginv(paperv) 
ssterm<-t(y) %*% invv %*% (diag(nspec)-designx %*% VCMginv(t(designx) %*% invv %*% designx) %*% t(designx) %*% invv) %*% y;

lik<-rhoterm - (nontopactivenodes-nhinodes)*log(ssterm)/2;

return(lik)

 } # finish getlik;

# =================== End of findmax and getlik ====================================
#=================== Now define merge.formulae.ag and merge.formulae.test.ag ==========


merge.formulae.ag <- function(form1, form2, ...){

# adapted and retyped and amended from merge.formula of
# http://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/


         # put character strings of responses into LHS1 and LHS2
         LHS1<-deparse(form1[[2]])
         LHS2<-deparse(form2[[2]])
	if(LHS2 !='.') stop('The response side of an update-style specification of the test terms must be "."')

         # put character strings of RHSs into RHS1 and RHS2
         RHS2<-strsplit(deparse(form1[[3]], "\\+"))[[1]]
         RHS2<-strsplit(deparse(form2[[3]], "\\+"))[[1]]
         
	# here, assume the environment bit is still right.
	
	out<-update(old=form1,new=form2,evaluate=FALSE)
	
	environment(out)<-parent.frame()

	return(out)
	
}  # finish merge.formulae.ag



merge.formula.test.ag <- function(form1, form2, ...){

# adapted, amended and retypes from merge.formula of
# http://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/

         # Get response1 into character string into LHS1
         	LHS1 <- deparse(form1[[2]])
	
	# put character strings for RHSs into RHS1 and RHS2
	
	RHS1<-strsplit(deparse(form1[[3]], "\\+"))[[1]]
	RHS2<-strsplit(deparse(form2[[3]], "\\+"))[[1]]
	ndots<-length(grepRaw(".", RHS2, all=TRUE, fixed=TRUE))  ## grepRaw put in on 2014_01_31
	if (ndots!=0) stop ('There must be no "." term in an "adding terms" style specification of the test terms.')

         # Put together the separate sides
	RHS <- c(RHS1, RHS2)
	LHS <- LHS1

	# put the two sides together with the  
	# reformulate function
	# use "reformulate" to turn the sides into a single formula
	
	out <- reformulate(RHS, LHS)
	
	# set the environment of the formula
	
	environment(out)<-parent.frame()

          return(out)
          
}


# =============== The body of inphyreg ====================

 contm<-"none"
  ret<--1
 if (is.character(intest)) {outtest<-as.formula(paste(deparse(cont[[2]]),as.character(cont[[1]]),deparse(cont[[3]]),"+",intest)); contm<-"char"} 
 else
   if (intest[[1]]=='~') if (length(intest)==3) {outtest<-merge.formulae.ag(cont, intest);contm<-"update"} else
      if (length(intest)==2) {outtest<-merge.formula.test.ag(cont, intest); contm<-"model"};
  if (contm=="none") {ret<-1 ; stop('The test terms were not properly specified. Please consult the documentation.')};
  
  nspec<-nrow(dataframe)
  nomspuse<-nspec-sum(insubset);  # How many species are omitted by user choice
    
  specdatax<-model.frame(cont,dataframe,na.action=na.omit)
  specsetactive<-intersect(attr(specdatax,"row.names"),subset(1:nspec,insubset==1))
  
  spu<-c(rep(0,nspec))
  spu[specsetactive]<-1
  nspecactive<-sum(spu)
  
  nommiss<-sum(insubset-spu) #  How many of these permitted by the user are then omitted for missing values
 
  y<-specdatax[[1]]
  designx<-model.matrix(cont,specdatax)

  specdataxz<-model.frame(outtest,dataframe,na.action=na.omit)

  designxz<-model.matrix(outtest,specdataxz)
  
  ## I now take a moment to check if the control model has zero RSS or the control+test model has zero RSS. This
  ##  can be tested non-phylogenetically as that property is in theory the same whatever the vcm.
  
  test1<-try(lm(cont,dataframe),silent=TRUE)
  if (class(test1)=="try-error") { print("lm on control model fails, about to quit"); browser(); stop()} else
    if (sum(test1$residuals^2)<=tolerance) stop("The control model fits perfectly, so no test possible: nothing left for the test variables to explain.")
  test1<-try(lm(outtest,dataframe),silent=TRUE)
  if (class(test1)=="try-error") { print("lm on control and test model fails, about to quit"); browser(); stop()} else
    if (sum(test1$residuals^2)<=tolerance) stop("The control+test model fits perfectly, so no test possible: perhaps F=infinity.")
  

  myphyls<-getphytxp(phyvar,taxmat,spu)
  phy<-myphyls$phy
  txp<-myphyls$txp
  on<-myphyls$on
  
  # Check all these at some stage for usage and try to reduce notational clutter
  nphyhinodes<-length(phy)-nspec+1
  nontopnodes<-length(txp)
  nnodes<-nontopnodes+1
  nhinodes<-nnodes-nspec
  nshnodes<-nhinodes
  nontopactivenodes<-nontopnodes-nspec+nspecactive;   
    
  specset<-c(1:nspec)
  nontopset<-c(1:nontopnodes)

  nontopsetactive<-union(specsetactive,(nspec+1):nontopnodes)

  spn<-c(spu,c(rep(1,nhinodes)))
  
  # Now for the heights, and the mode depends on arguments
  
   if (is.null(inheights)) hmode<-2 else if (is.null(taxmat)) hmode<-1 else hmode<-3
    
  # for convenience while levhts is delivered as a result
  levhts<- -1
  
#    print("just about to do hmode/levhts")
 #    browser()
     
  if (hmode==2) {heights<-pathlenfig2(nspec,spu,txp); hreport<-paste('Path lengths derived from "Figure 2" default method');levhts<-3} else {
    if (hmode==1)  {if (length(inheights)==length(phy)+1) levhts<-0} else
    if (hmode==3) {if (length(inheights)==ncol(taxmat)) levhts<-1 else if (length(inheights)==length(phy)+1) levhts<-0}
    
  if (levhts==-1L) {
    fop<-paste0("NOTE: inappropriate heights were supplied.\n")
    fop<-paste0(fop,"There must be one for each node including the root, or one for each taxonomic vector.\n")
    fop<-paste0(fop,"length(phy) = ", length(phy)," and length(inheights) = ",length(inheights),"\n")
    fop<-paste0(fop,"hmode = ",hmode,"\n")
    fop<-paste0(fop,('  Correct the weights, or specify "heightsdata=NULL" to use default "Figure 2" heights\n'))
    cat(fop)
    stop("Irrecoverable error.")}
    
    
 #   {if (nrow(inheights)==ncol(taxmat)) levhts<-1} else if (nrow(inheights)==nrow(phy)+1) levhts<-0 else levhts<-2
 
     heights<-rep(-1,length(txp)+1)
    if (hmode==1) {heights<-inheights[on]; hreport<-"Heights of each node were supplied separately"} else {
      if (levhts==1) {inheights <-c(0,inheights); heights<-inheights[1+levels[on]]; hreport<-"Heights of taxonomic levels was supplied"} else {
        if (levhts==0) {heights<-inheights[on]; hreport<-"Heights of each node were supplied separately" } else {
         stop('A heights dataset was supplied but the number of rows did not mach the number of taxonomic levels or total number of nodes')}}}
         
     }  # end of the hmode, levhts and heights-setting code.
     
     
## The above code reflects the possible combinations of hmode and levhts. If hmode is 2, then Fig 2 is used as no heights are
##  supplied. If hmode is 1, then there is no taxmat. So the only good possibility is that the supplied heights are one per node
##  If hmode is 3, then there is taxmat. Now the user may supply one height per taxonomic level OR one height per node

 heights<-scalehts(heights);
 
 hmatmarca<-array(0,c(nspec,nspec)); # to pass into functions. Its NOW the *un-rho-modulated height of* the most recent 
                                                                           # common ancestor function as a matrix. And default changed to 0 so powers work OK
 for (ii in 1:nspec) if (spu[ii]) for (jj in 1:nspec) if (spu[jj]) hmatmarca[ii,jj]<-heights[mrca(ii,jj,txp)] 

 papers<-array(0, c(nontopnodes,nhinodes)) # papers is S from the paper, and has a 1 if a node has another as a descendant
 for (ii in (nspec+1):nnodes) papers[getdesc(ii,setdiff(nontopsetactive,ii),txp),ii-nspec]<-1

 # This prepares for getlik
 paperv<- array(0,c(nspec,nspec))

paperprel<- array(0,c(nnodes,nspec))
paperl <-  array(0,c(nontopnodes,nspec))
papercvec<-c(rep(0,nontopnodes))

edge<--1;    # This block gets the rho for control model, with its likelihood; fitting rho unless user specifies a value
 if(rho<=0) {
   bestrho<-1; bestlik<-79;
   foundmax<- findmax(lorho,hirho,minrho,errrho) 
   bestrho<-foundmax$rho
   bestlik<-foundmax$lik
   edge<-foundmax$edge
       } else {
   bestrho<-rho
   bestlik<-getlik(bestrho)
   edge<-2
 }

onesv<-array(1,c(nspec,1))
paperv<- 1 - hmatmarca^bestrho
invv<-VCMginv(paperv)

for (ii in 1:nnodes) {  
 if (spn[ii]) {
  desci<-getdesc(ii,specsetactive,txp)
  lci<-getlc(desci,paperv)
 paperprel[ii,desci]<-t(lci)
 } };

for (ii in 1:nontopnodes) {  
 if (spn[ii]) paperl[ii,]<-paperprel[ii,]-paperprel[txp[ii],] 
 }


 for (ii in 1:nontopnodes) {
  if (spn[ii]) {
  desct<-getdesc(ii,specsetactive,txp)
  papercvec[ii]<-getsigm(desct,paperv)-(1-heights[txp[ii]]^bestrho) 
  } }

paperc<-diag(papercvec)

regweights<-rep(0,nontopnodes)
for (ii in 1:nontopnodes)  if (spn[ii])  regweights[ii]<-1/papercvec[ii];

longy<-paperl%*%y;    # Creates the long dataset --- note that dummy variables now have non 0,1 values
longdesignx<-data.matrix(paperl%*%designx)
longdesignx<-data.matrix(data.matrix(longdesignx)[,-1])  # Idea is to drop the first column of all zeroes
longdesignxz<-data.matrix(paperl%*%designxz)
longdesignxz<-data.matrix(data.matrix(longdesignxz)[,-1])  # Idea is to drop the first column of all zeroes

invc<-array(0,c(nontopnodes,nontopnodes)); for (ii in 1:nontopnodes) if (papercvec[[ii]]>0L) invc[[ii,ii]]<-1/papercvec[[ii]]
if (ncol(longdesignx)>=1) {    # Is there anything in the control model?
 papere<-(diag(nontopnodes)-as.matrix(longdesignx)%*%VCMginv(t(as.matrix(longdesignx))%*%invc%*%as.matrix(longdesignx))%*%t(as.matrix(longdesignx))%*%invc)%*%longy
}  else {
  papere<-longy} #  /* We are only controlling for the intercept */

# Just moved these long regressions up a few lines on 2014_01_31, so I can get a good value for longxems. I *should* be
#  doing the long before the short.
mylmlongx<-lm.wfit(x=longdesignx, y=longy, w=regweights)
mylmlongxz<-lm.wfit(x=longdesignxz, y=longy, w=regweights)

# longrss<-t(papere)%*%invc%*%papere;
longrss<-sum(mylmlongx$weights*mylmlongx$residuals^2)   ## Note that regweights was wrong whenever it has zero values, as lm.wfit omits those datapoints
longrdf<-mylmlongx$df.residual

if (longrss<=tolerance) {print("No variation in long regression, so test is impossible. About to exit.");browser(); stop()} # Fixes y=0 case but also more widely

papertau<-array(0,c(nontopnodes))
ecm1e<-array(0,c(nhinodes))
paperg<-array(0,c(nshnodes,nontopnodes))
for (ii in 1:nontopnodes) if (spn[ii]) ecm1e[txp[ii]-nspec]<-ecm1e[txp[ii]-nspec]+papere[ii]^2/paperc[ii,ii]

 ecm1e<-sqrt(ecm1e);

longxems<-longrss/longrdf; # /* Doesn't this make rank(X) matter not just numx? */
                            #  /* and it should surely be nspecactive */
      # CONSIDER fixing this mistake, if it is one, after checking I get the right numbers out the end
      #  Maybe its used only for the tolerance calculation
      # I've replaced it with the surely wholly rational appeal to lm.wfit's output (mulmlongx calculation just moved to a few lines up)
      #  and it looks as though my shregweights calculation is unchanged if regweights is multiplied by an arbitrary positive constant

# Then we get shregweights, which will be 1 or 0. Its done by an OR on all the daughter nodes having a residual that's
#  large enough according to the tolerance. So long as one daughter does, the higher node is included.

shregweights<-rep(0,nhinodes)
for (ii in 1:nontopnodes){ 
 if (spn[ii]) shregweights[txp[ii]-nspec]<- (shregweights[txp[ii]-nspec] ||  abs(papere[ii]*sqrt(regweights[ii]))>=sqrt(longxems)*tolerance)}  

#  /* I calculated the missing nodes above, and I use them and record them here */ # though with an extra "0" that always later needs
#   accounting for.
missingnodes<-0;
 for (ii in 1:nontopnodes) {
  if (spn[ii]) {
    if (shregweights[txp[ii]-nspec]) {
     papertau[ii]<-papere[ii]/ecm1e[txp[ii]-nspec]
     paperg[txp[ii]-nspec,ii]<-papertau[ii]      
    }    else {
     missingnodes<-union(missingnodes,txp[ii])
     papertau[ii]<-0
     paperg[txp[ii]-nspec,ii]<-0
    }
 }}
  
 papergcm1l<-paperg%*%invc%*%paperl;  # So this gets the paper's GC^{-1}L, I think -- the contrasts to get from long to short

wtdmeansxz<-t(onesv)%*%invv%*%cbind(y,designxz)/(t(onesv)%*%invv%*%onesv)[[1,1]]

 shorty<-papergcm1l %*% y
 shortdesignx<-data.matrix(papergcm1l%*%(data.matrix(designx)[,-1]))
 shortdesignxz<-data.matrix(papergcm1l%*%(data.matrix(designxz)[,-1]))  ## from here search for papercgm1l and designx
 
 shornode<-on[(nspec+1):nnodes]

## Note I had huge time-consuming trouble because the weights were arrays. They need to begin life as
##  just rep(0,length), and this allows them to be reused appropriately in the calculation of x * sqrt(w) inside lm.wfit.

# Note also that lm.wfit omits datapoints with zero weight. Thus if shregweights has zeroes, mylmshortx$weights will still be all 1s, and
#   shorter than shregweights. Thus I can't mix shregweights with the outputs of lm.wfit. I guess the same is true for the long ones, and 
#   will need looking after.
#
## That turns out to be too simple. lm.wfit omits datapoints in a null model that have y=0. It does not omit datapoints in a non-null model
##  that have y=0. I think this is the whole principle. If so, I need to expand the output of null models to the same length, to do my RCO/RCT
##  calculations. It also means that users might have a problem interpreting mylmshortx -- and I'm not sure I'll do much about that now.
##  I wouldn't like to be confused myself, however...

mylmshortx<-lm.wfit(x=shortdesignx, y=shorty, w=shregweights)
mylmshortxz<-lm.wfit(x=shortdesignxz, y=shorty, w=shregweights)

### Well, now it's just output and a few fancy bits. The work is all done.

## First, some calculations of important numbers

 dflx<-mylmlongx$rank
 dflxz<-mylmlongxz$rank
 dfsx<-mylmshortx$rank
 dfsxz<-mylmshortxz$rank

 dfxlost=dflx-dfsx;
 dfxzlost=dflxz-dfsxz;

testdendf<-mylmshortxz$df.residual - addDF
testnumdf<-dfsxz-dfsx;
testdenss<-sum(mylmshortxz$weights*mylmshortxz$residuals^2)    # What about weights? I've added them now
testtotss<-sum(mylmshortx$weights*mylmshortx$residuals^2)          #  ditto (2014_02_04_1455)
testnumss<-testtotss-testdenss

testlongdf=dflxz-dflx;
testlost=testlongdf-testnumdf;

shorttotdf=sum(shregweights);
shortcondf=dfsx;

hilostomsp<-nphyhinodes-nhinodes;
hilostvar <- length(missingnodes)-1;

testf=(testnumss/testnumdf)/(testdenss/testdendf);
poff=1-pf(testf,testnumdf,testdendf);

det<-list()
det$poff<-poff; det$testf<-testf; det$testnumdf<-testnumdf; det$testdendf<-testdendf
det$nspec<-nspec; det$nommiss<-nommiss; det$nspecactive<-nspecactive
det$nomspuse<-nomspuse; 

det$nphyhinodes<-nphyhinodes; det$nhinodes<-nhinodes; det$hilostomsp<-hilostomsp;
det$hilostvar<-hilostvar; det$shorttotdf<-shorttotdf

det$dflx<-dflx; det$dfxlost<-dfxlost; det$dfsx<-dfsx; det$testlongdf<-testlongdf; det$testlost<-testlost
det$shortcondf<-shortcondf; det$addDF<-addDF

det$rho<-bestrho; det$lik<-bestlik; det$edge<-edge

det$missingnodes<-missingnodes


## the plot stuff. But RCO may need expanding if shorty=0 and designx is null. The code below is inefficient but seems to work:
##  I first expand RCO if necessary, then calculate DR, multiplying by WSH (zero for omitted nodes), and then shrink both
##  RCT and DR back down. I could shrink RCT in the first place if appropriate, omit WSH weights, and Bob would be my
##  uncle...

 WSH<-shregweights  # Note this should equal lmshortxz$weights but not necessarily lmshortx$weights (which are always 1)
 NN<-shornode[ WSH==1L]
  
  if (length(mylmshortxz$residuals)>length(mylmshortx$residuals)) {
       RCO<-rep(0,length(mylmshortxz$residuals))
      RCO[(shregweights!=0)]<-mylmshortx$residuals/(mylmshortx$df.residual-det$addDF) } else
    RCO<-mylmshortx$residuals/(mylmshortx$df.residual-det$addDF)
  RCT<-mylmshortxz$residuals/(mylmshortxz$df.residual-det$addDF)
  if (length(RCO)!=length(RCT)) {cat(paste0("RCO and RCT of different lengths! Alarm!! About to browse...")); browser()}
  DR<- WSH*( RCO^2- RCT^2)
  RCT<- RCT[( WSH==1L)]
  DR<- DR[( WSH==1L)]
 
 det$influence<-data.frame(NN=NN, RCT=RCT/sqrt(mean(RCT^2)), DR=100*DR/sum(DR))
 ## DR so that the total change is 100, RCT so its standardised so its absolute value is like the absolute value of a standardised variable
 
 
# ========================   det should now contain all the information for the stuff *I* have to format

freq<-list(); preq<-list()

preq$parmx<-parmx; preq$parmxz<-parmxz; preq$means<-means; preq$opfunccall<-opfunccall

if (!missing(oppf)) {if (oppf>0L) freq$oppf<-oppf} 
if (!missing(opdf)) {if (opdf>0L) freq$opdf<-opdf} 
if (!missing(oprho)) {if (oprho>0L) freq$oprho<-oprho} 
if (!missing(dfwarning)) { if (dfwarning>0L) freq$dfwarning<-dfwarning}
if (!missing(plot)) { if (plot>0L) freq$plot<-plot} 

#  freq has the flags for the things I may have to format the output for.

# Now I print the information I have to format
  if (sum(as.numeric(freq))>0L) opsi(details=det, requests=freq, H0model=cont, HAmodel=outtest, originalIDs=on)

# Now I print the requested information that I use print() for
if (parmx) {cat("\n"); print(coef(mylmlongx))}
if (parmxz) {cat(paste("\n"));print(coef(mylmlongxz))}
if(opfunccall) {cat(paste("\n\nCall of phyreg:\n\n",sep="")); print(funccall)}
if(means) {cat(paste("\n")); print(wtdmeansxz)}

# anova(mylmshortx,mylmshortxz,test="F")  ## I may want to imitate this (or alternatively mimic lm to be able to get from lm.fit to lm, thence to anova

 #  the always included outputs go in first, then the optional ones conditionally
 ## Curiously, printing ipr afterwards shows that HAmodel comes with an environment, which no other element does. A problem?
 ipr<-list(H0model=cont, HAmodel=outtest, spu=spu,nomspuse=nomspuse,  longrss=longrss, shornode=shornode, details=det, means=wtdmeansxz, parmx=coef(mylmlongx), parmxz=coef(mylmlongxz), funccall=funccall, fullphy=phy, usedphy=txp, originalIDs=on)
 
 #ipr$rho<-list()
 #ipr$rho$bestrho<-bestrho
 #ipr$rho$bestloglik<-bestlik
 #ipr$rho$edge<-c("Error state","Internal maximum of likelihood","Minimum value of rho","Set by user")[2+edge]
 
 if(lmshortx | lmshortxz | lmlongx | lmlongxz)  ipr$lm<-list()
 if (lmshortx) {ipr$lmshortx<-mylmshortx}
 if (lmshortxz) {ipr$lmshortxz<-mylmshortxz}
 if (lmlongx) {ipr$lmlongx<-mylmlongx}
 if (lmlongxz) {ipr$lmlongxz<-mylmlongxz}

 if (hinput) {ipr$hinput<-heights}
 
 if (linputs) {ipr$linputs<-list(); ipr$linputs$y<-longy ; ipr$linputs$design<-longdesignxz; ipr$linputs$w<-regweights}
 if (sinputs) {ipr$sinputs<-list(); ipr$sinputs$y<-shorty;  ipr$sinputs$design<-shortdesignxz; ipr$sinputs$w<-shregweights}
 if (paper) { ipr$paper<-list()
     ipr$paper$c<-papercvec; ipr$paper$l<-paperl; ipr$paper$gcm1l<-papergcm1l; 
     ipr$paper$e<-papere; ipr$paper$g<-paperg; ipr$paper$tau<-papertau; ipr$paper$s<-papers
     }
     
     
   # And now for the PGLS FVs.
 
 FVx<-wtdmeansxz[1]
 if (ncol(longdesignx)>=1L) for (kk in 1:length(mylmlongx$coefficients)) FVx<-FVx+mylmlongx$coefficients[kk]*(designx[,kk+1]-wtdmeansxz[2+kk])
 ipr$pglsFVx<-FVx
 
 FVxz<-wtdmeansxz[1]
 if (ncol(longdesignxz)>=1L)  for (kk in 1:length(mylmlongxz$coefficients)) FVxz<-FVxz+mylmlongxz$coefficients[kk]*(designxz[,kk+1]-wtdmeansxz[2+kk])
 ipr$pglsFVxz<-FVxz

     
   ## Now to deal with sdiagnostics and plot  ### Now done along other output using opsi, above, with automatic storage of plotted variables
   
   
#  WSH<-mylmshortx$weights
 # NN<-shornode[ WSH==1L]
#  RCO<-mylmshortx$residuals/(mylmshortx$df.residual-ipr$det$addDF)
#  RCT<-mylmshortxz$residuals/(mylmshortxz$df.residual-ipr$det$addDF)
#  DR<- WSH*( RCO^2- RCT^2)
#  RCT<- RCT[( WSH==1L)]
#  DR<- DR[( WSH==1L)]
 
#  if (sdiagnostics) {
#  ipr$sdiagnostics<-list()
#  ipr$sdiagnostics$NN<-NN
#  ipr$sdiagnostics$RCO<RCO
#  ipr$sdiagnostics$RCT<-RCT
#  ipr$sdiagnostics$WSH<-WSH
#  ipr$sdiagnostics$DR<-DR
#   }
  
#   if (plot) plot(RCT,DR,main="Diagnostics of single contrasts",
#                                        xlab="Standardised residuals in control+test model",
#                                        ylab="Improvement in squared residual from adding test variables",
#                                        sub="How each higher node's single contrast contributes to error") 
 
   ## and we're done here...
   
return(ipr)
    
  }  # finish inphyreg
  
getphytxp<-function(phyvar,taxmat,spu) {  # called by inphyreg

# First get the phylogeny for all nspec species and take a number

mmp<-"No taxmat"      # Is this used anywhere, mmp?
  if (is.null(phyvar) && is.null(taxmat)) stop("No phylogeny available -- specify either phyvar or taxmat")
 if (!is.null(phyvar)) {phy<-phyvar;phyreport<-paste("Phylogeny read in from variable ",phyvar)}
 if (!is.null(taxmat)) {mmp<-makephy(taxmat); phy<-mmp$phy;phyreport<-paste("Phylogeny created from taxonomic variables in ",taxmat)}

 nspec<-min(phy)-1;

 # Then deal with missing species and get txp
 
 myuw<-uw(phy,spu,nspec)
 
 return(list(phy=phy,txp=myuw$txp,on=myuw$on,mmp=mmp))
 
 } # end of getphytxp
 
 # getdesc lists all the descendants of innode, according to phy, that are also in mask
getdesc<-
  function(innode,mask,phy) {
  ing<-c(innode)
  if(innode>1) {  # This IF resolves the problem noted just below
    for (i in seq(innode-1,1,-1) ) { # This needs avoiding in R -- SAS just didn't do the loop when innode-1=0 -- R fails
      par<-phy[i]
      if(phy[i] %in% ing) ing<-union(i,ing)  # I've replaced c(i) with just plain i in the union -- it seems to work
    }}             # I've just changed ing,i to i,ing in the union, to hope to get things in the right order
   ing<-intersect(ing,mask)
   return(ing)
    }

# getlc takes a matrix paperv and a subset of nodes desc (for descendants of a given node), and works out the 
#  linear contrast for taking the weighted mean for the given node. (I think!)
getlc<-function(desc, paperv) {
   v<-as.matrix(paperv[desc,desc])   # as.matrix added so that R will deal with the singleton case 
   ones<-rep(1, dim(v)[1])
   sigma<-1/sum(solve(v, ones))
   lc<-sigma*solve(v, ones)
   return(lc)
   } # finish getlc
   
# getsigm takes a matrix paperv and a subset of nodes desc (for descendants of a given node), and works out the 
#  sampling variance of the optimal weighted contrast for the mean for the given node. (I think!)
getsigm<-function(desc, paperv) {
   v<-as.matrix(paperv[desc,desc])
   ones<-rep(1, dim(v)[1])
   sigma<-1/sum(solve(v, ones))
   return(sigma)
   }  # finish getsigm

### On 24 January 2014 I've brought in newuw from uw_workings, and renamed it uw, at the same time renaming
###  the old uw into olduw. But the old one is so horrible, I should simply delete it, as its a patchwork of
###  illogicalities. Only sentiment and excessive caution prevent me deleting it now altogether...

## No, uw is now defined at top level

   
pathlenfig2<- function(nspec,spu,txp) {

  b<-rep(0,length(txp)+1)
  
  nspecactive<-0
  for (ii in 1:nspec) {b[ii]<-spu[ii]; nspecactive<-nspecactive+1}
    
  for (ii in 1:length(txp)) if (txp[ii]>0) b[txp[ii]]<-b[txp[ii]]+b[ii] 
  # b now has the number of species under each node, with 0 for omitted nodes
  b<-b-1 
  # b now has the number of species MINUS ONE under each node, with -1 for omitted nodes
  #### for (ii in 1:nspec) b[ii]<-b[nspec]+1-spu[ii] 
  # I think this should have b[nspec] replaced with b[ii]. Point of line is to put 0 for -1 for omitted?
  #  Works for the included species (i.e. leaves them at zero) provided b[nspec]==0 i.e. final species is included
  # Wouldn't work if it wasn't, for then b[nspec]==-1, and included species become -1, excluded become 0
  # As this has support from the GLIM macro csl_, I'm implementing it., thus:
  for (ii in 1:nspec) b[ii]<-b[ii]+1-spu[ii] 
  
  b<-b/(nspecactive-1) 
  # This relativises the heights to 1
  
return(b)

  } # finish pathlenfig2;

#========= Note this is a new version of makephy, replacing the old one on 2014_01_31. See bottom of definition

#                         for more details.
makephy<-function(taxmat) {

# This (previously makephyc) is a version that works with character variables for the elements of the matrix

nspec<-nrow(taxmat)
firstrow<-seq(1,nspec)

taxmat<-cbind(species=firstrow,taxmat)

# print(taxmat)

a<-array(0,c(nspec,nspec))
c<-seq(1,nspec);

nvectax<-ncol(taxmat);

for (jj in 1:(nspec-1)) {
  for (ii in (1+jj):nspec) { # print("ii then jj");print(ii); print(jj);
    mm<-nvectax; while (taxmat[[ii,mm]] == taxmat[[jj,mm]]) mm<-mm-1
    a[ii,jj]<-mm
   }}

nextnodenumber<-nspec+1;
orignodevec<-c(-1)
origparvec<-c(-1)
origlevvec<-c(-1)
for (vecdiff in 1:nvectax) {
 d<-rep(0,nspec) 
 for  (jj in 1:(nspec-1)) {
  notstarted<-1
  sistersofjj<--1
  if (d[jj]==0) {
   for (ii in (jj+1):nspec) {
    if (a[[ii,jj]]<vecdiff) {
     if (sistersofjj[1]==-1) sistersofjj<-ii else sistersofjj<-c(sistersofjj, ii)}
    if (d[ii]==0) {
     if (notstarted) {
      if (a[[ii,jj]]==vecdiff) {
       notstarted<-0
       curnodenumber<-nextnodenumber
       nextnodenumber<-nextnodenumber+1
       daughters<-c(c[jj], c[ii])
       c[jj]<-curnodenumber
       c[ii]<-curnodenumber
       d[ii]<-jj
       } # end; /* a[]=vecdiff */
     } # end; /* notstarted THEN */
     else { # do; /* started */
      if (a[[ii,jj]]==vecdiff) {
       daughternumber<-c[ii]
       daughters<-c(daughters, daughternumber)
       c[ii]<-curnodenumber
       d[ii]<-jj;
      } # end;/* the if a[ii,jj]=vecdiff */
     } # end; /* the else for notstarted*/
    } # end;   /* the d[ii]=0 then */
   } # end;    /* do ii */
   if (notstarted==0) {
    if (sistersofjj[1] != -1) {
      c[sistersofjj]<-c[jj]
      d[sistersofjj]<-jj
     }
     daughters<-unique(sort(daughters))
     parvec<-rep(curnodenumber, length(daughters))
     levvec<- rep(vecdiff, length(daughters)) 
     orignodevec<-c(orignodevec, daughters)
     origparvec<-c(origparvec, parvec)
     origlevvec<-c(origlevvec, levvec)
    } # end; # of if (notstarted==0), which seems right
   } # end;     /* do jj */ -- currently and wrongly if d[jj]==0... ## I must have swapped these in the SAS
  } # end;     /* if d[jj]=0 */ -- currently and wrongly do jj...    ## comments, and so wasted my time doing R...
} # end;      /* do vecdiff */

orignodevec<-orignodevec[-1];
origparvec<-origparvec[-1];
origlevvec<-origlevvec[-1];

phy<-rep(-1,length(origparvec))
levels<-rep(-1, length(origparvec))
for (ii in 1:length(phy)) {
  phy[orignodevec[ii]]<-origparvec[ii]
  levels[orignodevec[ii]]<-origlevvec[ii]
  }

# levels gives the level of the top of the segment of which ii is the bottom
# I'm going to convert it now before adding the top level

blevels<-levels
levels[phy]<-levels

nextlevel<-max(levels)+1
levels<-c(levels, nextlevel) # /* The root requires a level too. */

nextlevel<-max(blevels)+1
blevels<-c(blevels, nextlevel) # /* The root requires a level too. */


return(list(phy=phy,levels=blevels))

} # finish makephy (formerly makephyc);
## That was the "blevels" version from make_random.tex, which here I have made native by 
##  setting levels=blevels, and omitting the internal levels, which is just wrong.
## The version in make_random still at the moment has levels and blevels, because it's
##  embedded in another function that uses blevels.


 scalehts<-function( vec ) {
    return((vec-min(vec))/(vec[length(vec)]-min(vec)))
    } #  finish scalehts;


## ========================================================================
##========================================================================

#  NOW the body of phyreg begins...

funccall<-match.call()

cfile<-paste0(tempdir(),"/","curparms")

save_curparms<-function() {
  file.create(cfile)
  save(curparms,file=cfile)
}

get_curparms<-function(reset) {
 if (file.exists(cfile)) try(load(cfile),silent=TRUE)
 if (!exists("curparms") || reset){
   curparms<-list(control=NA, test=NA, subset=NULL,data=NA, phydata=NULL, heightsdata=NULL, addDF=0, rho=-1, lorho=0.3, hirho=0.6, errrho=0.02, 
      minrho=0.0001, tolerance=0.000001, oppf=5, opdf=0, parmx=0, parmxz=0, opfunccall=0, taxmatrix=NULL, linputs=FALSE, sinputs=FALSE, means=FALSE,
       lmshortx=FALSE, lmshortxz=FALSE, lmlongx=FALSE, lmlongxz=FALSE, hinput=FALSE, paper=FALSE, dfwarning=TRUE, oprho=FALSE, 
       plot=FALSE, reset=FALSE)}
  return(curparms)
}

curparms<-get_curparms(reset)

tempcurparms<-curparms

compargs<-c("control","test")
dataargs<-c("data","phydata","heightsdata","taxmatrix")
optdefargs<-c("subset","addDF","rho","lorho","hirho","errrho","minrho","tolerance","oppf","opdf","dfwarning","parmx","parmxz","opfunccall", "linputs", "sinputs", "means", "lmshortx", "lmshortxz", "lmlongx", "lmlongxz", "hinput", "paper", "oprho",  "plot", "reset")

strcomp<-"if (!((missing(myvar)))) tempcurparms$myvar<-myvar"
#stropt<-"if (!((missing(myvar)))) tempcurparms$myvar<-myvar"
stroptdef<-strcomp

# THEN I don't see that I need the distinction into the different types of parameters. But I guess I could maintain it in case I do later
#  realise why...

  for (ss in compargs) {eval(parse(text=gsub("myvar",ss,fixed=TRUE,x=strcomp)))}
#  for (ss in optargs) {eval(parse(text=gsub("myvar",ss,fixed=TRUE,x=stropt)))}
  for (ss in optdefargs) {eval(parse(text=gsub("myvar",ss,fixed=TRUE,x=stroptdef)))}
  
# Now for my new (2014_02_02) data argument handling, to store them as variable names (as strings) and not as the objects
     fcmf <- match.call(expand.dots = FALSE)
     targets<-c("data","phydata","taxmatrix","heightsdata")
     indx<-which(names(fcmf) %in% targets,arr.ind=TRUE)
     clist<-list(data="",phydata="",taxmatrix="",heightsdata="")
    
### We need this recursive function to get dollar-ed symbols into character variables
dolc<-function(x) {
 if (length(x)==1) {if  (is.symbol(x)) return (paste0(x)) else return("--mistake, l1, not symbol--")} else
  { if (x[[1]]=="$") return (c(dolc(x[[2]]),"$",dolc(x[[3]]))) else return ("-- mistake (residual) --")}}
  
# We give it the argument, it gives us a character vector, which we paste0(...,collapse="") to get into what we want.

#print("just about to do clist/dolc")
#browser()

     if(length(indx)>=1L) for (ii in 1:length(indx))
            if (is.null(fcmf[[indx[[ii]]]])) clist[names(fcmf)[[indx[[ii]]]]]<-"NULL" else clist[names(fcmf)[[indx[[ii]]]]]<-paste0(dolc(fcmf[[indx[[ii]]]]),collapse="")
            
     for (kk in 1:length(clist)) if (clist[[kk]] != "") tempcurparms[names(clist)[kk]]<-clist[kk]
     
 

if (!is.character(tempcurparms$data)) stop("No data supplied -- please use 'data=' argument of phyreg()")
if (!exists(tempcurparms$data)) stop("The data=  argument to phyreg does not exist as a symbol -- check spelling?")

## The || condition here (and the setting to "NULL" rather than NULL of the elements of clist) are needed because while it is possible to 
##  to set a variable to NULL, and to define a list with some NULL elements, I can't find a way to assign an element in a list to NULL
##  But this solution keeps those complications out of inphyreg. (Note the || evaluates only as necessary left to right.)

if(is.null(tempcurparms$data) || tempcurparms$data=="NULL") {dataarg<-NULL; tempcurparms$data<-"NULL"}  else 
  dataarg<-eval(parse(text=tempcurparms$data))
if(is.null(tempcurparms$phydata) || tempcurparms$phydata=="NULL") {phyvararg<-NULL; tempcurparms$phydata<-"NULL"}  else 
  phyvararg<-eval(parse(text=tempcurparms$phydata))
if(is.null(tempcurparms$taxmatrix) || tempcurparms$taxmatrix=="NULL") {taxmatrixarg<-NULL; tempcurparms$taxmatrix<-"NULL"} else 
  taxmatrixarg<-eval(parse(text=tempcurparms$taxmatrix))
if(is.null(tempcurparms$heightsdata) || tempcurparms$heightsdata=="NULL") {heightsdataarg<-NULL; tempcurparms$heightsdata<-"NULL"} else 
  heightsdataarg<-eval(parse(text=tempcurparms$heightsdata))
  
  
if ((is.null(taxmatrixarg)+is.null(phyvararg))!=1) stop('Remove either taxmatrix or phyvar with "parameter=NULL" -- only one is allowed.')
if (!is.null(phyvararg) && !is.null(heightsdataarg)) 
   if (length(phyvararg)+1!=length(heightsdataarg)) stop("Your heightsdata should be one element longer than your phydata")
#if (!is.null(taxmatrixarg) && !is.null(heightsdataarg)) 
 # if (length(heightsdataarg)!=ncol(taxmatrixarg))           ### This has to be done inside because we need to get phy from taxvar to see if the heights
  #  stop("Your heightsdata should one element for each column of your taxmatrix")        ##   are right that way.
  
  if (!is.null(heightsdataarg) && any(is.na(heightsdataarg))) stop("Some of your supplied heights are NA (i.e. missing)")    ### a bad thing...
  if (!is.null(phyvararg) && any(is.na(phyvararg))) stop("Some of your supplied phydata values are NA (i.e. missing)")    ### a bad thing...
  if (!is.null(taxmatrixarg) && any(is.na(taxmatrixarg))) stop("Some of your supplied taxmatrix values are NA (i.e. missing)")    ### a bad thing...
  
    # Now to deal with subset...
  
   if (is.null(tempcurparms$subset))  tempcurparms$subset<-rep(1,nrow(dataarg)) # This is the only non-scalar parameter with a default, 
                                                                                                                                                   #     so it needs separate setting here
    if(any(!is.numeric(tempcurparms$subset))) stop("Your subset variable should not have missing (NA) values")
                                                                          ## Note that the "data variables" are stored by name (and so can be changed on the fly) but subset is stored by value
 
   subset<-NA  ## This is so that within inphyreg, R will use its base function subset, instead of fretting about the parameter of phyreg of the same name.


    curparms<-tempcurparms # The double-act and late assignment ensure that erroneous variable sets don't get into curparms
                                                    #    but we don't save until afterwards anyway.

 # name conversions from phyreg to inphyreg are cont/control, intest/test, phyvar/phydata, taxmat/taxmatrix, inheights/heightsdata 

  ipr<-inphyreg(cont=curparms$control, intest=curparms$test, dataframe=dataarg, insubset=curparms$subset,  phyvar=phyvararg, taxmat=taxmatrixarg, inheights=heightsdataarg, rho=curparms$rho, lorho=curparms$lorho, hirho=curparms$hirho, errrho=curparms$errrho, minrho=curparms$minrho, tolerance=curparms$tolerance, oppf=curparms$oppf, opdf=curparms$opdf, dfwarning=curparms$dfwarning, parmx=curparms$parmx, parmxz=curparms$parmxz, opfunccall=curparms$opfunccall, addDF= curparms$addDF, means=curparms$means, sinputs=curparms$sinputs, linputs=curparms$linputs, lmshortx=curparms$lmshortx, lmshortxz=curparms$lmshortxz, lmlongx=curparms$lmlongx, lmlongxz=curparms$lmlongxz, hinput=curparms$hinput, paper=curparms$paper,oprho=curparms$oprho, plot=curparms$plot )
  
  result<-save_curparms()  ## save only if the call to inphyreg finishes properly (DOCUMENT THIS)

 class(ipr)<-"phyreglm"

 return(invisible(ipr))
 
   }
