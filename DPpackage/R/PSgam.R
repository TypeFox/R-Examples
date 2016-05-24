### PSgam.R                   
### Fit generalized additive linear mixed model with B splines 
###
### Copyright: Alejandro Jara, 2007-2012.
###
### Last modification: 15-09-2010.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author' contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###

"PSgam"<-
function(formula,family=gaussian(),offset,n,prior,mcmc,state,status,
         ngrid=20,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PSgam")


"PSgam.default"<-
function(formula,
         family=gaussian(),
         offset=NULL,
         n=NULL,
         prior,
         mcmc,
         state,
         status, 
         ngrid=50,
         data=sys.frame(sys.parent()),
         na.action=na.fail)
{

       #########################################################################################
       # internal functions
       #########################################################################################

       #################################################
       # Construct the matrix for n-th differences
       #################################################
         ndiff <- function(n, d = 1) 
         {
           if (d == 1)
             {D <- diff(diag(n))}
           else
            {D <- diff(ndiff(n, d - 1))}
           D
         }

       #################################################
       # Truncated p-th power function
       #################################################
         tpower <- function(x, t, p)
         {
              (x - t) ^ p * (x > t)
         }

       #################################################
       # Construct B-spline basis and Penalty matrix.
       #################################################
         b.base <- function(x, xl, xr, ndx, deg, pord, names)
         {
            dx <- (xr - xl) / ndx
            knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
            P <- outer(x, knots, tpower, deg)
            n <- dim(P)[2]
            D <- ndiff(n, deg + 1) / (gamma(deg + 1) * dx ^ deg)
            B <- (-1) ^ (deg + 1) * P %*% t(D)
            p <- dim(B)[2]

            D <- ndiff(p, pord)
            K <- t(D)%*%D

            colnames(B) <- paste(names,seq(1,p),sep="-")
            colnames(K) <- paste(names,seq(1,p),sep="-")
            rownames(K) <- paste(names,seq(1,p),sep="-")
            out<-list(B=B,knots=knots,K=K) 
         }

       #################################################################
       interpret.PSgam <- function (gf)
       #################################################################   
       # interprets a PSgam formula of the generic form:
       #   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
       # and returns:
       # 1. a model formula for the parametric part: pf (and pfok indicating whether it has terms)
       # 2. a list of descriptors for the smooths: smooth.spec
       # 3. a full version of the formulae with all terms expanded in full
       {
          p.env<-environment(gf) # environment of formula
          tf<-terms.formula(gf,specials="ps") # specials attribute indicates which terms are smooth   

          terms<-attr(tf,"term.labels") # labels of the model terms 
          nt<-length(terms) # how many terms?
      
          if (attr(tf,"response")>0)  # start the replacement formulae
          { response<-as.character(attr(tf,"variables")[2])
            pf<-rf<-paste(response,"~",sep="")
          }
          else pf<-rf<-"~"
          sps<-attr(tf,"specials")$ps     # array of indices of smooth terms 

          off<-attr(tf,"offset") # location of offset in formula

          vtab <- attr(tf,"factors") # cross tabulation of vars to terms
          if (length(sps)>0) for (i in 1:length(sps)) 
          {
            ind <- (1:nt)[as.logical(vtab[sps[i],])]
            sps[i] <- ind # the term that smooth relates to
          }

          ns <- length(sps) # number of smooths
          k <- ks <- kp <- 1 # counters for terms in the 2 formulae
          len.sps <- length(sps)
 
          smooth.spec <- list()
          if (nt)
          for (i in 1:nt) # work through all terms
          { 
            if (k<=ns&&((ks<=len.sps&&sps[ks]==i))) # it's a smooth
            { 
              st <- eval(parse(text=terms[i]),envir=p.env)
              if(k>1||kp>1)
              {
                 rf <- paste(rf,"+",st$full.call,sep="") # add to full formula
              }   
              else
              {
                 rf <- paste(rf,st$full.call,sep="")
              }   
              smooth.spec[[k]] <- st
              ks <- ks+1  # counts ps() terms
              k <- k+1    # counts smooth terms 
            } 
            else          # parametric
            { 
              if (kp>1) pf <- paste(pf,"+",terms[i],sep="") # add to parametric formula
              else pf<-paste(pf,terms[i],sep="")
              if (k>1||kp>1) rf <- paste(rf,"+",terms[i],sep="") # add to full formula
              else rf<-paste(rf,terms[i],sep="")
              kp <- kp+1    # counts parametric terms
            }
          }    
      
          if (!is.null(off)) # deal with offset
          { if (kp>1) pf <- paste(pf,"+",sep="")
            if (kp>1||k>1) rf <- paste(rf,"+",sep="")
            pf <- paste(pf,as.character(attr(tf,"variables")[1+off]),sep="")
            rf <- paste(rf,as.character(attr(tf,"variables")[1+off]),sep="")
            kp <- kp+1          
          }
          if (attr(tf,"intercept")==0) 
          {pf <- paste(pf,"-1",sep="");rf<-paste(rf,"-1",sep="");if (kp>1) pfok <- 1 else pfok <- 0}
          else { pfok <- 1;if (kp==1) { pf<-paste(pf,"1"); if (k==1) rf<-paste(rf,"1",sep="");}}
      
          fake.formula<-pf
          if (length(smooth.spec)>0) 
          for (i in 1:length(smooth.spec))
          { nt<-length(smooth.spec[[i]]$term)
            ff1<-paste(smooth.spec[[i]]$term[1:nt],collapse="+")
            fake.formula<-paste(fake.formula,"+",ff1)
          }
          fake.formula<-as.formula(fake.formula,p.env)
          ret<-list(pf=as.formula(pf,p.env),
                    pfok=pfok,
                    psok=k-1,
                    smooth.spec=smooth.spec,
                    full.formula=as.formula(rf,p.env),
                    fake.formula=fake.formula,
                    response=response)
          class(ret)<-"split.PSgam.formula"
          ret
       }

       normalize<-function(a,b)
       {
          a-b     
       }

       #################################################
       # Construct the Z matrix
       #################################################
         z.matrix <- function(object,datas,ngrid)
         {
             nevalps <- object$psok
             nsmooth <- 0
             z <- NULL
             znew <- NULL
             zmean <- NULL
             xreal <- NULL
             pmatrix <- NULL
             pordv <- NULL
             n.smoothers <- NULL
             nr.smoothers <- NULL
             
             count <- 0
             starts <- NULL
             ends <- NULL
             nknots <- NULL
             knots <- list()

             for(i in 1:nevalps)
             {
                 obj <- object$smooth.spec[[i]]
                 dimen <- obj$dim
                 ndx <- obj$bs.dim
                 deg <- obj$s.degree
                 pord <- obj$pord
                 
                 for(j in 1:dimen)
                 {
                    namer <- obj$term[j]
                    n.smoothers <- c(n.smoothers,paste("ps(",namer,")",sep=""))
                    nr.smoothers <- c(nr.smoothers,namer)
                    x <- datas[colnames(datas)==namer][[1]]
                    nrec <- length(x)
                    xl <- min(x)
                    xr <- max(x)
                    xmax <- xr + 0.01 * (xr - xl)
					xmin <- xl - 0.01 * (xr - xl)
                    
                    xgrid <- seq(xmin,xmax,len=ngrid)
                    xreal <- cbind(xreal,xgrid)

                    bbase <- b.base(c(x,xgrid,mean(x)), xmin, xmax, ndx, deg, pord, namer)
                    
                    count <- count+1
                    
                    pordv <- c(pordv,pord) 
                    
                    if(count==1)
                    {
                       starts <- 1
                       ends <- ndx+deg
                       send <- ends
                       pmatrix <- bbase$K
                    }
                    else
                    {
                       starts <- c(starts,send+1)
                       ends <- c(ends,send+ndx+deg)
                       send <- send+ndx+deg
                       pp <- dim(bbase$K)[2]+dim(pmatrix)[2]
                       
                       #cat("pp=",pp,"\n")
                       #cat("pmatrix=",dim(pmatrix)[2],"\n")
                       #cat("ends=",ends[count-1],"\n")
                       #cat("K=",dim(bbase$K)[2],"\n")
                       
                       pmatrixt <- matrix(0, nrow=pp,ncol=pp) 
                       pmatrixt[(1:ends[count-1]),(1:ends[count-1])] <- pmatrix
                       pmatrixt[(ends[count-1]+1):pp,(ends[count-1]+1):pp] <- bbase$K
                       pmatrix <- pmatrixt
                    }
                    
                    B <- bbase$B
                    zw1 <- B[1:nrec,]
                    zw2 <- B[(nrec+1):(nrec+ngrid),]
                    zw3 <- B[(ngrid+1),]

                    #zw1 <- t(apply(zw1,1,normalize,b=zw3))
                    #zw2 <- t(apply(zw2,1,normalize,b=zw3))
                    
                    zmean <- cbind(zmean,zw3)
                    znew <- cbind(znew,zw2)
					z <- cbind(z,zw1)
                    nknots <- c(nknots,length(bbase$knots))
                    knots[[count]]<-bbase$knots
                 }
                 nsmooth <- nsmooth + dimen
             }
             out <- list(z=z,starts=starts,ends=ends,nknots=nknots,
                         knots=knots, nsmooth=nsmooth,
                         znew=znew, xreal=xreal, zmean=zmean,
                         pmatrix= pmatrix, pordv=pordv,
                         n.smoothers=n.smoothers,
                         nr.smoothers=nr.smoothers)
             return(out)
         }

       #########################################################################################
       # call parameters
       #########################################################################################
         m <- mcall <- cl <- match.call()
         nm <- names(m)[-1]
         keep <- is.element(nm, c("data", "na.action","offset"))
         for (i in nm[!keep]) m[[i]] <- NULL
         
         gp<-interpret.PSgam(formula) # interpret the formula 
         
         if (!is.null(n))
         {
            allvars<-c(NULL,"n")
            m$formula<-as.formula(paste(paste(deparse(gp$fake.formula,backtick=TRUE),collapse=""),
                           "+",paste(allvars, collapse = "+")))
         }                       
         else 
         {
            m$formula<-gp$fake.formula
         }   

         m$family<-NULL
         m$drop.unused.levels<-TRUE
         m[[1]]<-as.name("model.frame")
         pmf <- m
         mf <- eval.parent(m)
         if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
         
       #########################################################################################
       # checking input
       #########################################################################################
         if(is.null(n))
         {
              collapsed <- 0         
         }
         else
         {
              collapsed <- 1
              if(family$family!="binomial" || family$link!="logit")
              { 
                   stop("The total number of cases only can be used in the logit link .\n")     
              }
         }

       #########################################################################################
       # data structure
       #########################################################################################
         nrec <- dim(mf)[1]
         resp <- mf[,1]
         roffset <- model.offset(mf)
         if (is.null(roffset)) roffset <-rep(0,nrec)

         ntrials <- rep(1,nrec)
          
         if(gp$pfok==1)
         {
            x<-model.matrix(gp$pf ,data=mf)
            namesx <- colnames(x)
            if(is.null(n))
            {
               p<-dim(x)[2]
               if(is.null(p))
               {
                  p <- 1
                  namesx <- "(Intercept)"
               }
               nfixed <- p
            }
            else
            {
               p <- dim(x)[2]
               ntrials <- x[,p]
               x <- x[,-p]
               namesx <- namesx[-p]
               p <- p-1
               nfixed <- p
            }
         }
         else
         {
            nfixed <- 0
            p <- 1
            x <- matrix(0,nrow=nrec,ncol=1)
         }

         if(gp$psok==0)
         {
            nsmooth <- 0
            nsmooths <- 1
            q <- 1
            z <- matrix(0,nrow=nrec,ncol=1)
            starts <- rep(0,1)
            ends <- rep(0,1)
            
            znew <- matrix(0,nrow=ngrid,ncol=1)
            xreal <- matrix(0,nrow=ngrid,ncol=1)
            zmean <- matrix(0,nrow=1,ncol=1)
            pmatrix <- matrix(0,nrow=1,ncol=1)
            pordv <- rep(0,nsmooths) 
            possfp <- rep(0,nsmooths)
            
            nznh <- 1
            iapm <- rep(0,q+1)
            japm <- rep(0,nznh)
            apm <-  rep(0,nznh)
            
         }
         else
         { 
            tmpZ <- z.matrix(gp,mf,ngrid)
            z <- tmpZ$z
            q <- dim(z)[2]
            nsmooth <- tmpZ$nsmooth
            nsmooths <- nsmooth
            starts <- tmpZ$starts
            ends <- tmpZ$ends  
            znew <- tmpZ$znew
            zmean <- tmpZ$zmean
            xreal <- tmpZ$xreal
            pmatrix <- tmpZ$pmatrix
            pordv <- tmpZ$pordv
            n.smoothers <- tmpZ$n.smoothers
            nr.smoothers <- tmpZ$nr.smoothers
            
            if(nfixed>0)
            {
               possfp <- rep(0,nsmooths)
               for(i in 1:nsmooth)
               {
                   for(j in 1:length(namesx))
                   {
                       if(nr.smoothers[i]==namesx[j])possfp[i]<-j
                   }
               }
            }
            else
            {
               possfp <- rep(0,nsmooths)
            }
            rm(tmpZ)
            
            #############################
            # Sparse format for pmatrix
            #############################
            maxnh <- 2000
            indnh <- matrix(0,nrow=maxnh,ncol=3)
            nznh <- 0
            nr <- q
            
            # Hash table
            for(i in 1:nr)
            {
                for(j in i:nr)
                {
                    if(pmatrix[i,j]!=0)
                    {
                        y <- pmatrix[i,j]
                        foo <- .Fortran("hashm",
							y          =as.double(y),
							j1         =as.integer(i),
							k1         =as.integer(j),
							ind        =as.double(indnh),
							m          =as.integer(maxnh),
							nr         =as.integer(nznh),
							PACKAGE    ="DPpackage")	
			   
						indnh <- matrix(foo$ind,nrow=maxnh,ncol=3)
						nznh <- foo$nr
                        
                        if(i!=j)
                        {
                           foo <- .Fortran("hashm",
							y          =as.double(y),
							j1         =as.integer(j),
							k1         =as.integer(i),
							ind        =as.double(indnh),
							m          =as.integer(maxnh),
							nr         =as.integer(nznh),
							PACKAGE    ="DPpackage")	

							indnh <- matrix(foo$ind,nrow=maxnh,ncol=3)
							nznh <- foo$nr
						}   
                    }		
                }
            }
            
			# Linked list
              tmpv <- rep(0,nznh)
              iapm <- rep(0,q+1) 
              japm <- rep(0,nznh)
              apm <- rep(0,nznh)
              
              foo <- .Fortran("hashiajaa",
					x          =as.double(indnh),
					nhash      =as.integer(maxnh),
					n          =as.integer(nr),
					ia         =as.integer(iapm),
					ja         =as.integer(japm),
					a          =as.double(apm),
					m          =as.integer(nznh),
					tmp        =as.integer(tmpv),
					PACKAGE    ="DPpackage")	               

              iapm <- foo$ia
              japm <- foo$ja
              apm <- foo$a

        } 
        
        #pm2<-matrix(0,q,q)
        #for(i in 1:q)
        #{
        #    for(j in iapm[i]:(iapm[i+1]-1))
        #    {
        #        pm2[i,japm[j]] <- apm[j]
        #    }
        #}
        #print(pmatrix)
        #print(pm2)
        
        rm(indnh)
        rm(pmatrix)

       #########################################################################################
       # elements for Pseudo Countour Probabilities' computation
       #########################################################################################

         form <- gp$fake.formula
         Terms <- terms(form,data=mf)

         possiP <- NULL
         if(nfixed>0)
         {
            mat <- attr(Terms,"factors")
            namfact <- colnames(mat)
            nvar <- dim(mat)[1]
            nfact <- dim(mat)[2]
            possiP <- matrix(0,ncol=2,nrow=nfact)
            if (missing(data)) dataF <- model.frame(formula=form,xlev=NULL)
               dataF <- model.frame(formula=form,mf,xlev=NULL)
            namD <- names(dataF)
            isF <- sapply(dataF, function(x) is.factor(x) || is.logical(x))
            nlevel <- rep(0,nvar)
            for(i in 1:nvar)
            {
                if(isF[i])
                {
                   nlevel[i] <- length(table(dataF[[i]]))
                }
                else
                {
                   nlevel[i] <- 1
                }
            }
            startp <- 1+attr(Terms, "intercept")
            for(i in 1:nfact)
            {
                tmp1 <- 1
                for(j in 1:nvar)
                {
                    if(mat[j,i]==1 && isF[j])
                    {
                       tmp1 <- tmp1*(nlevel[j]-1)
                    }
                }
                endp <- startp+tmp1-1
                possiP[i,1] <- startp    
                possiP[i,2] <- endp
                startp <- endp+1
            }
            dimnames(possiP) <- list(namfact,c("Start","End"))
         }   

         
       #########################################################################################
       # prior information
       #########################################################################################

         if(nfixed==0)
         {
            prec <- matrix(0,nrow=1,ncol=1)
            bet0 <- rep(0,1)
            sb <- rep(0,1)
         }
         else
         {
            bet0 <- prior$beta0
            prec <- solve(prior$Sbeta0)
            sb <- prec%*%bet0

            if(length(bet0)!=p)
            { 
                   cat(p,"\n") 
                   stop("Error in the dimension of the mean of the normal prior for the fixed effects.\n")     
            }

            if(dim(prec)[1]!=p || dim(prec)[2]!=p)
            { 
                   stop("Error in the dimension of the covariance of the normal prior for the fixed effects.\n")     
            }

         }

         if(family$family=="Gamma" || family$family=="gaussian")
         {
            tau1 <- prior$tau1
            tau2 <- prior$tau2
            tau <- c(tau1,tau2)
            if(tau1<0 || tau2<0)
            { 
               stop("The parameters of the Gamma prior for the dispersion parameter must be possitive.\n")     
            }
         }   


         if(nsmooth>0)
         {
            taub1 <- prior$taub1
            taub2 <- prior$taub2
            taub <- c(taub1,taub2)
            if(taub1<0 || taub2<0)
            { 
               stop("The parameters of the Gamma prior for the penalty parameters must be possitive.\n")     
            }         
         }

       #########################################################################################
       # mcmc specification
       #########################################################################################
         if(missing(mcmc))
         {
            nburn <- 1000
            nsave <- 1000
            nskip <- 0
            ndisplay <- 100
            mcmcvec <- c(nburn,nskip,ndisplay)
            tune1 <- 1.1
         }
         else
         {
            mcmcvec <- c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
            nsave <- mcmc$nsave
            if(is.null(mcmc$tune1))
            {
               tune1 <- 1.1
            }
            {
               tune1 <- mcmc$tune1
            }
         }

       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2)
         dispp <- 0
         if(family$family=="Gamma")dispp <- 1
         if(family$family=="gaussian")dispp <- 1
         mc <- rep(0,5)
         randsave <- matrix(0,nrow=nsave,ncol=q)
         thetasave <- matrix(0,nrow=nsave,ncol=(p+dispp+nsmooth))
         pssave <- matrix(0,nrow=nsmooths*ngrid,ncol=2)

       #########################################################################################
       # parameters depending on status
       #########################################################################################
    	 if(status==TRUE)
		 {
			resp2 <- resp
			if(family$family=="binomial")
			{
				if(family$link=="logit")
				{
					if(!is.null(n)) resp2 <- cbind(resp,ntrials-resp)
				}
			}   

	        if(nfixed==0)
			{
				beta <- matrix(0,nrow=1,ncol=1)
				sigmab <- rep(0,nsmooths)
				b <- rep(0,q)
				if(nsmooth>0)
				{
					sigmab<-rep(1,nsmooths)
				}   
	        }

	        if(nfixed>0)
			{
				fit0 <- glm.fit(x, resp2, family= family,offset=roffset)   
				beta <- coefficients(fit0)
				b <- rnorm(q,0,0.1) 
				sigmab	<-	rep(0,nsmooths)
				if(nsmooth>0)
				{
					sigmab <- rep(1,nsmooths)
				}   
			}
			if(family$family=="Gamma") disp <- 1.1
			if(family$family=="gaussian") disp <- 1.1
		}	
		if(status==FALSE)
		{
			if(nfixed>0)
			{
	           beta <- state$beta
	        }
	        else
	        {
	           beta <- rep(0,p)
	        }
	        
			if(nsmooth>0)
			{
				b <- state$b
	            sigmab <- state$sigmab
	        }
	        else
	        {
	           b <- rep(0,q)
	           sigmab <- rep(0,q)
	        }
	        
	        if(family$family=="Gamma")disp <- 1/state$phi
	        if(family$family=="gaussian")disp <- 1/state$phi
		}

       #########################################################################################
       # calling the fortran code
       #########################################################################################

         if(family$family=="binomial")
         {

            if(family$link=="probit")
            {
                # specific working space

                  iflagp<-rep(0,p) 
                  betac<-rep(0,p)
                  workmhp1<-rep(0,p*(p+1)/2) 
                  workvp1<-rep(0,p) 
                  xtx<-matrix(0,nrow=p,ncol=p)
                  xty<-rep(0,p) 

                  maxq <- max(ends-starts+1)

                  iflagq <- rep(0,maxq) 
                  bc <- rep(0,maxq)
                  theta <- rep(0,maxq)
                  workmhq1 <- rep(0,maxq*(maxq+1)/2) 
                  workvq1 <- rep(0,maxq) 
                  ztz <- matrix(0,nrow=maxq,ncol=maxq)
                  zty <- rep(0,maxq) 
                  ztzinv <- matrix(0,nrow=maxq,ncol=maxq)

                  betasave <- rep(0,p)
                  bsave <- rep(0,q)
                  y <- rep(0,nrec)

                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)

                  workvps <- rep(0,nsmooths*ngrid)


                # fitting the model

                  foo <- .Fortran("psgamprob",
						nrec       =as.integer(nrec),
						nfixed     =as.integer(nfixed),
						p          =as.integer(p),
						nsmooth    =as.integer(nsmooth),
						q          =as.integer(q),
						starts     =as.integer(starts),
						ends       =as.integer(ends),
						nsmooths   =as.integer(nsmooths),
						maxq       =as.integer(maxq),
						x          =as.double(x),
						z          =as.double(z),
						yr         =as.integer(resp),
						roffset    =as.double(roffset),
                        ngrid      =as.integer(ngrid),
                        znew       =as.double(znew),
                        xreal      =as.double(xreal),
                        possfp     =as.integer(possfp),
						prec       =as.double(prec),	 
						sb         =as.double(sb),	  		
						taub       =as.double(taub),
						iapm       =as.integer(iapm),
						japm       =as.integer(japm),
						apm        =as.double(apm),
						nznh       =as.integer(nznh),
						pordv      =as.double(pordv),
                        beta       =as.double(beta),   
                        b          =as.double(b),   
                        sigmab     =as.double(sigmab),   
                        y          =as.double(y),   
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
						pssave     =as.double(pssave),
						seed       =as.integer(seed),
						iflagp     =as.integer(iflagp),
                        workmhp1   =as.double(workmhp1),
						workvp1    =as.double(workvp1),
						xtx        =as.double(xtx),
						xty        =as.double(xty),
						iflagq     =as.integer(iflagq),
						bc         =as.double(bc),
                        workmhq1   =as.double(workmhq1),
						workvq1    =as.double(workvq1),
						ztz        =as.double(ztz),
						zty        =as.double(zty),
						ztzinv     =as.double(ztzinv),
						theta      =as.double(theta),
						mc         =as.double(mc), 		
                        betasave   =as.double(betasave),
                        bsave      =as.double(bsave),
                        workvps    =as.double(workvps),
						PACKAGE    ="DPpackage")	
            }


            if(family$link=="logit")
            {
                # specific working space

                  acrate <- rep(0,1+nsmooth) 

                  iflagp <- rep(0,p) 
                  betac <- rep(0,p)
                  workmp1 <- matrix(0,nrow=p,ncol=p)
                  workmp2 <- matrix(0,nrow=p,ncol=p)
                  workmhp1 <- rep(0,p*(p+1)/2) 
                  workvp1 <- rep(0,p) 
                  xtx <- matrix(0,nrow=p,ncol=p)
                  xty <- rep(0,p) 

                  maxq <- max(ends-starts+1)

                  iflagq <- rep(0,maxq) 
                  bc <- rep(0,maxq)
                  theta <- rep(0,maxq)
                  workmq1 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmq2 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmhq1 <- rep(0,maxq*(maxq+1)/2) 
                  workvq1 <- rep(0,maxq) 
                  ztz <- matrix(0,nrow=maxq,ncol=maxq)
                  zty <- rep(0,maxq) 
                  ztzinv <- matrix(0,nrow=maxq,ncol=maxq)

                  betasave <- rep(0,p)
                  bsave <- rep(0,q)

                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)

                  resp <- cbind(resp,ntrials)
                  
                  workvps <- rep(0,nsmooths*ngrid)


                # fit the model

                  foo <- .Fortran("psgamlogit",
						nrec       =as.integer(nrec),
						nfixed     =as.integer(nfixed),
						p          =as.integer(p),
						nsmooth    =as.integer(nsmooth),
						q          =as.integer(q),
						starts     =as.integer(starts),
						ends       =as.integer(ends),
						nsmooths   =as.integer(nsmooths),
						maxq       =as.integer(maxq),
						x          =as.double(x),
						z          =as.double(z),
						y          =as.integer(resp),
						roffset    =as.double(roffset),
                        ngrid      =as.integer(ngrid),
                        znew       =as.double(znew),
                        xreal      =as.double(xreal),
                        possfp     =as.integer(possfp),
						bet0       =as.double(bet0),
						prec       =as.double(prec),	 
						sb         =as.double(sb),	  		
						taub       =as.double(taub),
						iapm       =as.integer(iapm),
						japm       =as.integer(japm),
						apm        =as.double(apm),
						nznh       =as.integer(nznh),
						pordv      =as.double(pordv),
                        beta       =as.double(beta),   
                        b          =as.double(b),   
                        sigmab     =as.double(sigmab),   
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
                        acrate     =as.double(acrate),   
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
						pssave     =as.double(pssave),
						seed       =as.integer(seed),
						iflagp     =as.integer(iflagp),
						betac      =as.double(betac),
						workmp1    =as.double(workmp1),
						workmp2    =as.double(workmp2),
                        workmhp1   =as.double(workmhp1),
						workvp1    =as.double(workvp1),
						xtx        =as.double(xtx),
						xty        =as.double(xty),
						iflagq     =as.integer(iflagq),
						bc         =as.double(bc),
						workmq1    =as.double(workmq1),
						workmq2    =as.double(workmq2),
                        workmhq1   =as.double(workmhq1),
						workvq1    =as.double(workvq1),
						ztz        =as.double(ztz),
						zty        =as.double(zty),
						ztzinv     =as.double(ztzinv),
						theta      =as.double(theta),
						mc         =as.double(mc), 		
                        betasave   =as.double(betasave),
                        bsave      =as.double(bsave),
                        workvps    =as.double(workvps),
						PACKAGE    ="DPpackage")	

            }
         }

         if(family$family=="poisson")
         {
            if(family$link=="log")
            {
                # specific working space

                  acrate <- rep(0,1+nsmooth) 

                  iflagp <- rep(0,p) 
                  betac <- rep(0,p)
                  workmp1 <- matrix(0,nrow=p,ncol=p)
                  workmp2 <- matrix(0,nrow=p,ncol=p)
                  workmhp1 <- rep(0,p*(p+1)/2) 
                  workvp1 <- rep(0,p) 
                  xtx <- matrix(0,nrow=p,ncol=p)
                  xty <- rep(0,p) 

                  maxq <- max(ends-starts+1)

                  iflagq <- rep(0,maxq) 
                  bc <- rep(0,maxq)
                  theta <- rep(0,maxq)
                  workmq1 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmq2 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmhq1 <- rep(0,maxq*(maxq+1)/2) 
                  workvq1 <- rep(0,maxq) 
                  ztz <- matrix(0,nrow=maxq,ncol=maxq)
                  zty <- rep(0,maxq) 
                  ztzinv <- matrix(0,nrow=maxq,ncol=maxq)

                  betasave <- rep(0,p)
                  bsave <- rep(0,q)

                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)

                  workvps <- rep(0,nsmooths*ngrid)


                # fit the model

                  foo <- .Fortran("psgampois",
						nrec       =as.integer(nrec),
						nfixed     =as.integer(nfixed),
						p          =as.integer(p),
						nsmooth    =as.integer(nsmooth),
						q          =as.integer(q),
						starts     =as.integer(starts),
						ends       =as.integer(ends),
						nsmooths   =as.integer(nsmooths),
						maxq       =as.integer(maxq),
						x          =as.double(x),
						z          =as.double(z),
						y          =as.integer(resp),
						roffset    =as.double(roffset),
                        ngrid      =as.integer(ngrid),
                        znew       =as.double(znew),
                        xreal      =as.double(xreal),
                        possfp     =as.integer(possfp),
						bet0       =as.double(bet0),
						prec       =as.double(prec),	 
						sb         =as.double(sb),	  		
						taub       =as.double(taub),
						iapm       =as.integer(iapm),
						japm       =as.integer(japm),
						apm        =as.double(apm),
						nznh       =as.integer(nznh),
						pordv      =as.double(pordv),
                        beta       =as.double(beta),   
                        b          =as.double(b),   
                        sigmab     =as.double(sigmab),   
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
                        acrate     =as.double(acrate),   
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
						pssave     =as.double(pssave),
						seed       =as.integer(seed),
						iflagp     =as.integer(iflagp),
						betac      =as.double(betac),
						workmp1    =as.double(workmp1),
						workmp2    =as.double(workmp2),
                        workmhp1   =as.double(workmhp1),
						workvp1    =as.double(workvp1),
						xtx        =as.double(xtx),
						xty        =as.double(xty),
						iflagq     =as.integer(iflagq),
						bc         =as.double(bc),
						workmq1    =as.double(workmq1),
						workmq2    =as.double(workmq2),
                        workmhq1   =as.double(workmhq1),
						workvq1    =as.double(workvq1),
						ztz        =as.double(ztz),
						zty        =as.double(zty),
						ztzinv     =as.double(ztzinv),
						theta      =as.double(theta),
						mc         =as.double(mc), 		
                        betasave   =as.double(betasave),
                        bsave      =as.double(bsave),
                        workvps    =as.double(workvps),
						PACKAGE    ="DPpackage")	
            }
         }  

         if(family$family=="Gamma")
         {
            if(family$link=="log")
            {
                # specific working space

                  acrate <- rep(0,2+nsmooth) 

                  iflagp <- rep(0,p) 
                  betac <- rep(0,p)
                  workmp1 <- matrix(0,nrow=p,ncol=p)
                  workmp2 <- matrix(0,nrow=p,ncol=p)
                  workmhp1 <- rep(0,p*(p+1)/2) 
                  workvp1 <- rep(0,p) 
                  xtx <- matrix(0,nrow=p,ncol=p)
                  xty <- rep(0,p) 

                  maxq <- max(ends-starts+1)

                  iflagq <- rep(0,maxq) 
                  bc <- rep(0,maxq)
                  theta <- rep(0,maxq)
                  workmq1 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmq2 <- matrix(0,nrow=maxq,ncol=maxq)
                  workmhq1 <- rep(0,maxq*(maxq+1)/2) 
                  workvq1 <- rep(0,maxq) 
                  ztz <- matrix(0,nrow=maxq,ncol=maxq)
                  zty <- rep(0,maxq) 
                  ztzinv <- matrix(0,nrow=maxq,ncol=maxq)

                  betasave <- rep(0,p+1)
                  bsave <- rep(0,q)

                  seed1 <- sample(1:29000,1)
                  seed2 <- sample(1:29000,1)
                  seed <- c(seed1,seed2)

                  workvps <- rep(0,nsmooths*ngrid)


                # fit the model

                  foo <- .Fortran("psgamgam",
						nrec       =as.integer(nrec),
						nfixed     =as.integer(nfixed),
						p          =as.integer(p),
						nsmooth    =as.integer(nsmooth),
						q          =as.integer(q),
						starts     =as.integer(starts),
						ends       =as.integer(ends),
						nsmooths   =as.integer(nsmooths),
						maxq       =as.integer(maxq),
						x          =as.double(x),
						z          =as.double(z),
						y          =as.double(resp),
						roffset    =as.double(roffset),
                        ngrid      =as.integer(ngrid),
                        znew       =as.double(znew),
                        xreal      =as.double(xreal),
                        possfp     =as.integer(possfp),
						bet0       =as.double(bet0),
						prec       =as.double(prec),	 
						sb         =as.double(sb),	  		
						taub       =as.double(taub),
						iapm       =as.integer(iapm),
						japm       =as.integer(japm),
						apm        =as.double(apm),
						nznh       =as.integer(nznh),
						pordv      =as.double(pordv),
						tau        =as.double(tau),
                        beta       =as.double(beta),   
                        b          =as.double(b),   
                        sigmab     =as.double(sigmab),   
                        disp       =as.double(disp),   
						mcmc       =as.integer(mcmcvec),
						nsave      =as.integer(nsave),
						tune1      =as.double(tune1),
                        acrate     =as.double(acrate),   
						cpo        =as.double(cpo),
						randsave   =as.double(randsave),
						thetasave  =as.double(thetasave),
						pssave     =as.double(pssave),
						seed       =as.integer(seed),
						iflagp     =as.integer(iflagp),
						betac      =as.double(betac),
						workmp1    =as.double(workmp1),
						workmp2    =as.double(workmp2),
                        workmhp1   =as.double(workmhp1),
						workvp1    =as.double(workvp1),
						xtx        =as.double(xtx),
						xty        =as.double(xty),
						iflagq     =as.integer(iflagq),
						bc         =as.double(bc),
						workmq1    =as.double(workmq1),
						workmq2    =as.double(workmq2),
                        workmhq1   =as.double(workmhq1),
						workvq1    =as.double(workvq1),
						ztz        =as.double(ztz),
						zty        =as.double(zty),
						ztzinv     =as.double(ztzinv),
						theta      =as.double(theta),
						mc         =as.double(mc), 		
                        betasave   =as.double(betasave),
                        bsave      =as.double(bsave),
                        workvps    =as.double(workvps),
						PACKAGE    ="DPpackage")	
            }
         }   


         if(family$family=="gaussian")
         {
            # specific working space

              iflagp <- rep(0,p) 
              betac <- rep(0,p)
              workmhp1 <- rep(0,p*(p+1)/2) 
              workvp1 <- rep(0,p) 
              xtx <- matrix(0,nrow=p,ncol=p)
              xty <- rep(0,p) 

              maxq <- max(ends-starts+1)

              iflagq <- rep(0,maxq) 
              bc <- rep(0,maxq)
              theta <- rep(0,maxq)
              workmhq1 <- rep(0,maxq*(maxq+1)/2) 
              workvq1 <- rep(0,maxq) 
              ztz <- matrix(0,nrow=maxq,ncol=maxq)
              zty <- rep(0,maxq) 
              ztzinv <- matrix(0,nrow=maxq,ncol=maxq)

              betasave <- rep(0,p+1)
              bsave <- rep(0,q)

              seed1 <- sample(1:29000,1)
              seed2 <- sample(1:29000,1)
              seed <- c(seed1,seed2)

              workvps <- rep(0,nsmooths*ngrid)


            # fit the model

              foo <- .Fortran("psgamgauss",
					nrec       =as.integer(nrec),
					nfixed     =as.integer(nfixed),
					p          =as.integer(p),
					nsmooth    =as.integer(nsmooth),
					q          =as.integer(q),
					starts     =as.integer(starts),
					ends       =as.integer(ends),
					nsmooths   =as.integer(nsmooths),
					maxq       =as.integer(maxq),
					x          =as.double(x),
					z          =as.double(z),
					y          =as.double(resp),
					roffset    =as.double(roffset),
					ngrid      =as.integer(ngrid),
					znew       =as.double(znew),
					xreal      =as.double(xreal),
					possfp     =as.integer(possfp),
					prec       =as.double(prec),	 
					sb         =as.double(sb),	  		
					taub       =as.double(taub),
					iapm       =as.integer(iapm),
					japm       =as.integer(japm),
					apm        =as.double(apm),
					nznh       =as.integer(nznh),
					pordv      =as.double(pordv),
					tau        =as.double(tau),
					beta       =as.double(beta),   
					b          =as.double(b),   
					sigmab     =as.double(sigmab),   
					disp       =as.double(disp),   
					mcmc       =as.integer(mcmcvec),
					nsave      =as.integer(nsave),
					cpo        =as.double(cpo),
					randsave   =as.double(randsave),
					thetasave  =as.double(thetasave),
					pssave     =as.double(pssave),
					seed       =as.integer(seed),
					iflagp     =as.integer(iflagp),
					workmhp1   =as.double(workmhp1),
					workvp1    =as.double(workvp1),
					xtx        =as.double(xtx),
					xty        =as.double(xty),
					iflagq     =as.integer(iflagq),
					bc         =as.double(bc),
					workmhq1   =as.double(workmhq1),
					workvq1    =as.double(workvq1),
					ztz        =as.double(ztz),
					zty        =as.double(zty),
					ztzinv     =as.double(ztzinv),
					theta      =as.double(theta),
					mc         =as.double(mc), 		
					betasave   =as.double(betasave),
					bsave      =as.double(bsave),
					workvps    =as.double(workvps),
					PACKAGE    ="DPpackage")	
         }   


       #########################################################################################
       # save state
       #########################################################################################
        
         if(family$family=="Gamma" || family$family=="gaussian")
         {
            phi <- 1/foo$disp
         }
         else
         {
            phi <- 0
         }
         mc <- foo$mc
         names(mc) <- c("Dbar", "Dhat", "pD", "DIC","LPML")
         
         dimen <- p+dispp+nsmooth
         thetasave <- matrix(foo$thetasave,nrow=nsave, ncol=dimen)
         randsave <- matrix(foo$randsave,nrow=nsave, ncol=q)
         pssave <- matrix(foo$pssave,nrow=nsmooths*ngrid, ncol=2)

         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)         
         cpo <- cpom[,1]         
         fso <- cpom[,2]

         if(nfixed==0)
         {
            pnames1 <- NULL
         }
         if(nfixed>0)
         {
            pnames1 <- namesx
         }
         
         if(nsmooth==0)
         {
            pnames2 <- NULL
         }
         if(nsmooth>0)
         {
            pnames2 <- n.smoothers
         }
         
         if(family$family=="Gamma" || family$family=="gaussian")
         {
         	colnames(thetasave) <- c(pnames1,"phi",pnames2)
         }
         else
         {
                colnames(thetasave)<-c(pnames1,pnames2)
         }

		 model.name <- "Bayesian semiparametric generalized additive model using P-Splines"
	 
         coeff <- apply(thetasave,2,mean)		

		 state <- list(	b=foo$b,
						beta=foo$beta,
						sigmab=foo$sigmab,
						phi=phi)

		 save.state <- list(thetasave=thetasave,randsave=randsave,
							pssave=pssave)

         acrate <- foo$acrate

		 out <- list(modelname=model.name,
		 			 coefficients=coeff,
					 call=cl,
					 prior=prior,
					 mcmc=mcmc,
					 state=state,
					 save.state=save.state,
					 nrec=foo$nrec,
					 nfixed=foo$nfixed,
					 nsmooth=foo$nsmooth,
					 q=foo$q,
					 dispp=dispp,
					 cpo=cpo,
					 fso=fso,
					 prior=prior,
					 x=x,
					 z=z,
					 mf=mf,
					 dimen=dimen,
					 acrate=acrate,
					 possiP=possiP,
					 mc=mc,
					 starts=starts,
					 ends=ends,
					 xreal=xreal,
					 n.smoothers=n.smoothers,
					 nr.smoothers=nr.smoothers,
					 possfp=possfp,
					 formula=formula)
                 
         cat("\n\n")        

         class(out) <- c("PSgam")
         return(out) 

}


###                    
### Tools: anova, print, summary, plot
###
### Copyright: Alejandro Jara, 2007
### Last modification: 23-07-2007.


"anova.PSgam"<-function(object, ...)
{

######################################################################################
cregion<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a simultaneous credible region for a vector 
#  parameter from the MCMC sample
# 
#  Reference: Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995)
#             Bayesian computation and stochastic systems (with Discussion)
#             Statistical Science, vol. 10, 3 - 66, page 30
#  and        Held, L. (2004) Simultaneous inference in risk assessment; a Bayesian 
#             perspective In: COMPSTAT 2004, Proceedings in Computational 
#             Statistics (J. Antoch, Ed.) 213 - 222, page 214
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are computed.
######################################################################################
{
    #Basic information
     nmonte<-dim(x)[1]
     p<-dim(x)[2]
     
    #Ranks for each component
     ranks <- apply(x, 2, rank, ties.method="first")
     
    #Compute the set S={max(nmonte+1-min r_i(t) , max r_i(t)): t=1,..,nmonte}
     left <- nmonte + 1 - apply(ranks, 1, min)
     right <- apply(ranks, 1, max)
     S <- apply(cbind(left, right), 1, max)
     S <- S[order(S)]
    
    #Compute the credible region
     k <- floor(nmonte*probs)     
     tstar <- S[k]
     out<-list()
     for(i in 1:length(tstar))
     {
        upelim <- x[ranks == tstar[i]]
        lowlim <- x[ranks == nmonte + 1 - tstar[i]]    
        out[[i]] <- rbind(lowlim, upelim)
        rownames(out[[i]]) <- c("Lower", "Upper")
        colnames(out[[i]]) <- colnames(x)
     }
     names(out) <- paste(probs)
     return(out)
}

######################################################################################
cint<-function(x,probs=c(0.90,0.975))
######################################################################################
#  Function to compute a credible interval from the MCMC sample
#
#  Arguments 
#  sample : a data frame or matrix with sampled values (one column = one parameter).
#  probs  : probabilities for which the credible regions are to be computed.
######################################################################################
{
    #Compute the credible interval
     delta<-(1-probs)/2
     lprobs<-cbind(delta,probs+delta) 
     out<-matrix(quantile(x,probs=lprobs),ncol=2)
     colnames(out) <- c("Lower","Upper")
     rownames(out) <- paste(probs)
     return(out)
}

######################################################################################
hnulleval<-function(mat,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     npar<-dim(mat)[2]   
     lower<-rep(0,npar)
     upper<-rep(0,npar)
     for(i in 1:npar)
     {
        lower[i]<-mat[1,i]< hnull[i]
        upper[i]<-mat[2,i]> hnull[i]
     }
     total<-lower+upper
     out<-(sum(total==2) == npar)
     return(out)
}

######################################################################################
hnulleval2<-function(vec,hnull)
######################################################################################
#  Evaluate H0
#  AJV, 2006
######################################################################################
{
     lower<-vec[1]< hnull
     upper<-vec[2]> hnull

     total<-lower+upper
     out<-(total==2)
     return(out)
}


######################################################################################
pcp<-function(x,hnull=NULL,precision=0.001,prob=0.95,digits=digits)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Region)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-rep(0,dim(x)[2])
    if (dim(x)[2]!=length(hnull)) stop("Dimension of x and hnull must be equal!!")

    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cregion(x,probs=probsf)

    is.hnull <- hnulleval(cr[[2]],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval(cr[[length(cr)]],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval(cr[[i]],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1 - probs[first]
       }
    }
    output <- list(cr=cr[[1]], prob=prob, pval=pval,hnull=hnull)
    return(output)
}


######################################################################################
pcp2<-function(x,hnull=NULL,precision=0.001,prob=0.95)
######################################################################################
#  Function to compute Pseudo Countour Probabilities (Interval)
#  AJV, 2006
######################################################################################
{
    if(is.null(hnull))hnull<-0
    probs <- seq(precision, 1-precision, by=precision)
    neval <- length(probs)
    probsf <- c(prob,probs)
    cr <-  cint(x,probs=probsf)

    is.hnull <- hnulleval2(cr[2,],hnull)
    if(is.hnull)
    {
       pval <- 1-precision
    }   
    else
    {
       is.hnull <- hnulleval2(cr[(neval+1),],hnull)
       if (!is.hnull) 
       {
         pval <- precision
       }  
       else
       {
         is.hnull<-rep(0,neval+1)
         for(i in 1:(neval+1))
         {
            is.hnull[i] <- hnulleval2(cr[i,],hnull)
         }   
         is.hnull <- is.hnull[-1]
         first <- neval - sum(is.hnull) + 1
         pval <- 1-probs[first]
       }
    }
    output <- list(cr=cr[1,], prob=prob, pval=pval,hnull=hnull)
    return(output)
}

######################################################################################
######################################################################################
######################################################################################

	 P <- NULL
	 df <- NULL
	 tmp <- NULL


	 possiP <- object$possiP
	 nfact <- nrow(possiP)
        
	 nsmooth <- object$nsmooth

     if(object$nfixed>1)
     { 

	 for(i in 1:nfact)
	 {
		   if(sum(object$possiP[i,1]==object$possfp)==0)
		   {
				tmp <- c(tmp,rownames(possiP)[i])
				if((possiP[i,2]-possiP[i,1])>0)
				{ 
					df <- c(df,(possiP[i,2]-possiP[i,1])+1)
					x <- matrix(object$save.state$thetasave[,possiP[i,1]:possiP[i,2]])
					foo <- pcp(x=x) 
					P <- c(P,foo$pval)
				}else
				{
					df <- c(df,1)
					x <- object$save.state$thetasave[,possiP[i,1]:possiP[i,2]]
					foo <- pcp2(x=x) 
					P <- c(P,foo$pval)
				}
			}
	   }

	 for(i in 1:nfact)
	 {
		   if(sum(object$possiP[i,1]==object$possfp)>0)
		   {
				tmp <- c(tmp,rownames(possiP)[i])
				if((possiP[i,2]-possiP[i,1])>0)
				{ 
					df <- c(df,(possiP[i,2]-possiP[i,1])+1)
					x <- matrix(object$save.state$thetasave[,possiP[i,1]:possiP[i,2]])
					foo <- pcp(x=x) 
					P <- c(P,foo$pval)
				}else
				{
					df <- c(df,1)
					x <- object$save.state$thetasave[,possiP[i,1]:possiP[i,2]]
					foo <- pcp2(x=x) 
					P <- c(P,foo$pval)
				}
			}
	   }
       }

       possiQ <- cbind(object$starts,object$ends)
       for(i in 1:nsmooth)
       {
           if((possiQ[i,2]-possiQ[i,1])>0)
           { 
			  df <- c(df,(possiQ[i,2]-possiQ[i,1])+1)
              x <- matrix(object$save.state$randsave[,possiQ[i,1]:possiQ[i,2]])
              foo <- pcp(x=x) 
              P <- c(P,foo$pval)
           }
           else
           {
			  df <- c(df,1)
              x <- object$save.state$randsave[,possiQ[i,1]:possiQ[i,2]]
              foo <- pcp2(x=x) 
              P <- c(P,foo$pval)
           }
       }

       table <- data.frame(df,P) 

	   tmp <- c(tmp,object$n.smoothers)

       dimnames(table) <- list(tmp, c("Df","PsCP"))
       structure(table, heading = c("Table of Pseudo Contour Probabilities\n", 
        paste("Response:", deparse(formula(object$formula)[[2]]))), class = c("anovaPsCP",
        "data.frame"))
}



"print.PSgam"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) {
        cat("Posterior Inference of Parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if(is.null(x$acrate)) cat(" ")
    else
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   
   
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"summary.PSgam"<-function(object, hpd=TRUE, ...) 
{
    stde<-function(x)
    {
    	n<-length(x)
    	return(sd(x)/sqrt(n))
    }

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

    thetasave<-object$save.state$thetasave

    dimen1<-object$dimen-object$nsmooth
    dimen2<-dim(thetasave)[2]-dimen1

### Parametric part of the model

    if(dimen1==1)
    {
       mat<-matrix(thetasave[,1:dimen1],ncol=1) 
    }
    else
    {
       mat<-thetasave[,1:dimen1]
    }

    coef.p<-object$coefficients[1:dimen1]
    coef.m <-apply(mat, 2, median)    
    coef.sd<-apply(mat, 2, sd)
    coef.se<-apply(mat, 2, stde)

    if(hpd)
    {             
         limm <- apply(mat, 2, hpdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }
    else
    {
         limm<-apply(mat, 2, pdf)
         coef.l <- limm[1,]
         coef.u <- limm[2,]
    }

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    if(hpd)
    {
         dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
    }
    else
    {
         dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
    }

    ans <- c(object[c("call", "modelname")])
    ans$coefficients<-coef.table


### Penalties

    if(dimen2==1)
    {
       mat<-matrix(thetasave[,(dimen1+1):(dimen1+dimen2)],ncol=1) 
    }
    else
    {
       mat<-thetasave[,(dimen1+1):(dimen1+dimen2)]
    }

    coef.p<-object$coefficients[(dimen1+1):(dimen1+dimen2)]
    coef.m <-apply(mat, 2, median)    
    coef.sd<-apply(mat, 2, sd)
    coef.se<-apply(mat, 2, stde)

    if(hpd)
    {             
         limm<-apply(mat, 2, hpdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }
    else
    {
         limm<-apply(mat, 2, pdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    if(hpd)
    {
         dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
    }
    else
    {
         dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
    }

    ans$pen<-coef.table

### CPO
    ans$cpo<-object$cpo

### Model comparison
   
    coef.table<-matrix(object$mc,nrow=1,ncol=5)
    dimnames(coef.table) <- list(" ", c("Dbar", "Dhat", "pD", "DIC","LPML"))
    ans$mc<-coef.table
    
    ans$nrec<-object$nrec
    ans$acrate<-object$acrate

    class(ans) <- "summaryPSgam"
    return(ans)
}


"print.summaryPSgam"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 

    cat("\nModel's performance:\n")
    print.default(format(x$mc, digits = digits), print.gap = 2, 
    quote = FALSE)

    if (length(x$coefficients)) {
        cat("\nParametric component:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$pen)) {
        cat("\nPenalty parameters:\n")
        print.default(format(x$pen, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No smoothers\n")

    if(is.null(x$acrate)) cat("")
    else
    {
       cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
    }   
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"plot.PSgam"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
{

	fancydensplot1<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
	# Author: AJV, 2006
	#
	{
		dens <- density(x,n=npts)
		densx <- dens$x
		densy <- dens$y

		meanvar <- mean(x)
		densx1 <- max(densx[densx<=meanvar])
		densx2 <- min(densx[densx>=meanvar])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ymean <- densy1 + ((densy2-densy1)/(densx2-densx1))*(meanvar-densx1)
        

        if(hpd==TRUE)
		{
			alpha<-0.05
			alow<-rep(0,2)
        	aupp<-rep(0,2)
        	n<-length(x)
			a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(x),
		                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
			xlinf<-a$alow[1]            
			xlsup<-a$aupp[1]            
		}
		else
		{
			xlinf <- quantile(x,0.025)
			xlsup <- quantile(x,0.975)
		}

		densx1 <- max(densx[densx<=xlinf])
		densx2 <- min(densx[densx>=xlinf])
		densy1 <- densy[densx==densx1]
		densy2 <- densy[densx==densx2]
		ylinf <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlinf-densx1)

		densx1 <- max(densx[densx<=xlsup])
		densx2 <- min(densx[densx>=xlsup])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ylsup <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlsup-densx1)

        plot(0.,0.,xlim = c(min(densx), max(densx)), ylim = c(min(densy), max(densy)),
             axes = F,type = "n" , xlab=xlab, ylab=ylab, main=main, cex=1.2)

        
        xpol<-c(xlinf,xlinf,densx[densx>=xlinf & densx <=xlsup],xlsup,xlsup)
        ypol<-c(0,ylinf,densy[densx>=xlinf & densx <=xlsup] ,ylsup,0)
             
        polygon(xpol, ypol, border = FALSE,col=col)
        
        lines(c(min(densx), max(densx)),c(0,0),lwd=1.2)
        
        segments(min(densx),0, min(densx),max(densy),lwd=1.2)
        
        lines(densx,densy,lwd=1.2)
             
        segments(meanvar, 0, meanvar, ymean,lwd=1.2)
        segments(xlinf, 0, xlinf, ylinf,lwd=1.2)
        segments(xlsup, 0, xlsup, ylsup,lwd=1.2)

		axis(1., at = round(c(xlinf, meanvar,xlsup), 2.), labels = T,pos = 0.)
        axis(1., at = round(seq(min(densx),max(densx),length=15), 2.), labels = F,pos = 0.)
        axis(2., at = round(seq(0,max(densy),length=5), 2.), labels = T,pos =min(densx))
	}

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }


    if(is(x, "PSgam"))
    {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))

           for(i in 1:n)
           {
               title1 <- paste("Trace of",pnames[i],sep=" ")
               title2 <- paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           ngrid <- nrow(x$xreal)
           nsmooth <- ncol(x$xreal)
           meanss <- x$save.state$pssave[,1]
           sdss <- sqrt(x$save.state$pssave[,2])

		   #limm <- apply(x$save.state$pssave[,1], 2, hpdf)
		   #linf <- limm[1,]
		   #lsup <- limm[2,]

           linf <- meanss-sdss
           lsup <- meanss+sdss

           layout(matrix(seq(1,1,1), nrow=1 , ncol=1 ,byrow=TRUE))           
           for (i in 1:nsmooth)
           {
               begs <- (i-1)*ngrid+1
               ends <- i*ngrid
               
               ymin <- min(linf[begs:ends],lsup[begs:ends])
               ymax <- max(linf[begs:ends],lsup[begs:ends])
               
               dff <- diff(c(ymin,ymax))
               ymin <- ymin-0.07*dff
               ymax <- ymax+0.07*dff
               
               xreal <- x$mf[colnames(x$mf)==x$nr.smoothers[i]][[1]]
 
               plot(x$xreal[,i],meanss[begs:ends],type="l",
                    xlab=x$nr.smoothers[i],ylab=x$n.smoothers[i],lwd=1,
                    ylim=c(ymin,ymax))    
               lines(x$xreal[,i],linf[begs:ends],lwd=1,lty=2)    
               lines(x$xreal[,i],lsup[begs:ends],lwd=1,lty=2)    
               rug(xreal)

           }
          
        }
   
        else
        {
            coef.p <- x$coefficients
			n <- length(coef.p)
			pnames <- names(coef.p)
			poss <- 0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }
            if (poss==0) 
			{
				stop("This parameter is not present in the original model.\n")
			}
	    
			par(ask = ask)
			layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
            title1 <- paste("Trace of",pnames[poss],sep=" ")
            title2 <- paste("Density of",pnames[poss],sep=" ")       
            plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
            fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)

        }
   }

}

