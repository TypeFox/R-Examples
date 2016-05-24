      ##############################################################################
      #### subsection 1: Function to get full tables for  association studies #########
      ###############################################################################
      fulltable <- function(twocoldata)
      {
        # m is the # rows in the data
        m <-nrow(twocoldata)
        # construct 3-3 full table
        tables<- array(0, dim=c(3,3,m))

        controlsum<-sum(twocoldata[,1])
        treatsum<-sum(twocoldata[,2])

         # get 2-by-2 table for each gene
          tables[1,1,]=twocoldata[,1]
          tables[1,2,]=twocoldata[,2]
          tables[2,1,]=controlsum - twocoldata[,1]
          tables[2,2,]= treatsum -twocoldata[,2]

         # get 3rd row and column for each gene
         tables[1,3,]= tables[1,1,]+ tables[1,2,]
         tables[2,3,]= tables[2,1,]+ tables[2,2,]

           tables[3,1,]=rep(controlsum,m)
           tables[3,2,]=rep(treatsum,m)
           tables[3,3,]=rep(treatsum+controlsum,m)

         return(tables)
      }


      ################################################################# ###########################
      #### Subsection 3: Function to Get marginals and cellcounts from generated counts
      #####             for association stuies
      ################################################################# ###########################

      getcellcountsandmarginals <- function(countveccontrol,countvectreat)
      {
        m <- length(countveccontrol)
        twocolsimdata <- cbind(countveccontrol,countvectreat)
        simtables <- fulltable(twocolsimdata)

        simallcellcounts <- vector("list",m)
        for (i in 1:m) {simallcellcounts[[i]] <- simtables[1:2,1:2,i]}

        simallmarginals <- vector("list",m)
        for (i in 1:m) {simallmarginals[[i]] <- c(simtables[1,3,i],simtables[2,3,i],simtables[3,1,i])}

         y2 <- list(simallcellcounts,simallmarginals)
          return(y2)
      }

      ###################################################################################
      #### Subsection 4: Function to Get supports under true null
      #####                when test statistic is non-central hypergeometric
      ###################################################################################

      ## Function to get support of pvalues,  marginals is a vector as a list component
      # pvalue support is always computed under true null, so set pis=1
      pvalueSupport <- function(marginal)
      {   # if dim(marginals)!=c(3,3),report error
          masstmp <- dnoncenhypergeom(NA,marginal[1],marginal[2],marginal[3],1)
          mass <- masstmp[,2]
          lgh2 = length(mass)
          # change 11/21: matrix into double
           temp <- double(lgh2)
          # temp <- matrix(0,nrow=lgh2,ncol=1)
          for (i in 1:lgh2)
          {
            # temp[i] is the pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*mass)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)

          # save mean to the fist entry
          support<- c(meantmp,psupport)
          return(support)
          }

      ################################################################# ######### ######### ######### #########
      #### Subsection 5: Function to compuete deltas (lambda - F(lambda) for adjusted estimator         #########
      ################################################################# ######### #########  ######### #########
      # works for each support
      deviations <- function(lambda,asupport)
      {
            # remove the lst entry of asupport since it is the mean
            bigornot <- lambda >= asupport[-1]
            bignum <- sum(bigornot)

             ### if lambda smaller than the smallest pvalue, then delta[i]=lambda
            if ( bignum ==0)    { delta=lambda}

            # if lambda falls into an interval formed by successive points in support
            if ( bignum > 0)
            {      fall= max(which(bigornot))
             delta=lambda- asupport[-1][fall] }
          return(delta)
      }

    ################################################
    # following devation function specifically designed
    # the only different is that the first two entries of asupport
    # are the mean under null and pvalue under null, then comes the
    # support
    deviationsExtPvalSupp <- function(lambda,asupport)
    {
          # remove the lst entry of asupport since it is the mean
          # remove the 2nd entry of asuppot since it is the pvalue
          bigornot <- lambda >= asupport[-1][-1]
          bignum <- sum(bigornot);
    
           ### if lambda smaller than the smallest pvalue, then delta[i]==lambda
          if ( bignum ==0)    { delta=lambda}
    
          # if lambda falls into an interval formed by successive points in support
          if ( bignum > 0)
          {      fall= max(which(bigornot))
           delta=lambda- asupport[-1][-1][fall] }
        return(delta)
    }

      ###################################################################################
      #### Subsection 6: Function to Compute Storey's and Adjusted estimators    #########
      ###################################################################################

      twoestimators <- function(lambda,epsilon,pvector,deltavec)
      {
          # m is the coherent length of the argument
           m <- length(pvector)

          over <- m*(1-lambda)
          stunmin <- sum(pvector > lambda)/over 
          st <- min(1,stunmin)
          if (is.nan(st))  print("Storey's estimator of pi0 is NaN")
          if (is.na(st))  print("Storey's estimator of pi0 is NA")
          adjst <- max(0,min(1,stunmin - sum(epsilon*deltavec)/over)) 
                                                                       
          if (is.nan(adjst))  print("Generalized estimator of pi0 is NaN")
          if (is.na(adjst))  print("Generalized estimator of pi0 is NA")
          
          propets <- c(st,adjst)
          return(propets)
      }


      #################################################################################
      ####   Subsection 10: Distribution of the p-values  (Generally Applicable)  #########
      #################################################################################

       pvalueDist <- function(peva,support)
          {
            lgh=length(support)
            lgh3 <- length(peva)
            y <- double(lgh3)
            for (i in 1:lgh3)  {
              cpv <- peva[i] >= support
              # the key is to compare where t falls in the support
              if (sum(cpv)==0) {y[i]=0} # because t is strictly less than minimum of support
              else if (sum(cpv)==lgh) {y[i]=1} # because t is no less than maximum of support
              else # t falls in one interval
              {  y[i] <- support[max(which(cpv))]   }
            }
            return(y)
            }


      ################################################################################### #######################
      #### Subsection 14: pvalue supports, under true null  Test Stat ~ Bino for equality of two Poisson   #########
      ###################################################################################   #####################

      pvalueByBinoSupport <- function(marginal)
      {
          #  under null that two poisson rates are equal, the UMP test follows Binomial with p=0.5
          # get all possible masses
          mass <- dbinom(seq(from=0,to=marginal[2]),marginal[2],0.5)

          ## compute pvalue
          realizedmass <- dbinom(marginal[1],marginal[2],0.5)

           # temp[i] is the two-sided pvalue for i-th table given marginals
           pvalue <- sum(mass[which(mass <= realizedmass)])
           if (is.na(pvalue)) print("pvalue is NaN")

          ### compute p value support
          lgh2 = length(mass)
          temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(mass[which(mass <= mass[i])])
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*mass)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)

          # save mean to the fist entry
          support<- c(meantmp,pvalue, psupport)
          return(support)
          }


      ###################################################################################
      #### Subsection 17: pvalue supports, under true null  negative binomial      #########
      ###################################################################################

      pvalueByNegativeBinoSupport <- function(marginalandize)
      {

          #  under null that two negative binomials are identical in size and prob, the robinson-smyth test is used
          # y_c | y_c + y_t
          mass1 <- dnbinom(seq(from=0,to=marginalandize[3]),mu=marginalandize[5],size=marginalandize[4])
          mass2 <- dnbinom(seq(from=marginalandize[3],to=0),mu=marginalandize[5],size=marginalandize[4])
          lgh2 = length(mass1)
          #11/21/2011 comment: I do not think the two sizes should be the same

          # produc of two masses, satisfying y_c + y_t = marginal[3]
          prodmass <- mass1*mass2
          if (max(prodmass) ==0) print("Caution joint mass function conditional on total is zero")

         ######## compute p value
           # realized masses
          realmasscontrol <- dnbinom(marginalandize[1],mu=marginalandize[5],size=marginalandize[4])
          realmasstreat <- dnbinom(marginalandize[2],mu=marginalandize[5],size=marginalandize[4])

          realprod <- realmasscontrol*realmasstreat

            # temp[i] is the two-sided pvalue for i-th table given marginals
            pvalue <- sum(prodmass[which(prodmass <= realprod)])/sum(prodmass)
            if (is.na(pvalue) | is.nan(pvalue)) print("pvalue is NA or NaN")

          ###### compute p value support
           temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(prodmass[which(prodmass <= prodmass[i])])/sum(prodmass)
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*prodmass)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)

          # save mean to the fist entry, pvalue to the second
          support<- c(meantmp,pvalue,psupport)
          return(support)
          }


    #######################################################################
    #### Subsection 18:  Storey's and Chen's procedure               #########
    #######################################################################
    # 
      StAndChenFDREstimatorApp <- function (simpvalues,FDRlevel,twopisAtBstPair)
      {
         m <- length(simpvalues)
        # step length in search depends on minimal difference between p values
         sortedpvalvec <- sort(simpvalues)
         diffsortedpvals <- diff(sortedpvalvec)
         # since diff does not compute the differece between minimal p value and zero
         # it should be added, so step length in thresh is the minimum of
         # minimal diff and min p value
         threshstep <- min(min(diffsortedpvals[diffsortedpvals > 0]),min(sortedpvalvec))

        # it is known that FDR threshold is no less than bonferroni, but no larger than
        # FDRlevel. the step length should be sensitive to distances bewteen p values to
        # make practical change in R(t) as t changes
        # sometimes threshstep can be so small that R itself can not handle it

        # Step 1: crude search first   do not start from min(sortedpvalvec) when m is small
        threshveccrude <- seq(from= 0,to=1,by=10^-3)
        ecdfcrude0 <- max(sum(sortedpvalvec <= threshveccrude[1]),1)/m
        # initialize to estimator with the bonferroni a

        chenFDPcrude <- twopisAtBstPair[2]*threshveccrude[1]/ecdfcrude0
        chenFDPcrude0 <- chenFDPcrude

         ### chen's fdr estimator
        jchencrude <- 1
        while (chenFDPcrude <= FDRlevel)
        {
           jchencrude <- jchencrude +1;
           # check if mat threshol (t=1) is reach, is so, break
           if (jchencrude > length(threshveccrude)) {print("max threshold (t=1) in crude search for generalized estimators exhausted");break}
          cheneCDFcrude <- max(sum(sortedpvalvec <= threshveccrude[jchencrude]),1)/m
          chenFDPcrude <- twopisAtBstPair[2]*threshveccrude[jchencrude]/cheneCDFcrude
          if (is.na(chenFDPcrude) | is.nan(chenFDPcrude)) print("NA occured as FDP of generalized estimators")

        }
        if (jchencrude == 2)  {    chensupthreshcrude <- threshveccrude[1]
           chenFDPcrudeAtSThresh <- chenFDPcrude0   }     else {
           chensupthreshcrude <- threshveccrude[jchencrude-1]
        ecdfchensupthreshcrude <- max(sum(sortedpvalvec <= chensupthreshcrude),1)/m
        chenFDPcrudeAtSThresh <- twopisAtBstPair[2]*chensupthreshcrude/ecdfchensupthreshcrude  }
      #  cat("chenFDPcrudeAtSThresh is",chenFDPcrudeAtSThresh,"at sup thresh",chensupthreshcrude,"\n")

        # use while loop for storey's
        istcrude <-1
        stFDPcrude <- twopisAtBstPair[1]*threshveccrude[1]/ecdfcrude0
         stFDPcrude0 <- stFDPcrude
        while (stFDPcrude <= FDRlevel)
        {
          istcrude <- istcrude+1;
          # check if mat threshol (t=1) is reach, is so, break
          if (istcrude > length(threshveccrude)){print("max threshold (t=1) in crude search for Storey's estimators exhausted");break}
          steCDFcrude <- max(sum(sortedpvalvec <= threshveccrude[istcrude]),1)/m
          stFDPcrude <- twopisAtBstPair[1]*threshveccrude[istcrude]/steCDFcrude
          if (is.na(stFDPcrude) | is.nan(stFDPcrude)) print("NA or NaN occured as FDP of Storey's estimators")

        }
        if (istcrude == 2) { stsupthreshcrude <- threshveccrude[1]
             stFDPcrudeAtSThresh <- stFDPcrude0    }     else {
              stsupthreshcrude <- threshveccrude[istcrude-1]
             ecdfstsupthreshcrude <- max(sum(sortedpvalvec <= stsupthreshcrude),1)/m
        stFDPcrudeAtSThresh <- twopisAtBstPair[1]*stsupthreshcrude/ecdfstsupthreshcrude  }
      #  cat("stFDPcrudeAtSThresh is",stFDPcrudeAtSThresh,"at sup thresh",stsupthreshcrude,"\n")

        ############################################
        #### step 2: after a crude threshold has been located, do a finer search from the crude sup threshold
        ### chen's fdr estimator
        # check if mat threshol (t=1) is reach, is so, break
        if (chensupthreshcrude == 1) { print("max threshold (t=1) for generalized estimators exhausted, no fine search in second round")
            chensupthresh <-chensupthreshcrude
        chenFDPAtSThresh <- chenFDPcrudeAtSThresh  }     else {
         threshvecchen <- seq(from =chensupthreshcrude,to=threshveccrude[jchencrude],
                                  by= min((threshveccrude[jchencrude]-chensupthreshcrude)/m,FDRlevel/m))
        ecdf0chen <- max(sum(sortedpvalvec <= threshvecchen[1]),1)/m
        chenFDP <- twopisAtBstPair[2]*threshvecchen[1]/ecdf0chen
        chenFDP0 <- chenFDP
        jchen <- 1
        # use while loop for chen's
        while (chenFDP <= FDRlevel)
        {
          jchen <- jchen+1;
          if (jchen > length(threshvecchen)) break
          cheneCDF <- max(sum(sortedpvalvec <= threshvecchen[jchen]),1)/m
          chenFDP <- twopisAtBstPair[2]*threshvecchen[jchen]/cheneCDF
          if (is.na(chenFDP) | is.nan(chenFDP)) print("NA or nan occured as FDP of generalized estimators in fine search")
          }

        # get the suprema of the two thresholds fro two procedures
        if (jchen==2) { chensupthresh <- threshvecchen[1]
                         chenFDPAtSThresh <- chenFDP0
        } else {     chensupthresh <- threshvecchen[jchen-1]
          ecdfchensupthresh <- max(sum(sortedpvalvec <= chensupthresh),1)/m
        chenFDPAtSThresh <- twopisAtBstPair[2]*chensupthresh/ecdfchensupthresh }
        #  cat("chenFDPAtSThresh is",chenFDPAtSThresh,"at sup thresh",chensupthresh,"\n")
        } # end of else


         # use while loop for storey's
         if (stsupthreshcrude == 1) {  print("max threshold (t=1) for Storey's estimators exhausted, no fine search in second round")
            stsupthresh <- stsupthreshcrude
         stFDPAtSThresh <- stFDPcrudeAtSThresh} else {
        threshvecst <- seq(from = stsupthreshcrude,to=threshveccrude[istcrude],
                               by= min((threshveccrude[istcrude]-stsupthreshcrude)/m,FDRlevel/m))
         ecdf0st <- max(sum(sortedpvalvec <= threshvecst[1]),1)/m
         stFDP <- twopisAtBstPair[1]*threshvecst[1]/ecdf0st
          stFDP0 <- stFDP
         ist <-1
        while (stFDP <= FDRlevel)
        {
          ist <- ist+1;
          if (ist > length(threshvecst)) break
          steCDF <- max(sum(sortedpvalvec <= threshvecst[ist]),1)/m
          stFDP <- twopisAtBstPair[1]*threshvecst[ist]/steCDF
          if (is.na(stFDP) | is.nan(stFDP)) print("NA or NaN occured as FDP for Storey's estimators in fine search")
        }
        if (ist==2) {stsupthresh <- threshvecst[1]
          stFDPAtSThresh <- stFDP0     } else {
        stsupthresh <- threshvecst[ist-1]
         ecfstsupthresh <- max(sum(sortedpvalvec <= stsupthresh),1)/m
        stFDPAtSThresh <- twopisAtBstPair[1]*stsupthresh/ecfstsupthresh }
       #  cat("stFDPAtSThresh",stFDPAtSThresh,"at sup thresh",stsupthresh,"\n")
      } # end of esle


        # get discoveries

        stDiscovery <- which(simpvalues <= stsupthresh)
        if (length(stDiscovery)==0) print("No discoveries found by Storey's procedures")
        chenDiscovery <- which(simpvalues <= chensupthresh)
        if (length(chenDiscovery)==0) print("No discoveries found by generalized estimators")

        StChenDiscoveries <- list(stDiscovery,chenDiscovery,c(stsupthresh,chensupthresh,stFDPAtSThresh,chenFDPAtSThresh))
        return(StChenDiscoveries)

      }

    ###################################################
    ####### Subsection 19: BH FDR procedure        ########
    ######################################################
    # BH procedure requires input of pvalues and their indices
      # BH procedure requires input of pvalues and their indices
      BHFDRApp <- function(Pvals,FDRlevel)
      {   
        if (any(is.na(Pvals)) | any(is.nan(Pvals))) { cat("^^NA or NAN in pvalues... ","\n")
        print(Pvals[is.na(Pvals)])
        }
        
        if (any(is.na(FDRlevel)) | any(is.nan(FDRlevel))) cat("^^NA or NAN FDRlevel... ","\n")
           
       if (is.vector(Pvals))  lgh3 <- length(Pvals)
        if (is.matrix(Pvals) | is.data.frame(Pvals))  lgh3 <- dim(Pvals)[1] 
        PvalAndIdx <- cbind(Pvals,seq(1:lgh3))
        PvalAndIdxOrd <- PvalAndIdx[order(PvalAndIdx[,1]),]
        
     #   cat("^^Dims of PvalAndIdxOrd is:",as.vector(dim(PvalAndIdxOrd)),"\n")
    
        BHstepups <- seq(1:lgh3)*(FDRlevel/lgh3)
        cmp3 <- PvalAndIdxOrd[,1] <= BHstepups
        scmp3 <- sum(cmp3)
    
       # collect rejections if any
        if (scmp3 == 0) {   
        print ("No rejections by BH")        # when there are no rejections
        rejAndTresh <- list(matrix(numeric(0), ncol=2,nrow=0),0)       
        }  else   {  r <- max(which(cmp3))        
            #  cat("^^^ Minimax index in BH is:",r,"\n")
             # cat("^^^ Minimax threshold in BH is:",BHstepups[r],"\n")
                   if (r==1) {
                    # when r =1, a row is chosen, so should change it into a matrix of two columns
                     BHrej = as.matrix(t(PvalAndIdxOrd[1:r,]))  
                   #  print(BHrej)
                   } else {
                      BHrej <- PvalAndIdxOrd[1:r,]
                   }
               rejAndTresh <- list(BHrej,BHstepups[r])    
              }
            return(rejAndTresh)
      }


     ############################################################################################
    ######## Subsection 20: use NBPSeq to get quantities for computing pvalue distribution    ########
    ############################################################################################
         GetPseudoCountsLibsizeMeansSizesFromNBPpack <- function(obj,grp1,grp2) {
            y = obj$pseudo.counts
            mall = obj$pseudo.lib.sizes
            phi = obj$phi
            alpha = obj$alpha
            grp.ids = obj$grp.ids
            n = dim(y)[1]
            rel.counts = y/(matrix(1, n, 1) %*% mall) # each entry of pseudocounts matrix is divided by m to get rel.counts
            grp12.ids = grp.ids %in% c(grp1, grp2)
            grp12.size = sum(grp12.ids)
            pie = rel.counts[, grp12.ids] %*% matrix(1/grp12.size, grp12.size,1)
             # pie is use directly to get the mean, mean is the used to get size in the parametrization
                # matrix(1/grp12.size, grp12.size,1) is a column of 6 entries each being 1/6, this matrix product produces equivalently
           # sum(rel.counts[i,]/6) = sum(rel.counts[i,])/6 as the i-th entry
 
           if (any(pie == 0)) cat("^^NBPSeq: zero estimated frequency obtaiend ...","\n")  
            # warning: zero pie's produced   #del 3/18/2013
    
            grp1.ids = grp.ids %in% grp1
            grp2.ids = grp.ids %in% grp2
            y1 = y[, grp1.ids]
            y2 = y[, grp2.ids]
            s1 = apply(y1, 1, sum)
            s2 = apply(y2, 1, sum)
            r1 = dim(y1)[2]
            r2 = dim(y2)[2]    # r1 = r2 for balanced design
            m1 = mall[grp1.ids][1]
            m2 = mall[grp2.ids][1]
            phi1 = phi2 = phi
            mu1 = pie * m1
            # caution: zero means can be produced
            theta11 = mu1^(2 - alpha)/phi1  # caution: zero size can be produced
    
            p1 = mu1/(mu1 + theta11)  # caution: because of zero mu1 and theta11, 0/0 division occurs
    
            theta1 = theta11 * r1
            mu2 = pie * m2     # mu2 == mu1 as vectors, when r1 = r2
            theta21 = mu2^(2 - alpha)/phi2  #  theta21 ==  theta11, when r1 = r2
            p2 = mu2/(mu2 + theta21)
    
            theta2 = theta21 * r2  #   theta2 == theta1 as vectors, when r1 = r2
           
           # ready to return
           quantitiesForPvalDists <- list(s1,s2,theta1,p1) # under balanced desing theta1=theta2, p1=p2
            return(quantitiesForPvalDists)
        }
      
    ###################################################################################
    #### subsection 21: pvalue supports, under true null  negative binomial      #########
    ###################################################################################
     pvalueByNegativeBinoSupportApp <- function(marginalandize)
      {
        # check line
        if (any(is.na(marginalandize)) | any(is.nan(marginalandize)))  cat("^^NA or NAN found in a row of input ...","\n")
        # end check line
        if (length(marginalandize) != 4) cat("^^ Length of processed row as input should be 4 ..","\n")
        s1 <- marginalandize[1]
        s2 <- marginalandize[2]
      
      # under the null, two lines changed 3.21.2013
        s_com <- marginalandize[3]; p_com <- marginalandize[4]
        
        #### embed comput pvalue function from NBP package, except otherwise noted   
        # start commenting   
      #  pr.obs = dnbinom(s1, size = theta1, prob = 1-p1) * dnbinom(s2, size = theta2, prob = 1-p2);
      #  s = s1 + s2;
        # from NBPSeq: compute probs under null, this is wrong unless, theta2=theta1, p1 = p2.
        # NBPSeq gives p1=p2, mu1=mu1 for balanced desing, since pie is pulled and library size is the same. 3/18/2013  
        
      #  pr = dnbinom(0:s, size = theta1, prob = 1-p1) * dnbinom(s:0, size = theta2, prob = 1-p2);   
        # end of commenting
                   
          pr.obs = dnbinom(s1, size = s_com, prob = 1-p_com) * dnbinom(s2, size = s_com, prob = 1-p_com)
          s = s1 + s2      
          pr = dnbinom(0:s, size = s_com, prob = 1-p_com) * dnbinom(s:0, size = s_com, prob = 1-p_com)   
             
       # cat("^^ product probs are:",pr,"\n")
         if (is.na(pr.obs) | is.nan(pr.obs)) cat("^^Observed prob is NA or NAN...","\n")
       
         if (any(is.na(pr))) cat("^^NA in product probs inside pval compute-function...","\n")
           if (any(is.nan(pr))) cat("^^NAN in product probs inside pval compute-function...","\n")
        
        if ((sum(pr)) == 0 & !any(is.na(pr)) & !any(is.nan(pr))) 
          cat("^^Sum of joint probs is zero inside pval compute-function...","\n")
        
        id.extreme = (pr <= pr.obs);
      
        ## thetas = theta1 + theta2;
        ## ps = dnbinom(s, thetas, 1-p);        
        pval = sum(pr[id.extreme])/sum(pr)
        
        if (is.na(pval)) cat("^^NA in pvalue obtained inside pval compute-function...","\n")
           if (is.nan(pval)) cat("^^NAN in pvalue obtained compute-function...","\n")  
      
          lgh2 = length(pr)
      
          ###### compute p value support
           temp <- double(lgh2)
          for (i in 1:lgh2)
          {
            # temp[i] is the two-sided pvalue for i-th table given marginals
            temp[i] <- sum(pr[pr <= pr[i]])/sum(pr)
          }
          # compute the expectation of pvalue
           meantmp <- sum(temp*pr)
           # sort pvalue support
          psupport <- sort(temp,decreasing=FALSE)
      
          # save mean to the fist entry, pvalue to the second
          support<- c(meantmp,pval,psupport)
          return(support)
          }
    
##############################################################
#  subsection 22: function to extract normalization process from edgeR
##############################################################
CountsNormalizerInEdgeR <- function (object, pair = 1:2, dispersion = NULL, prior.count.total = 0.5)
{
   cat("^^edgeR: start normalizing ...","\n")
   
    if (!is(object, "DGEList"))
        stop("Currently only supports DGEList objects as the object argument.")
    if (length(pair) != 2)
        stop("Pair must be of length 2.")
    group <- as.factor(object$samples$group) 
    levs.group <- levels(group)
    if (is.numeric(pair))  {
        pair <- levs.group[pair] } else  { pair <- as.character(pair)}
        
    if (!all(pair %in% levs.group))
        stop("At least one element of given pair is not a group.\n Groups are: ",
            paste(levs.group, collapse = " "))

    # change: only allow CommonDispersion, but it seems mglmOneGroup can take vector of dispersion 3/20/2013
    
    ## start of normalizing counts, codes from exactTest  ##
    ldisp <- length(dispersion)
    ntags <- nrow(object$counts)
    if (ldisp == 1)
        dispersion <- rep(dispersion, ntags)
    group <- as.character(group)
    j <- group %in% pair
    y <- object$counts[, j, drop = FALSE]
    lib.size <- object$samples$lib.size[j]
    norm.factors <- object$samples$norm.factors[j]
    group <- group[j]
    if (is.null(rownames(y)))
        rownames(y) <- paste("tag", 1:ntags, sep = ".")
    lib.size <- lib.size * norm.factors
    offset <- log(lib.size)
    lib.size.average <- exp(mean(offset))  
    
    abundance <- mglmOneGroup(y, dispersion = dispersion, maxit = 50, offset = offset)  # default maxit 50
        cat("^^EdgeR: length of abundance:",length(abundance),"\n")   # line added 3/19/2013
        
    logCPM <- (abundance + log(1e+06))/log(2)
    prior.count <- lib.size
    prior.count <- prior.count.total * prior.count/sum(prior.count)
    j1 <- group == pair[1]
    n1 <- sum(j1)
    if (n1 == 0)
        stop("No libraries for", pair[1])
    y1 <- y[, j1, drop = FALSE]
    abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], ntags,
        n1, byrow = TRUE), offset = offset[j1])
    j2 <- group == pair[2]
    n2 <- sum(j2)
    if (n1 == 0)
        stop("No libraries for", pair[2])
    y2 <- y[, j2, drop = FALSE]
    abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], ntags,
        n2, byrow = TRUE), offset = offset[j2])
    logFC <- (abundance2 - abundance1)/log(2)
   
    e <- exp(abundance)      # here for input.mean
    
    input.mean <- matrix(e, ntags, n1)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j1])
    y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean,
        dispersion = dispersion)
    input.mean <- matrix(e, ntags, n2)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * lib.size[j2])
    y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean,
        dispersion = dispersion)
    # end of counts normalization
    
    ### sart of codes from exactTestDoubleTail in edgeR
    ntags <- NROW(y1)
    n1 <- NCOL(y1)
    n2 <- NCOL(y2)
    if (n1 > 1)
        s1 <- round(rowSums(y1))
    else s1 <- round(y1)
    if (n2 > 1)
        s2 <- round(rowSums(y2))
    else s2 <- round(y2)
    if (length(dispersion) == 1)
        dispersion <- rep(dispersion, ntags)
    s <- s1 + s2
    mu <- s/(n1 + n2)
    mu1 <- n1 * mu
    mu2 <- n2 * mu
    ### end of codes from exactTestDoubleTail in edgeR

    # in order to get the conditional test, the sizes of the negative binomial distribution
    # should be identical, so it should be this when n1=n2
    # the original codes in function "exactTestDoubleTail" of edgeR sets size= size = (n1 + n2)/dispersion
    # for joint, but size1 <- n1/dispersion,  size2 <- n2/dispersion for the NB for each group
    
    size_com = 1/dispersion; size_gp1 = n1/dispersion;   p_com = mu/(mu + size_com) # added 3/21/2013
    # above line: because of the way  pvalueByNegativeBinoSupportApp, size_grp_com and pi_com are used
                               
    edgeR_norm_counts = cbind(s1,s2,size_gp1,p_com)  # changed, now 4 cols
    
    return(edgeR_norm_counts)

}  #end of extraction function
#### end of subsections


  ############################################################# ###########################
      #### Subsection 23: Function to Get marginals and cellcounts from generated counts
      #####             for differential expression
      ################################################################# ###########################

 getcellcountsandmarginals_DE <- function(data_in)
      {
        m <-nrow(data_in)
        # construct 3-3 full table
        tables<- array(0, dim=c(3,3,m))

        tables[1,1,] = data_in[,1]
        tables[2,1,] = data_in[,2]
        tables[1,2,] = data_in[,3] - data_in[,1]
        tables[2,2,] = data_in[,4] - data_in[,2]
        tables[1,3,]= data_in[,3]
        tables[2,3,]= data_in[,4]
        tables[3,3,] = data_in[,3] + data_in[,4]
        tables[3,1,] = data_in[,1] + data_in[,2]
        tables[3,2,] = tables[3,3,] - tables[3,1,]
        
        simallcellcounts <- vector("list",m)
        for (i in 1:m) {simallcellcounts[[i]] <- tables[1:2,1:2,i]}

        simallmarginals <- vector("list",m)
        for (i in 1:m) {simallmarginals[[i]] <- c(tables[1,3,i],tables[2,3,i],tables[3,1,i])}

         y2 <- list(simallcellcounts,simallmarginals)
          return(y2)
      }
