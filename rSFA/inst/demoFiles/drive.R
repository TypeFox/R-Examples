# % This demo reproduces the driving force example of the logistic map (Figure 2)
# % from Wiskott, L. (2003), "Estimating Driving Forces of Nonstationary Time
# % Series with Slow Feature Analysis",
# % http://arxiv.org/abs/cond-mat/0312317. 
# %
# % Here the driving force gam(t) consist of a slow part gamS(t) and a fast
# % part gamF(t). See Konen, W. (2009), "How slow is slow? SFA detects signals
# % that are slower than the driving force", http://arxiv.org/abs/0911.4397. 
# %  
# % INPUT:
# %       m               embedding dimension (def. 19)
# %       tau             delay (def. 1)
# %       nuf             base frequency, the larger the faster varies the 
# %                       driving force (def. 50)
# %       q               'chaos parameter': for q<0.3 the logistic map is
# %                       fully in the chaotic regime, for 0.3<q<3.9 there
# %                       are more and more non-chaotic parts in the time
# %                       series (def. 0.1)
# %       dographics      (def. 2) if >= 2, make some plots
# %       noiseperc       (def. 0) e.g. 0.05 means 'add 5% noise to the data w
# %       method          (def. 'SVDSFA') either 'SVDSFA', the new, numeric stable  
# %                       method, or 'GENEIG' (the old version of [Berkes03]). 
# %						  GENEIG is not implemented in the current version, since
# %						  R lacks the option to calculate generalized eigenvalues easily.
# %       mf              (def. 1) multiplication factor for modulation
# %                       frequeny gamS. If mf<1 the slow component gets
# %                       slower while the fast component stays the same
# % OUTPUT: Struct res with tags:
# %       cor             cor(1): correlation between gam(t) and y1(t)
# %                       cor(2): correlation between gamS(t) and y1(t)
# %       slowness        slowness for gam(t), a*y1(t)+b, original w(t),gamS(t)
# %       eta             alternative slowness
# %
# % EXAMPLE:
# %       drive1(19,1,50,1.9,2,0,method)
# % This will result in a wrong slow signal for method = 'GENEIG', but to a
# % well aligned slow signal for method = 'SVDSFA'.
# %
# % Author: Wolfgang Konen, May 2009 - Dec 2009
# %
drive<-function(m=19,tau=1,nuf=59,q=0.1,dographics=2, noiseperc=0.00,method="SVDSFA",mf=1){ 

    nlReg = 1;             #% if =1, do non-linear regression for the functions 
                            #% gam and gamS and show results in figure(3)
                            #% report in res.etaR the eta-value of the
                            #% fitted function
    ppType = "PCA";        #% preprocessing type
    
    #% create the input signal
    T0 = 5990;
    T = T0+m*tau;
    t0 = seq(0,T0-1,length=T0)
    t = seq(0,T-1,length=T)

                    #% REMARKABLE: If we set nuf=50, then the driving force is
                    #% no longer very 'slow'. However, it contains a slower
                    #% subcomponent. And - bingo - SFA detects the slower 
                    #% subcomponent and thus a signal which is more slowly 
                    #% than the driving force itself (has a smaller eta 
                    #% than gam(t), see below).
                    #% [nuf=50 requires however m=20 or higher to reconstruct 
                    #% the slow component with good accuracy]
    #%
    #% Variant 1: driving force gam is composed of two sine waves at
    #% frequency nuf*0.0047 and nuf*0.0005 (this is what we wrongly called 'beat'
    #% above). SFA nicely learns to detect either one of these components,
    #% depending on nuf
    #%gam = sin(nuf*0.0021*t).*cos(nuf*0.0026*t);        % 0.0026
    gam = (sin(nuf*0.0047*t) + sin(nuf*0.0005*t*mf))/2;    #% equivalent to line above for mf=1
    #%gam = (0*sin(nuf*0.0047*t) + sin(nuf*0.0005*t))/2;     % activate this line for arxiv2009_SFA2.m
    #%gam = (sin(nuf*0.0012*t) + sin(nuf*0.0005*t))/2;     
    gamS = + sin(nuf*0.0005*t*mf);  #% only the slow component
    #gamF = + sin(nuf*0.0047*t);     #% only the fast component
    #%
    #% Variant 2: driving force is composed of two nearby sine waves at
    #% frequency nuf*0.0021 and nuf*0.0026. Now this contains a modulation (beat) at
    #% frequency nuf*0.0005, but the modulation is hardly detected by SFA, no matter how
    #% high we put nuf.
    #%    gam = (sin(nuf*0.0021*t) + sin(nuf*0.0026*t))/2; 
    #%    gamS = + cos(nuf*0.0005*t);     % only the slow component
    #%    gamF = + sin(nuf*0.0047*t);     % only the fast component

    w = matrix(0,T,1)+0.1;
#%     ntrans=30;
#%     for i=1:ntrans
#%         %avoid transient signal
#%         w(i+1)=(4.0-q + gam(i)*0.1)*w(i)*(1-w(i));     % logistic map, Eq. (22)
#%     end    
#%     w(1:2)=w((ntrans-1):ntrans);
    for(i in 1:(T-1)){
        w[i+1]=(4.0-q + gam[i]*0.1)*w[i]*(1-w[i]);     #% logistic map, Eq. (22)
    }
    if (dographics>=2){ #first plot
       # % plot the input signal
        #if (~BIOMA), clf; end;
        pskip=2;
        ende=min(c(6000, customSize(t,2)));
        par(mfrow=c(3,1))  ; 
        plot(t[seq.int(from=1,by=pskip,to=ende)],w[seq.int(from=1,by=pskip,to=ende)],ylab="w(t)",xlab="t",col="red",pch=".",main="input signal w(t)");
        #% why '1:pskip:end'? - to make the PDF of the plot with pskip=2
        #% somewhat smaller (360kB instead of 600kB). Default is pskip=1.
    }    
    noise=rnorm(T)*noiseperc;
	print(paste("noise in % = ", mean(abs(noise))/mean(abs(w))*100))
    w=w+noise;

    x = matrix(0,T0,m);
    for(i in 1:m){
        low = 1+tau*(i-1); 
        hig = T0+tau*(i-1);
        x[,i]=w[low:hig];
	}

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%
    #% Slow Feature Analysis
    #%
    sfaList = sfa2(x,method,ppType);
	y=sfaList$y
	sfaList$y<-NULL
	
    #% plot the output of the slowest varying function
    y1 = y[,1];
    y2 = y[,2];
   
    if(dographics>=2){ #second plot
        plot(t0[seq.int(from=1,by=pskip,to=ende)], y1[seq.int(from=1,by=pskip,to=ende)],col="blue",pch=".",ylab="y1(t)",xlab="t",main="output of the slowest varying function y1(t)");
        #% why '1:pskip:end'? - to make the PDF of the plot with pskip=2
        #% somewhat smaller (360kB instead of 600kB). Default is pskip=1.
    }
    
    gm = t(gam[(1+floor(m/2)):(T0+floor(m/2))]);     #% why '+m/2'? If we pick up the driving 
                                            #% force here, there is no phase
                                            #% shift to the slowest SFA signal
    gS = t(gamS[(1+floor(m/2)):(T0+floor(m/2))]);
    #gF = t(gamF[1+floor(m/2):T0+floor(m/2)]); #gF or g_F is never used anywhere

    
    #% find the linear transformation [a b] which brings a*y(t)+b and gm(t) in
    #% optimal alignment: z = [a b] = X\gm is the vector which solves the 
    #% linear system X*z = gm optimally in the LS sense.
    X = cbind(y1, matrix(1,rev(customSize(y1)))); #TODO the process here might be better solved with lm ?
	z=sfaBSh(X,t(gm)); #user defined function in sfaHelperFunctions
	a=z[1]; b=z[2];
    z0=sfaBSh(X,t(gS)); 
	a0=z0[1]; b0=z0[2];
    
    if(dographics>=2){ #third plot
		#% Now plot a*y(t)+b as blue points and overlay gm(t) as green line:
		ende=min(c(500,length(gm),length(t0)));
		pskip3=4;
		plot(t0[1:ende],gm[1:ende],col="green",type="l",xlab="t",ylab="gamma(t)",
			main=expression(paste("true driving force ",
						gamma,
						"(t) [green], its slow part ",
						gamma[S],
						"(t) [red], alignment of y1(t) [blue]",
						sep="")),
			cex.main=1.2);
		lines(t0[1:ende],gS[1:ende],col="red",lty=2)
		lines(t0[seq.int(from=1,by=pskip3,to=ende)], a0*y1[seq.int(from=1,by=pskip3,to=ende)]+b0,col="blue")			

    }  #%% dographics

    if (nlReg==1){ 
        R = sfaNlRegress(sfaList,x,gm)$R;
        RS = sfaNlRegress(sfaList,x,gS)$R;
        if (dographics>=2){ #additional new plot
			ende=min(c(500,length(gm),length(t0)));
            dev.new()
            plot(t0[1:ende],gm[1:ende],col="red",type="l",
				main=expression(paste("NL regression to true driving force ",
						gamma,
						"(t) [black/red], and to slow part ",
						gamma[S],
						"(t) [blue/red]",
						sep="")),
				cex.main=0.9);
            lines(t0[1:ende],R[1:ende],col="black");    
            lines(t0[1:ende],gS[1:ende],col="red");  
            lines(t0[1:ende],RS[1:ende],col="blue");    
        }
		etaR=0
        etaR[1] = etaval(R,T0);        #% slowness of NL regression for gamma(t)
        etaR[2] = etaval(RS,T0);      #% slowness of NL regression for gamma_S(t)
    }
    #% correlation and slowness as indicative numbers:
    #% (slowness = mean of |finite differences| )
#	corr=cor(gm,y1); #TODO doesnt work yet
#	coefc[1]=corr[1,2];
#	corr=cor(gS,y1); 
#	coefc[2]=corr[1,2];
    # if dographics>=3
        ##% Plot |correlation| and slowness for *all* SFA output signals,
        ##% ordered by their eigenvalue size (Where are sudden jumps?)
        # C = abs(corrcoef([gS,y])); 
        # figure(4), subplot(2,1,1); plot(C(1,2:end));
        # for i=1:size(y,2)
            # etay(i)=etaval(y(:,i),T0);
        # end
        # subplot(2,1,2); plot(etay);
    # end
#    if(is.na(coefc[1])){        #% should not happen (but can, if nuf=0)
#        warning("correlation with gm contains NaN");
#    }
#    if(is.na(coefc[2])){        #% should not happen (but can, if nuf=0)
#        warning("correlation with gS contains NaN");
#    }
	slowness=0
    slowness[1] = mean(abs(diff(gm)));
    slowness[2] = mean(abs(diff(a*y1+b)));
    slowness[3] = mean(abs(diff(w)));
    slowness[4] = mean(abs(diff(gS)));
	eta=0
    eta[1] = etaval(gm,T0); #%sqrt(mean(diff(gm).^2))/sqrt(mean(gm.^2));
    eta[2] = etaval(a*y1+b,T0); #%sqrt(mean(diff(a*y1+b).^2))/sqrt(mean((a*y1+b).^2));
    eta[3] = etaval(w,T); #%sqrt(mean(diff(w).^2))/sqrt(mean(w.^2));
    eta[4] = etaval(gS,T0); #%sqrt(mean(diff(gS).^2))/sqrt(mean(gS.^2));
    res=list()
#    res$cor=coefc;
    res$eta=eta;
    if(nlReg==1){res$etaR=etaR;}
    res$slowness=slowness;
    res$gm=gm; 
    res$gS=gS;
    res$t0=t0;
    res$y1=y1;
    res$sfaList=sfaList;
	return(res)
}
