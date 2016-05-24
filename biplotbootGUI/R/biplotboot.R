

biplotboot <-function(x)
{
        #############################################################################
        #########	libraries
        #############################################################################
        
        transforma <- function(tipo, matriz)
        {    
                if (tipo=="Subtract the global mean"){                        
                        Xstd <- as.matrix(matriz)
                        media <- mean(Xstd)
                        Xstd <- Xstd-media        	
                }#end if (tipo=="Subtract the global mean")
                
                if (tipo=="Column centering"){
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 2, function(x){x-mean(as.matrix(x))})
                }#end if (tipo=="Column centering")
                
                if (tipo=="Standardize columns"){	
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 2, function(x){(x-mean(as.matrix(x)))/sqrt(var(as.matrix(x)))})
                }#end if (tipo=="Standardize columns")
                
                if (tipo=="Row centering"){		
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 1, function(x){x-mean(as.matrix(x))})
                        Xstd <- t(Xstd)
                }#end if (tipo=="Row centering")
                
                if (tipo=="Standardize rows"){
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 1, function(x){(x-mean(as.matrix(x)))/sqrt(var(as.matrix(x)))})
                        Xstd <- t(Xstd)
                }#end if (tipo=="Standardize rows")
                
                if (tipo=="Double centering"){
                        Xstd <- as.matrix(matriz)
                        mediac <- colMeans(Xstd)
                        mediaf <- rowMeans(Xstd)
                        globalm <- mean(Xstd)
                        
                        mediafm <- array(unlist(rep(rowMeans(Xstd), dim(Xstd)[2])), dim=dim(Xstd))
                        mediacm <- t(array(unlist(rep(colMeans(Xstd), dim(Xstd)[1])), dim=c(dim(Xstd)[2],dim(Xstd)[1])))
                        Xstd <- Xstd - mediafm - mediacm +globalm
                        
                }#end if (tipo=="Double centering")
                
                if (tipo=="Raw data"){
                        Xstd <- as.matrix(matriz)
                }#end if (tipo=="Raw data")
                
                rownames(Xstd) <- rownames(matriz)
                return(Xstd)        
                
        }# end transforma 
        
        
        ex.biplot <- function(Xpon, tipo, nejes, tChoice, tindi, tvar)
        {
                ##############################################################################
                #####        	Coordinates
                ##############################################################################
                ejes<-c()
                for (i in 1:nejes)
                {
                        ejes<-c(ejes, paste("Axis",i))
                }#end for (i in 1:nejes)
                
                Xpon<-transforma(tChoice,Xpon)
                descom <<- La.svd(Xpon)
                sumaRvalprop<-sum((descom$d)^2)	
                inerciatot<<-(descom$d[1:length(descom$d)])^2/sumaRvalprop
                descom<-La.svd(Xpon,nu=nejes,nv=nejes)
                bonajuste<-0
                if (tipo == "RMP"){
                        coindividuos<-descom$u%*%diag(descom$d[1:nejes])
                        covariables<-t(descom$v)				
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        suma2valprop<-sum((descom$d[1:nejes])^2)
                        inercia<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<-cumsum(inercia)		
                        bonajuste<-(suma2valprop/sumaRvalprop)*100
                        calcol<-nejes/length(inerciatot)*100
                        calfilas<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo == "RMP") 
                
                if (tipo == "CMP"){
                        coindividuos<-descom$u
                        covariables<-t(descom$v)%*%diag(descom$d[1:nejes])				
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        suma2valprop<-sum((descom$d[1:nejes])^2)
                        inercia<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<-cumsum(inercia)
                        bonajuste<-(suma2valprop/sumaRvalprop)*100
                        calfilas<-nejes/length(inerciatot)*100
                        calcol<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo == "CMP")
                
                if (tipo=="RCMP")
                {
                        coindividuos<-descom$u%*%diag(descom$d[1:nejes])
                        covariables<-t(descom$v)%*%diag(descom$d[1:nejes])
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        suma2valprop<-sum((descom$d[1:nejes])^2)
                        inercia<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<-cumsum(inercia)
                        calcol<-(suma2valprop/sumaRvalprop)*100
                        calfilas<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo=="RCMP")
                
                covartotal<<-covariables		
                coindividuosnam<-as.data.frame(coindividuos)
                rownames(coindividuosnam)<-tindi
                colnames(coindividuosnam)<-ejes
                
                covariablesnam<-as.data.frame(covariables)
                rownames(covariablesnam)<-tvar
                colnames(covariablesnam)<-ejes
                
                coindivcuad<-coindividuos^2
                CRTi<-rowSums(coindivcuad)
                CRTi<-(CRTi*1000)/suma2valprop
                CRTi<-as.data.frame(CRTi)
                rownames(CRTi)<-rownames(coindividuosnam)
                
                covarcuad<-covariables^2
                CRTj<-rowSums(covarcuad)
                CRTj<-(CRTj*1000)/suma2valprop
                CRTj<-as.data.frame(CRTj)
                rownames(CRTj)<-rownames(covariablesnam)
                
                ##############################################################################
                #####		Length of variables
                ##############################################################################
                longitudprin<<-array(dim=dim(CRTj))
                
                for (i in 1:dim(covartotal)[1])
                {
                        lonprin<-sqrt(covartotal[i,1]^2 + covartotal[i,2]^2)
                        longitudprin[i,1]<<-lonprin
                }#end for (i in 1:dim(covartotal)[1])
                
                longitudprin <<- as.data.frame (longitudprin)
                rownames(longitudprin)<-rownames(CRTj)
                
                CREiFq<-array(dim=dim(coindividuos))	
                CREjFq<-array(dim=dim(covariables))
                sumavar<-rowSums(covarcuad)
                CRFqEi<-coindivcuad
                sumaindi<-rowSums(coindivcuad)
                CRFqEj<-covarcuad
                
                for(i in 1:nejes)
                {
                        CREiFq[,i]<-((coindivcuad)[,i]*1000)/((descom$d[i])^2)
                        CREjFq[,i]<-((covarcuad)[,i]*1000)/((descom$d[i])^2)			  
                        CRFqEi[,i]<-((coindivcuad)[,i]*1000)/(sumaindi)
                        CRFqEj[,i]<-((covarcuad)[,i]*1000)/(sumavar)
                }#end for(i in 1:nejes)
                
                CREiFq<-as.data.frame(CREiFq)
                rownames(CREiFq)<-rownames(coindividuosnam)
                colnames(CREiFq)<-ejes
                
                CREjFq<-as.data.frame(CREjFq)
                rownames(CREjFq)<-rownames(covariablesnam)
                colnames(CREjFq)<-ejes
                
                CRFqEi<-as.data.frame(CRFqEi)
                rownames(CRFqEi)<-rownames(coindividuosnam)
                colnames(CRFqEi)<-ejes
                
                CRFqEj<-as.data.frame(CRFqEj)
                rownames(CRFqEj)<-rownames(covariablesnam)
                colnames(CRFqEj)<-ejes
                
                ######c?lculo de ?ngulos entre variables
                n <- dim(covartotal)[1]
                o <- c(0,0)
                ang <- array(dim = c(n,n)) 
                
                for (z in 1:n)
                {
                        for(j in z:n)
                        {
                                ang[z,j] <- (atan2(covartotal[z,1], covartotal[z,2]) - atan2(covartotal[j,1], covartotal[j,2])) * 180/pi
                                if (ang[z,j] < 0)
                                        ang[z,j] <- -ang[z,j]
                                if (ang[z,j] > 180)
                                        ang[z,j] <- -(ang[z,j] - 360)
                                ang[j,z] <- ang[z,j]
                        }#end for(j in z:n)
                }#end for (z in 1:n)
                
                
                ang<-as.data.frame(ang)
                rownames(ang)<-rownames(CRTj)
                colnames(ang)<-rownames(CRTj)
                
                ######c?lculo de ?ngulos entre variables y ejes
                n <- dim(covartotal)[1]
                o <- c(0,0)
                axisx <- c(1,0)
                axisy <- c(0,1)
                angax <- array(dim = c(n,2)) 
                
                for (z in 1:n)
                {
                        angax[z,1] <- (atan2(covartotal[z,1], covartotal[z,2]) - atan2(axisx[1], axisx[2])) * 180/pi
                        angax[z,2] <- (atan2(covartotal[z,1], covartotal[z,2]) - atan2(axisy[1], axisy[2])) * 180/pi
                        
                        for (j in 1:2)
                        {   
                                if (angax[z,j] < 0)
                                        angax[z,j] <- -angax[z,j]
                                if (angax[z,j] > 180)
                                        angax[z,j] <- -(angax[z,j] - 360)
                                if (angax[z,j] > 90)
                                        angax[z,j] <- -(angax[z,j] - 180)
                        }#end for (j in 1:2)
                }#end for (z in 1:n)
                
                
                angax<-as.data.frame(angax)
                rownames(angax)<-rownames(CRTj)
                colnames(angax)<-c("Axis 1", "Axis 2")
                
                resultados<- list(ejes=ejes,descom=descom, coindividuos=coindividuos, covariables=covariables, suma2valprop=suma2valprop,
                                  inercia=inercia, cuminer=cuminer, bonajuste=bonajuste, calcol=calcol, calfilas=calfilas, covartotal=covartotal,
                                  coindividuosnam=coindividuosnam, covariablesnam=covariablesnam, CRTi=CRTi, CRTj=CRTj, longitudprin=longitudprin,
                                  CREiFq=CREiFq, CREjFq=CREjFq, CRFqEi=CRFqEi, CRFqEj=CRFqEj, ang=ang, angax=angax, valpro=descom$d)
                return(resultados)
        }# function ex.biplot
        
        
        resample_boot<-function(X)
        {
                indices <- sample(1:dim(X)[1], replace = T)
                Xres <- X[indices,]
                return(list(Xres, indices))
        }# end resample_boot<-function(X)
        
        cal.ic <- function (muestra, liminf, limsup, valorobs, muestrajack, niter)
        {
                c.mean <- mean (muestra)
                se <- sd(muestra)
                sesgo <- c.mean - valorobs
                t.ic <- se * (-qt(liminf,(length(muestra)-1)))
                ic.t <- c(c.mean - t.ic, c.mean + t.ic)
                ic.p <- quantile (muestra,c(liminf, limsup), na.rm=TRUE)
                z0 <- qnorm(length(muestra[which(muestra<valorobs)])/as.numeric(niter))
                dent <- mean(muestrajack)- muestrajack
                acc <- sum(dent * dent * dent)/(6 * (sum(dent * dent))^1.5)
                alpha1 <- qnorm(liminf)
                alpha2 <- qnorm(limsup)
                zalpha1 <- pnorm(z0 + (z0 + alpha1)/(1 - acc * (z0 + alpha1)))
                zalpha2 <- pnorm(z0 + (z0 + alpha2)/(1 - acc * (z0 + alpha2)))
                ic.bca <- quantile (muestra,c(zalpha1, zalpha2), na.rm=TRUE)
                return(c(c.mean, se,sesgo,ic.t, ic.p, ic.bca))
        }#end cal.ic <- function (muestra, liminf, limsup, valorobs)
        
        tclRequire("BWidget")
        
        mientorno <- new.env()
        cpdfVal<-"Color pdf"
        cepsVal<-"Color eps"
        Xpon <- NULL
        symbols <- c("*",".", "o","O","0","+","-","|","%","#")
        tipo<-"RCMP" 
        nejes<-3  
        dim1<-1
        dim2<-2
        dim1ant<-1
        dim2ant<-2
        dim3<-3 
        niter<-100
        alphaic <- 95
        indicei<-NULL
        Namei<-NULL
        Cexi<-1
        NameCexi<-NULL
        colori<-NULL
        simChoicei<-NULL
        colores<-c()  
        indicev<-NULL
        Namev<-NULL
        NameValv<-NULL
        Cexv<-1
        NameCexv<-NULL  
        colorv<-NULL
        colorder<-NULL
        Nameder<-NULL
        Cexder<- NULL
        descom<-NULL
        inerciatot<-NULL
        msginertia<-NULL
        nejes<-NULL
        cbVal<-NULL
        cgfVal <- NULL
        ccrVal <- NULL
        cciVal <- NULL
        ceiVal <- NULL
        cavarVal <- NULL
        cavarjVal <- NULL
        #cvcooVal <- NULL
        ccrtiVal <- NULL
        ccreifqVal <- NULL
        ccrfqeiVal <- NULL
        ccrtjVal <- NULL
        ccrejfqVal <- NULL
        ccrfqejVal <- NULL
        clvVal <- NULL
        anteriorx <- NULL
        anteriory <- NULL
        xCoords <- NULL
        yCoords <- NULL
        zCoords <- NULL
        datos <- NULL
        textos <- NULL
        indexClosest <- NULL
        indexLabeled <- NULL                           
        indexLabeledaux <- NULL
        parPlotSize <- NULL
        usrCoords <- NULL
        tChoice <- "Raw data"
        img <- NULL
        imgbar <- NULL
        ejes <- NULL
        descom <- NULL
        coindividuos <- NULL
        covariables <- NULL
        suma2valprop <- NULL
        inercia <- NULL
        cuminer <- NULL
        bonajuste <- NULL
        calcol <- NULL
        calfilas <- NULL
        covartotal <- NULL
        coindividuosnam <- NULL
        covariablesnam <- NULL
        CRTi <- NULL
        CRTj <- NULL
        longitudprin <- NULL
        CREiFq <- NULL
        CREjFq <- NULL
        CRFqEi <- NULL
        CRFqEj <- NULL
        ang <- NULL
        angax <- NULL
        clb <- "normal"
        sumaRvalprop <- NULL
        
        muestras <- c()
        muestrasstd <- c()
        covar <- c()
        covarb <- c()
        coind <- c()
        angulo <- c()
        anguloeje <- c()
        autovalores <- c()
        CRTjv <- c()
        CREjFqv <- c()
        CRFqEjv <- c()
        bondades <- c()
        calidadcolumnas <- c()
        indices <- NULL
        proj<-"normal"
        Choiceproj<- 0
        colvariablesp <- c()
        tit_graph <- "Graph"
        nclust <-"3"
        niterclust <- "10"
        nstart <- "1"
        Limix1 <-tclVar("0")
        Limix2 <-tclVar("0")
        Limiy1 <-tclVar("0")
        Limiy2 <-tclVar("0")
        hescale <- "1.5"
        vescale <- "1.5"
        textvariables <- NULL
        textindividuos <- NULL
        pertb <- NULL
        pert <- NULL
        dcolor <- NULL
        hc <- NULL
        pertcentr <- NULL
        distchoos <- NULL
        linkchoos <- NULL
        algorchoos <- NULL
        pamx <- NULL
        entry.limix1 <- NULL
        entry.limix2 <- NULL
        entry.limiy1 <- NULL
        entry.limiy2 <- NULL
        
        
        colorescoor <- c("skyblue","red","green","blue","yellow","pink","orange", "navyblue",
                         "violet", "brown", "grey", "navyblue", "darkgreen", "papayawhip", "paleturquoise", "purple",
                         "seagreen", "azure", "coral", "springgreen", "steelblue", "plum", "orchid", 
                         "lemonchiffon", "lavender", "honeydew", "gold", "deeppink", "darksalmon", "darkmagenta")
        
        ##############################################################################
        ####	We create vectors of the colors
        ##############################################################################
        if(!is.data.frame(x))
        {
                msg<-("ERROR: data must be a data frame.")
                tkmessageBox(message=msg)
                stop()
        }
        
        if(dim(x)[1]!=dim(na.omit(x))[1])
        {
                msg<-("ERROR: data have missing values. Analysis can not be executed.")
                tkmessageBox(message=msg)
                stop()
        }
        colvariables<-rep("blue",times = dim(x)[2])		
        colindividuos<-rep("green",times = dim(x)[1])
        
        
        ##############################################################################
        ####	We create vectors of the colors
        ##############################################################################
        
        textvariables<-colnames(x)
        textindividuos<-rownames(x)
        
        ##############################################################################
        ####	We create vectors of the character size
        ##############################################################################
        
        cexvariables<-rep(1,times = dim(x)[2])		
        cexindividuos<-rep(1,times = dim(x)[1])
        
        
        ##############################################################################
        ####	We create vectors of the symbols
        ##############################################################################
        
        simvariables<-rep(" ",times = dim(x)[2])		
        simindividuos<-rep("+",times = dim(x)[1])
        
        
        
        #         for (i in 1:dim(x)[1])
        #         {
        #                 colindividuos[i]<-paste("colori",i, sep = "")
        #                 assign(colindividuos[i],"black", envir = mientorno)
        #                 textindividuos[i]<-paste("labeli",i, sep = "")
        #                 assign(textindividuos[i],rownames(x)[i], envir = mientorno)
        #                 cexindividuos[i]<-paste("cexi",i, sep = "")
        #                 assign(cexindividuos[i],1, envir = mientorno)
        #                 simindividuos[i]<-paste("simi",i, sep = "")
        #                 assign(simindividuos[i],"+", envir = mientorno)
        #         }#end for (i in 1:dim(x)[1])
        #         
        #         
        #         for (i in 1:dim(x)[2])
        #         {
        #                 colvariables[i]<-paste("colorv",i, sep = "")
        #                 assign(colvariables[i],"black", envir = mientorno)
        #                 textvariables[i]<-paste("labelv",i, sep = "")
        #                 assign(textvariables[i],colnames(x)[i], envir = mientorno)
        #                 cexvariables[i]<-paste("cexv",i, sep = "")
        #                 assign(cexvariables[i],1, envir = mientorno)
        #                 simvariables[i]<-paste("simv",i, sep = "")
        #                 assign(simvariables[i]," ", envir = mientorno)
        #         }#end for (i in 1:dim(x)[2])
        
        
        #############################################################################
        ### Informative window
        #############################################################################
        
        winfor<-tktoplevel()
        tkwm.title(winfor,"Classical Biplots")
        #### Frames
        
        framewi<-tkframe(winfor, relief = "flat", borderwidth = 2, background = "white")
        framewi1<-tkframe(framewi, relief = "ridge", borderwidth = 2, background = "white")
        framewi2<-tkframe(framewi, relief = "ridge", borderwidth = 2, background = "white")
        
        framewi21<-tkframe(framewi2, relief = "ridge", borderwidth = 2, background = "white")
        framewi21i<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
        framewi21a<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
        framewi21f<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
        
        framewi21ft<-tkframe(framewi21f, relief = "ridge", borderwidth = 2, background = "white")
        framewi21fb<-tkframe(framewi21f, relief = "ridge", borderwidth = 2, background = "white")
        
        framewi21fl<-tkframe(framewi21fb, relief = "ridge", borderwidth = 2, background = "white")
        framewi21fr<-tkframe(framewi21fb, relief = "ridge", borderwidth = 2, background = "white")
        
        
        framewi22<-tkframe(framewi2, relief = "flat", borderwidth = 2, background = "white")
        framewi22c<-tkframe(framewi22, relief = "flat", borderwidth = 2, background = "white")
        framewi22l<-tkframe(framewi22, relief = "flat", borderwidth = 2, background = "white")
        
        framewi221c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi222c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi223c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi224c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi225c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi226c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi227c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi228c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi229c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2210c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2211c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2212c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2213c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2214c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        framewi2215c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
        
        framewi221l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi222l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi223l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi224l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi225l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi226l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")	
        framewi227l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi228l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi229l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2210l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2211l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2212l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2213l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2214l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewi2215l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
        framewigr<-tkframe(winfor, relief = "flat", borderwidth = 2, background = "white")
        
        fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
        fontFixedWidth <- tkfont.create(family="courier",size=12)
        tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
        tkpack(tklabel(framewi1,text="BOOTSTRAP ON CLASSICAL BIPLOTS",font=fontHeading, foreground = "blue"), expand = "TRUE", side="left",expand="TRUE", fill = "both")
        tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
        
        ######Iterations   ###################################
        Niter <- tclVar(niter)
        entry.Niter <-tkentry(framewi21i,width=10,textvariable=Niter)
        tkconfigure(entry.Niter,textvariable=Niter)			
        
        tkpack(tklabel(framewi21i, text="Number of resamples"),entry.Niter, expand = "TRUE", side="left", fill = "both")
        
        ######alpha confidence intervals ###################################
        Nalpha <- tclVar(alphaic)
        entry.Nalpha <-tkentry(framewi21a,width=10,textvariable=Nalpha)
        tkconfigure(entry.Nalpha,textvariable=Nalpha)
        
        tkpack(tklabel(framewi21a, text="Confidence Level     "),entry.Nalpha, expand = "TRUE", side="left", fill = "both")
        
        tkpack(framewi21fl, framewi21fr, expand = "TRUE",side="left", fill="both")
        tkpack(framewi21ft, framewi21fb, expand = "TRUE",side="top", fill="both")
        tkpack(framewi21i, framewi21a, framewi21f, expand = "TRUE",side="top", fill="both")
        
        ######save histograms ###################################        
        tkpack(tklabel(framewi21ft,text="Save files as:"),
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        cpdf <- tkradiobutton(framewi21fl)
        cpdfb <- tkradiobutton(framewi21fl)
        rbpdfValue <- tclVar("Color pdf")
        tkconfigure(cpdf,variable=rbpdfValue,value="Color pdf")
        tkconfigure(cpdfb,variable=rbpdfValue,value="Black and white pdf")
        tkpack(tklabel(framewi21fl,text="Color pdf"),cpdf,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(framewi21fl,text="Black and white pdf"),cpdfb,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        
        
        ceps <- tkradiobutton(framewi21fr)
        cepsb <- tkradiobutton(framewi21fr)
        rbepsValue <- tclVar("Color eps")
        tkconfigure(ceps,variable=rbepsValue,value="Color eps")
        tkconfigure(cepsb,variable=rbepsValue,value="Black and white eps")
        tkpack(tklabel(framewi21fr,text="Color eps"),ceps,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(framewi21fr,text="Black and white eps"),cepsb,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        
        
        
        ###### Parameters to estimate   ###################################
        tkpack(tklabel(framewi221l, text="Calculate confidence intervals for:"), expand = "TRUE", side="left",expand="TRUE", fill = "both")
        tkpack(tklabel(framewi221c, text=" "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
        
        ##### Checkbox Godness of fit  #######
        
        cgf <- tkcheckbutton(framewi222c)
        cgfValue <- tclVar("0")
        tkconfigure(cgf,variable=cgfValue)
        
        ##### Checkbox quality of representation  #######
        
        ccr <- tkcheckbutton(framewi223c)
        ccrValue <- tclVar("0")
        tkconfigure(ccr,variable=ccrValue)
        
        ##### Checkbox quality of representation rows #######
        
        cci <- tkcheckbutton(framewi2215c)
        cciValue <- tclVar("0")
        tkconfigure(cci,variable=cciValue)
        
        ##### Checkbox eigenvalues  #######
        
        cei <- tkcheckbutton(framewi224c)
        ceiValue <- tclVar("0")
        tkconfigure(cei,variable=ceiValue)
        
        ##### Checkbox angles between variables  #######
        
        cavar <- tkcheckbutton(framewi225c)
        cavarValue <- tclVar("0")
        tkconfigure(cavar,variable=cavarValue)
        
        ##### Checkbox angles between variables  #######
        
        cavarj <- tkcheckbutton(framewi226c)
        cavarjValue <- tclVar("0")
        tkconfigure(cavarj,variable=cavarjValue)
        
        ##### Checkbox variable coordinates #######
        
        #cvcoo <- tkcheckbutton(framewi227c)
        #cvcooValue <- tclVar("0")
        #tkconfigure(cvcoo,variable=cvcooValue)
        
        ##### Checkbox CRTj #######
        
        ccrtj <- tkcheckbutton(framewi228c)
        ccrtjValue <- tclVar("0")
        tkconfigure(ccrtj,variable=ccrtjValue)
        
        ##### Checkbox CRTi #######
        
        ccrti <- tkcheckbutton(framewi2212c)
        ccrtiValue <- tclVar("0")
        tkconfigure(ccrti,variable=ccrtiValue)
        
        ##### Checkbox CREjFq #######
        
        ccrejfq <- tkcheckbutton(framewi229c)
        ccrejfqValue <- tclVar("0")
        tkconfigure(ccrejfq,variable=ccrejfqValue)
        
        ##### Checkbox CREFqEj #######
        
        ccrfqej <- tkcheckbutton(framewi2210c)
        ccrfqejValue <- tclVar("0")
        tkconfigure(ccrfqej,variable=ccrfqejValue)
        
        ##### Checkbox CREiFq #######
        
        ccreifq <- tkcheckbutton(framewi2213c)
        ccreifqValue <- tclVar("0")
        tkconfigure(ccreifq,variable=ccreifqValue)
        
        ##### Checkbox CREFqEi #######
        
        ccrfqei <- tkcheckbutton(framewi2214c)
        ccrfqeiValue <- tclVar("0")
        tkconfigure(ccrfqei,variable=ccrfqeiValue)
        
        ##### Checkbox length of variables #######
        
        clv <- tkcheckbutton(framewi2211c)
        clvValue <- tclVar("0")
        tkconfigure(clv,variable=clvValue)
        
        tkpack( tklabel(framewi222l, text="-Goodness of fit", anchor="nw"), 
                tklabel(framewi223l, text="-Quality of approximation for columns", anchor="nw"), 
                tklabel(framewi224l, text="-Eigenvalues", anchor="nw"),
                tklabel(framewi225l, text="-Angles between variables", anchor="nw"),
                tklabel(framewi226l, text="-Angles between variables and axes", anchor="nw"),
                tklabel(framewi228l, text="-Relative contribution to total variability of the column element j", anchor="nw"),
                tklabel(framewi229l, text="-Relative contribution of the column element j to the q-th factor", anchor="nw"),
                tklabel(framewi2210l, text="-Relative contribution of the q-th factor to column element j", anchor="nw"),
                tklabel(framewi2211l, text="-Variable lengths", anchor="nw"),
                tklabel(framewi2212l, text="-Relative contribution to total variability of the row element i", anchor="nw"),
                tklabel(framewi2213l, text="-Relative contribution of the row element i to the q-th factor", anchor="nw"),
                tklabel(framewi2214l, text="-Relative contribution of the q-th factor to row element i", anchor="nw"),
                tklabel(framewi2215l, text="-Quality of approximation for rows", anchor="nw"), 
                expand = "FALSE", side="top",expand="TRUE", fill = "both")
        
        tkpack(cgf, ccr, cei, cavar, cavarj, #cvcoo, 
               ccrtj, ccrejfq, ccrfqej, clv,
               ccrti, ccreifq, ccrfqei, cci,
               expand = "TRUE",side="top", fill="both")
        tkpack(framewi221l,framewi222l,framewi223l,framewi224l,framewi225l,framewi226l,framewi227l,framewi228l,framewi229l,framewi2210l,framewi2211l,framewi2212l,framewi2213l,framewi2214l,framewi2215l, expand = "TRUE",side="top", fill="both")
        tkpack(framewi221c,framewi222c,framewi223c,framewi224c,framewi225c,framewi226c,framewi227c,framewi228c,framewi229c,framewi2210c,framewi2211c,framewi2212c,framewi2213c,framewi2214c,framewi2215c, expand = "TRUE",side="top", fill="both")
        tkpack(framewi22c, framewi22l, expand = "TRUE",side="left", fill="both")
        
        
        OnOKinf <- function()
        {
                tkdestroy(winfor) 
                Xpon <<- array(data=unlist(x), dim=dim(x))
                niter <<- tclvalue(Niter)
                alphaic <<- tclvalue(Nalpha)
                cgfVal <<- as.character(tclvalue(cgfValue))
                ccrVal <<- as.character(tclvalue(ccrValue))
                cciVal <<- as.character(tclvalue(cciValue))
                ceiVal <<- as.character(tclvalue(ceiValue))
                cavarVal <<- as.character(tclvalue(cavarValue))
                cavarjVal <<- as.character(tclvalue(cavarjValue))
                #	cvcooVal <<- as.character(tclvalue(cvcooValue))
                ccrtjVal <<- as.character(tclvalue(ccrtjValue))
                ccrejfqVal <<- as.character(tclvalue(ccrejfqValue))
                ccrfqejVal <<- as.character(tclvalue(ccrfqejValue))
                ccrtiVal <<- as.character(tclvalue(ccrtiValue))
                ccreifqVal <<- as.character(tclvalue(ccreifqValue))
                ccrfqeiVal <<- as.character(tclvalue(ccrfqeiValue))
                cpdfVal <<- as.character(tclvalue(rbpdfValue))
                cepsVal <<- as.character(tclvalue(rbepsValue))
                clvVal <<- as.character(tclvalue(clvValue))
                
                ##############################################################################
                #####	Window to change labels and colors and select the biplot 
                ##############################################################################
                
                tt<-tktoplevel()
                tkwm.title(tt,"Options")
                
                #####Dropdown menu#############################
                
                topMenutt <- tkmenu(tt)
                tkconfigure(tt,menu=topMenutt)
                fileMenutt <- tkmenu(topMenutt,tearoff=FALSE)
                fileMenutrans <- tkmenu(topMenutt, tearoff=FALSE)
                
                
                tkadd(fileMenutt,"command",label="HJ-biplot",command=function() tipo<<-"RCMP")
                tkadd(fileMenutt,"command",label="GH-biplot",command=function() tipo<<-"CMP")
                tkadd(fileMenutt,"command",label="JK-biplot",command=function() tipo<<-"RMP")
                
                tkadd(topMenutt,"cascade",label="Biplot",menu=fileMenutt)
                
                tkadd(fileMenutrans,"command",label="Subtract the global mean",command=function() tChoice<<-"Subtract the global mean")
                tkadd(fileMenutrans,"command",label="Column centering",command=function() tChoice<<-"Column centering")
                tkadd(fileMenutrans,"command",label="Standardize columns",command=function() tChoice<<-"Standardize columns")
                tkadd(fileMenutrans,"command",label="Row centering",command=function() tChoice<<-"Row centering")
                tkadd(fileMenutrans,"command",label="Standardize rows",command=function() tChoice<<-"Standardize rows")
                tkadd(fileMenutrans,"command",label="Double Centering",command=function() tChoice<<-"Double centering")
                tkadd(fileMenutrans,"command",label="Raw data",command=function() tChoice<<-"Raw data")
                
                tkadd(topMenutt,"cascade",label="Transformations",menu=fileMenutrans)
                
                #### Frames
                
                framett<-tkframe(tt, relief = "flat", borderwidth = 2, background = "white")
                framett1<-tkframe(framett, relief = "ridge", borderwidth = 2, background = "white")
                framett2<-tkframe(framett, relief = "ridge", borderwidth = 2, background = "white")
                framett3<-tkframe(framett, relief = "ridge", borderwidth = 2)
                
                framet1<-tkframe(framett1, relief = "ridge", borderwidth = 2, background = "white")
                frametext1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                frameok1<-tkframe(framett1, relief = "ridge", borderwidth = 2, background = "white")
                
                framecol1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                framecol11<-tkframe(framecol1, relief = "flat", borderwidth = 2, background = "white")
                framecol12<-tkframe(framecol1, relief = "flat", borderwidth = 2, background = "white")
                
                framename1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                framename11<-tkframe(framename1, relief = "flat", borderwidth = 2, background = "white")
                framename12<-tkframe(framename1, relief = "flat", borderwidth = 2, background = "white")
                
                framecex1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                framecex11<-tkframe(framecex1, relief = "flat", borderwidth = 2, background = "white")
                framecex12<-tkframe(framecex1, relief = "flat", borderwidth = 2, background = "white")
                
                frames1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                frames11<-tkframe(frames1, relief = "flat", borderwidth = 2, background = "white")
                frames12<-tkframe(frames1, relief = "flat", borderwidth = 2, background = "white")
                
                framet2<-tkframe(framett2, relief = "ridge", borderwidth = 2, background = "white")
                frametext2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                frameok2<-tkframe(framett2, relief = "ridge", borderwidth = 2, background = "white")
                
                framecol2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                framecol21<-tkframe(framecol2, relief = "flat", borderwidth = 2, background = "white")
                framecol22<-tkframe(framecol2, relief = "flat", borderwidth = 2, background = "white")
                
                framename2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                framename21<-tkframe(framename2, relief = "flat", borderwidth = 2, background = "white")
                framename22<-tkframe(framename2, relief = "flat", borderwidth = 2, background = "white")
                
                framecex2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                framecex21<-tkframe(framecex2, relief = "flat", borderwidth = 2, background = "white")
                framecex22<-tkframe(framecex2, relief = "flat", borderwidth = 2, background = "white")
                
                framet3<-tkframe(framett3, relief = "flat", borderwidth = 2)#, background = "white")
                frametext3<-tkframe(framett3, relief = "flat", borderwidth = 2)#, background = "white")
                frameok3<-tkframe(framett3, relief = "flat", borderwidth = 2)#, background = "white")
                framett3auxgs<-tkframe(framett3, relief = "ridge", borderwidth = 2)#, background = "white")
                
                
                frames2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                #		frames21<-tkframe(frames2, relief = "flat", borderwidth = 2, background = "white")
                #		frames22<-tkframe(frames2, relief = "flat", borderwidth = 2, background = "white")
                
                framehvtitle<-tkframe(framett3auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                framehv<-tkframe(framett3auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                framehvnames<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                framehnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                framevnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                
                framehvtext<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                framehtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                framevtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                
                framegraphic<-tkframe(tt, relief = "flat", borderwidth = 2, background = "white")
                
                ##### Checkbox to show the axes or not #######
                
                cb <- tkcheckbutton(frames2)
                cbValue <- tclVar("0")
                tkconfigure(cb,variable=cbValue)
                
                
                ##############################################################################
                ##### 	List of individuals
                ##############################################################################
                
                scri <- tkscrollbar(framet1, repeatinterval=5, command=function(...)tkyview(tli,...))
                tli<-tklistbox(framet1,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scri,...),background="white")	
                tkpack(tklabel(frametext1,text="Individuals"),side="left",expand = "TRUE",fill="both")
                
                for (i in 1:(dim(Xpon)[1]))
                {
                        tkinsert(tli,"end",textindividuos[i])
                }#end for (i in 1:(dim(Xpon)[1]))
                
                tkselection.set(tli,0) #  Indexing starts at zero.
                
                OnOKi <- function()
                {
                        Choicei <- textindividuos[as.numeric(tkcurselection(tli))+1]
                        
                        ##### Color of the selected individual  #############
                        
                        indicei<<-as.numeric(tkcurselection(tli))+1
                        colori <<- colindividuos[indicei[1]]
                        tkconfigure(canvasi,bg=colori)
                        
                        ##### Text of the selected variable###############
                        
                        Namei <<- tclVar(textindividuos[indicei[1]])
                        tkconfigure(entry.Namei,textvariable=Namei)
                        
                        ##### Size of the selected variable###############
                        
                        Cexi <<- tclVar(cexindividuos[indicei[1]])
                        tkconfigure(entry.Cexi,textvariable=Cexi) 
                }#end OnOKi <- function()
                
                OK.buti <-tkbutton(frameok1,text="    OK    ",command=OnOKi)
                
                tkpack(tli,scri,expand = "TRUE", side="left", fill = "both")
                tkpack.configure(scri,side="left")
                tkpack(OK.buti,expand = "TRUE", side="left", fill = "both")
                
                #######Color#######################################
                
                indicei <<- as.numeric(tkcurselection(tli))+1
                colori<<- colindividuos[indicei[1]]
                canvasi <- tkcanvas(framecol11,width="57",height="20",bg=colori)
                
                ChangeColori <- function()
                {
                        colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colindividuos[indicei[1]],title="Choose a color"))
                        
                        if (nchar(colori)>0)
                        {
                                tkconfigure(canvasi,bg=colori)
                                colindividuos[indicei]<<-colori
                        }
                }#end ChangeColori <- function()
                
                ChangeColor.buttoni<- tkbutton(framecol12,text="Change Color",command=ChangeColori,width=4)
                tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")
                
                ######Labels   ###################################
                Namei <- textindividuos[indicei[1]]
                entry.Namei <-tkentry(framename11,width=10,textvariable=Namei)
                
                OnOKli <- function()
                {
                        NameVali <- tclvalue(Namei)
                        textindividuos[indicei[1]]<-NameVali
                        
                }#end OnOKli <- function()
                
                OK.butli <-tkbutton(framename12,text=" Change label",command=OnOKli,width=4)
                tkbind(entry.Namei, "<Return>",OnOKli)
                tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
                
                ###### Sizes   ###################################
                Cexi <- cexindividuos[indicei[1]]
                entry.Cexi <-tkentry(framecex11,width=10,textvariable=Cexi)
                
                OnOKci <- function()
                {
                        NameCexi <- tclvalue(Cexi)
                        cexindividuos[indicei]<<-NameCexi
                }#end OnOKci <- function()
                
                OK.butci <-tkbutton(framecex12,text=" Change size",command=OnOKci,width=4)
                tkbind(entry.Cexi, "<Return>",OnOKci)
                tkpack(entry.Cexi,OK.butci,expand = "TRUE", side="left", fill = "both")
                
                ######Symbols  ###################################
                
                comboBoxi <- tkwidget(frames11,"ComboBox",editable=FALSE,values=symbols, width=7)
                
                chang.symi <- function()
                {
                        simChoicei <- symbols[as.numeric(tclvalue(tcl(comboBoxi,"getvalue")))+1]
                        simindividuos[indicei]<<-simChoicei
                }#end chang.symi <- function()
                
                Change.symboli <-tkbutton(frames12,text="   Change symbol   ",command=chang.symi,width=4, height=1)
                tkpack(comboBoxi,Change.symboli,side="left",expand="TRUE", fill="both")
                
                ##### List of variables ###########################
                
                scrv <- tkscrollbar(framet2, repeatinterval=5, command=function(...)tkyview(tlv,...))
                tlv<-tklistbox(framet2,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scrv,...),background="white")
                tkpack(tklabel(frametext2,text="Variables"),side="left",expand = "TRUE",fill="both")
                
                for (i in 1:dim(Xpon)[2])
                {
                        tkinsert(tlv,"end",textvariables[i])
                }#end for (i in 1:dim(Xpon)[2])
                
                tkselection.set(tlv,0) #  Indexing starts at zero.
                
                OnOKv <- function()
                {
                        Choicev <- textvariables[as.numeric(tkcurselection(tlv))+1]
                        
                        ##### Color of the selected variable  #############
                        indicev<<-as.numeric(tkcurselection(tlv))+1
                        colorv <<- colvariables[indicev[1]]
                        tkconfigure(canvasv,bg=colorv)
                        
                        ##### Text of the selected variable  ###############
                        
                        Namev <<- tclVar(textvariables[indicev[1]])
                        tkconfigure(entry.Namev,textvariable=Namev)
                        
                        ##### Size of the selected variable  ###############
                        
                        Cexv <<- tclVar(cexvariables[indicev[1]])
                        tkconfigure(entry.Cexv,textvariable=Cexv)
                }#end OnOKv <- function()
                
                OK.butv <-tkbutton(frameok2,text="    OK    ",command=OnOKv)
                tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")
                tkpack(tlv,scrv,expand = "TRUE", side="left", fill = "both")
                tkpack.configure(scrv,side="left")
                tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")
                tkfocus(tt)
                
                #######Color#######################################
                indicev<-as.numeric(tkcurselection(tlv))+1
                colorv <- colvariables[indicev[1]]
                canvasv <- tkcanvas(framecol21,width="57",height="20",bg=colorv)
                
                ChangeColorv <- function()
                {
                        colorv <<- tclvalue(tcl("tk_chooseColor",initialcolor=colvariables[indicev[1]],title="Choose a color"))
                        if (nchar(colorv)>0)
                        {
                                tkconfigure(canvasv,bg=colorv)
                                colvariables[indicev]<<-colorv
                        }#end if (nchar(colorv)>0)
                }#end ChangeColorv <- function()
                
                ChangeColor.buttonv<- tkbutton(framecol22,text="Change Color",command=ChangeColorv,width=4)
                tkpack(canvasv,ChangeColor.buttonv,expand = "TRUE", side="left", fill = "both")
                
                ######Labels  ###################################
                
                Namev <- textvariables[indicev[1]]
                entry.Namev <-tkentry(framename21,width=10, textvariable=Namev)
                
                OnOKlv <- function()
                {
                        NameValv <- tclvalue(Namev)
                        textvariables[indicev[1]] <<-NameValv
                        
                }#end OnOKlv <- function()
                
                OK.butlv <-tkbutton(framename22,text=" Change label",command=OnOKlv,width=4)
                tkbind(entry.Namev, "<Return>",OnOKlv)
                tkpack(entry.Namev,OK.butlv,expand = "TRUE", side="left", fill = "both")
                
                ###### Sizes  ###################################
                
                Cexv <- cexvariables[indicev[1]]
                entry.Cexv <-tkentry(framecex21,width=10, textvariable=Cexv)
                
                OnOKcv <- function()
                {
                        NameCexv <<- tclvalue(Cexv)
                        cexvariables[indicev] <<-NameCexv
                }#end OnOKcv <- function()
                
                OK.butcv <-tkbutton(framecex22,text=" Change size",command=OnOKcv,width=4)
                tkbind(entry.Cexv, "<Return>",OnOKcv)
                tkpack(entry.Cexv,OK.butcv,expand = "TRUE", side="left", fill = "both")
                
                Graphics <- function()
                {
                        
                        barvp<-tktoplevel()
                        tkdestroy(tt)
                        tkwm.title(barvp,"Eigenvalues")
                        hescale <<- tclvalue(entryvalueh)
                        vescale <<- tclvalue(entryvaluev)
                        
                        plotbar<-function()
                        {
                                Xbar<-transforma(tChoice,Xpon)
                                descom <<- La.svd(Xbar)
                                sumaRvalprop<-sum((descom$d)^2)
                                inerciatot<<-(descom$d[1:length(descom$d)])^2/sumaRvalprop
                                barplot(descom$d, col="blue", xlab="", ylab="", names.arg=round(inerciatot, digits=2))
                                msginertia<<-"Proportion of inertia explained by each axis:"
                                
                                for (i in 1:length(descom$d))
                                {
                                        msginertia<<-paste(msginertia, "\n",i, "\t", round(inerciatot[i]*100, digits=2), "%")
                                }#end for (i in 1:length(descom$d))
                        }#end plotbar<-function()
                        
                        imgbar <<- tkrplot(barvp,fun=plotbar,hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                        tk2tip(imgbar, msginertia)
                        
                        Onaxis <- function()
                        {
                                nejes <- tclvalue(numaxis)
                                nejes <-as.numeric(nejes)
                                
                                if (nejes > length(descom$d))
                                {
                                        msg <- paste("The maximum number of dimensions is ",length(descom$d))
                                        tkmessageBox(message=msg)
                                }else{
                                        
                                        tkdestroy(barvp)
                                        nejes <- as.integer(nejes)
                                        
                                        
                                        resultados<- ex.biplot(Xpon, tipo, nejes, tChoice, tindi=textindividuos, tvar=textvariables)
                                        ejes<<-resultados$ejes
                                        descom<<-resultados$descom
                                        coindividuos<<-resultados$coindividuos
                                        covariables<<-resultados$covariables
                                        suma2valprop<<-resultados$suma2valprop
                                        inercia<<-resultados$inercia
                                        cuminer<<-resultados$cuminer
                                        bonajuste<<-resultados$bonajuste
                                        calcol<<-resultados$calcol
                                        calfilas<<-resultados$calfilas
                                        covartotal<<-resultados$covartotal
                                        coindividuosnam<<-resultados$coindividuosnam
                                        covariablesnam<<- resultados$covariablesnam
                                        CRTi<<-resultados$CRTi
                                        CRTj<<-resultados$CRTj
                                        longitudprin<<-resultados$longitudprin
                                        CREiFq<<-resultados$CREiFq
                                        CREjFq<<-resultados$CREjFq
                                        CRFqEi<<-resultados$CRFqEi
                                        CRFqEj<<-resultados$CRFqEj
                                        ang<<-resultados$ang
                                        angax<<-resultados$angax
                                        cat("File saved in:    ",file="Results.txt")
                                        cat(getwd(),file="temp.txt")        				
                                        file.append("Results.txt","temp.txt")	
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")		
                                        cat("CONTRIBUTIONS:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution to total variability of the row element i:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRTi, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution to total variability of the column element j:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRTj, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the row element i to the q-th factor:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CREiFq, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the column element j to the q-th factor:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CREjFq, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the q-th factor to row element i:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRFqEi, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the q-th factor to column element j:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRFqEj, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        
                                        if (tipo != "RCMP"){
                                                cat("\n",file="temp.txt")					
                                                file.append("Results.txt","temp.txt")	
                                                cat("Goodness of fit:  ",file="temp.txt")					
                                                file.append("Results.txt","temp.txt")					
                                                cat(round(bonajuste, digits=2),file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                cat(" %",file="temp.txt")					
                                                file.append("Results.txt","temp.txt")	
                                        }#end if (tipo != "RCMP")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Quality of approximation for rows:  ",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        cat(round(calfilas, digits=2),file="temp.txt")
                                        file.append("Results.txt","temp.txt")
                                        cat(" %",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Quality of approximation for columns:  ",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        cat(round(calcol, digits=2),file="temp.txt")
                                        file.append("Results.txt","temp.txt")
                                        cat(" %",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Individual coordinates:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(coindividuosnam, digits=2),file="temp.txt",sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Variable coordinates:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(covariablesnam, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Variable lengths:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(longitudprin, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Angles between variables:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(ang, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")		
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Angles between variables and axes:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(angax, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Eigenvalues: \n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(descom$d, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")
                                        file.append("Results.txt","temp.txt")
                                        cat("Proportion of inertia explained by each axis: \n",file="temp.txt")
                                        file.append("Results.txt","temp.txt")
                                        write.table(round(inercia, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        file.show("Results.txt")
                                        file.remove("temp.txt")
                                        
                                        
                                        datos<-rbind(coindividuos,covariables)
                                        textos<-datos
                                        limix <- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                        limiy <- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                        
                                        centro<-c(0,0)
                                        
                                        xCoords<-textos[,dim1]
                                        yCoords<-textos[,dim2]
                                        labelsVec <- c()
                                        sizesVec <- c()
                                        colores<-c()
                                        simbolos<-c()
                                        
                                        colores<-c(colindividuos, colvariables)
                                        labelsVec<-c(textindividuos, textvariables)
                                        sizesVec<-c(cexindividuos, cexvariables)
                                        simbolos<-c(simindividuos, simvariables)
                                        
                                        #                                         for (k in 1:length(colindividuos))
                                        #                                         {
                                        #                                                 colores<-c(colores, get(colindividuos[k], envir = mientorno))
                                        #                                                 labelsVec<-c(labelsVec, get(textindividuos[k], envir = mientorno))
                                        #                                                 sizesVec<-c(sizesVec, get(cexindividuos[k], envir = mientorno))
                                        #                                                 simbolos<-c(simbolos, get(simindividuos[k], envir = mientorno))
                                        #                                         }#end for (k in 1:length(colindividuos))
                                        #                                         
                                        #                                         for (k in 1:length(colvariables))
                                        #                                         {
                                        #                                                 colores<-c(colores, get(colvariables[k], envir = mientorno))
                                        #                                                 labelsVec<-c(labelsVec, get(textvariables[k], envir = mientorno))
                                        #                                                 sizesVec<-c(sizesVec, get(cexvariables[k], envir = mientorno))
                                        #                                                 simbolos<-c(simbolos, get(simvariables[k], envir = mientorno))
                                        #                                         }#end for (k in 1:length(colvariables))
                                        #                                         
                                        indexLabeled<-c(1:length(xCoords))
                                        indexLabeledaux<-c()
                                        labeledPoints <- list()
                                        ################Show axes or not
                                        cbVal <<- as.character(tclvalue(cbValue))
                                        
                                        wgr <- tktoplevel()
                                        tkwm.title(wgr,"Graph")
                                        
                                        
                                        normalLine <- function(variableele, colorele, coordenadas, colorcoor) 
                                        { 
                                                A <- variableele
                                                B <- centro
                                                
                                                slopeAB <- (B[[2]] - A[[2]])/(B[[1]] - A[[1]]) 
                                                slopeNorm <- -1/slopeAB 
                                                a <- A[[2]] - slopeAB * A[[1]]
                                                
                                                b<-c()
                                                xintersect<-c()
                                                yintersect<-c()
                                                
                                                for(i in 1:dim(coordenadas)[1])
                                                {
                                                        b[i] <- coordenadas[i,2] - slopeNorm * coordenadas[i,1] 
                                                        
                                                        xintersect[i] <- (b[i] - a)/(slopeAB - slopeNorm) 
                                                        yintersect[i] <- b[i] + slopeNorm * xintersect[i] 
                                                        
                                                }#end for (i in 1:dim(coordenadas)[1])
                                                
                                                abline(a =a,b=slopeAB, col=colorele, lwd=3)
                                                
                                                for(i in 1:dim(coordenadas)[1])
                                                        segments(xintersect[i], yintersect[i], coordenadas[i,1], coordenadas[i,2],lty=2, col=colorcoor[i]) 
                                        } #end normalLine
                                        
                                        plotFunctiond <- function(screen=TRUE)
                                        {
                                                tclvalue(Limix1) <- limix[1]
                                                tclvalue(Limix2) <- limix[2]
                                                tclvalue(Limiy1) <- limiy[1]
                                                tclvalue(Limiy2) <- limiy[2]
                                                
                                                labelsVec <- c()
                                                sizesVec <- c()
                                                colores<-c()
                                                simbolos<-c()
                                                
                                                colores<-c(colindividuos, colvariables)
                                                labelsVec<-c(textindividuos, textvariables)
                                                sizesVec<-c(cexindividuos, cexvariables)
                                                simbolos<-c(simindividuos, simvariables)
                                                
                                                xCoords<<-textos[,dim1]
                                                yCoords<<-textos[,dim2]
                                                params <- par(bg="white")
                                                plot(datos[,c(dim1,dim2)],main= tit_graph,type="n",xlab=paste("Axis", dim1, ":", round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste("Axis", dim2, ":", round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limix * 1.05, ylim = limiy * 1.05)
                                                
                                                if(proj=="normal")
                                                {
                                                        
                                                        if(clb=="normal")
                                                        {
                                                                points(coindividuos[,dim1],coindividuos[,dim2],pch=simbolos[1:length(simindividuos)],col=colores[1:length(colindividuos)])    
                                                                if (length(indexLabeled)>0)
                                                                        for (i in (1:length(indexLabeled)))
                                                                        {
                                                                                indexClosest <- indexLabeled[i]
                                                                                text(xCoords[indexClosest],yCoords[indexClosest], labels=labelsVec[indexClosest], col= colores[indexClosest],cex=as.numeric(sizesVec)[indexClosest])
                                                                        }#end for (i in (1:length(indexLabeled)))
                                                                
                                                        }else{
                                                                
                                                                for(i in 1:as.numeric(nclust))
                                                                {
                                                                        clusteri <-coindividuos[which(pertb==i),c(dim1,dim2), drop=FALSE]
                                                                        clusteri <- t(t(clusteri))
                                                                        hpts <- chull(clusteri)
                                                                        hpts <- c(hpts, hpts[1])
                                                                        if(length(hpts)>2)
                                                                        {
                                                                                lines(clusteri[hpts,], col=unique(pertb)[i]+1)                                                                                
                                                                        }
                                                                }
                                                                points(coindividuos[,dim1],coindividuos[,dim2],pch=simbolos[1:length(simindividuos)],col=pertb+1)    
                                                                colorescl<-c(pertb+1, colvariables)
                                                                if (length(indexLabeled)>0)
                                                                        for (i in (1:length(indexLabeled)))
                                                                        {
                                                                                indexClosest <- indexLabeled[i]
                                                                                text(xCoords[indexClosest],yCoords[indexClosest], labels=labelsVec[indexClosest], col= colorescl[indexClosest],cex=as.numeric(sizesVec)[indexClosest])
                                                                        }#end for (i in (1:length(indexLabeled)))
                                                                
                                                                if(clb=="km")
                                                                {
                                                                        points(hc$centers[unique(pertcentr),c(dim1,dim2)], pch=8, col=unique(pertb)+1)
                                                                }
                                                                
                                                                
                                                        }
                                                        
                                                        
                                                        arrows(centro[1],centro[2],covariables[,dim1],covariables[,dim2],col=colores[1+length(colindividuos):length(colores)],#lty="dotted",
                                                               length=0.10, lwd=0.08)
                                                        points(centro[1],centro[2],pch=18,col="black")
                                                        
                                                        if (cbVal=="1"){
                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                        }#end if (cbVal=="1")
                                                        
                                                        
                                                        
                                                }else{
                                                        colvariablesp <- rep("grey",dim(Xpon)[2])
                                                        colvariablesp[Choiceproj] <-colvariables[Choiceproj]
                                                        coloresp<-c(colindividuos, colvariablesp)
                                                        points(coindividuos[,dim1],coindividuos[,dim2],pch=simbolos[1:length(simindividuos)],col=colores[1:length(colindividuos)])
                                                        
                                                        if (length(indexLabeled)>0)
                                                                for (i in (1:length(indexLabeled)))
                                                                {
                                                                        indexClosest <- indexLabeled[i]
                                                                        text(xCoords[indexClosest],yCoords[indexClosest], labels=labelsVec[indexClosest], col= coloresp[indexClosest],cex=as.numeric(sizesVec)[indexClosest])
                                                                }#end for (i in (1:length(indexLabeled)))
                                                        
                                                        normalLine(covariables[Choiceproj,c(dim1,dim2)], colvariablesp[Choiceproj], coindividuos[,c(dim1,dim2)], colindividuos)
                                                        
                                                }
                                                parPlotSize <<- par("plt")
                                                usrCoords   <<- par("usr")
                                                par(params)
                                        }#end plotFunctiond <- function(screen=TRUE)
                                        
                                        g3d<-function()
                                        {
                                                if (nejes>2)
                                                { 
                                                        zCoords<-datos[,dim3]
                                                        bg3d("white")
                                                        aspect3d("iso")
                                                        lims <- par3d("bbox")
                                                        
                                                        if (cbVal=="1"){
                                                                axes3d()
                                                        }#end if (cbVal=="1")
                                                        
                                                        labelsVec <- c()
                                                        sizesVec <- c()
                                                        colores<-c()
                                                        simbolos<-c()
                                                        
                                                        colores<-c(colindividuos, colvariables)
                                                        labelsVec<-c(textindividuos, textvariables)
                                                        sizesVec<-c(cexindividuos, cexvariables)
                                                        simbolos<-c(simindividuos, simvariables)
                                                        
                                                        points3d(xCoords,yCoords,zCoords, color=colores)
                                                        texts3d(xCoords, yCoords, zCoords,labelsVec,color=colores, cex= as.numeric(sizesVec))
                                                        
                                                        for (i in 1:(dim(covariables)[1]))
                                                        {
                                                                linea<-rbind(covariables[i,c(dim1, dim2, dim3)],c(0,0,0))	
                                                                segments3d(linea[,1],linea[,2], linea[,3],color=colores[i+length(colindividuos)])
                                                        }#end for (i in 1:(dim(covariables)[1]))
                                                        
                                                        rgl.bringtotop()
                                                }else{
                                                        msg <- "You have selected less than 3 dimensions. 3D-graph not available"
                                                        tkmessageBox(message=msg)
                                                }#end if (nejes>2)
                                        }#end g3d<-function()
                                        
                                        #############################################################################
                                        ### Functions to save the graph
                                        #############################################################################
                                        SaveFileJPG <- function() {
                                                FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Jpeg files} {.jpg .jpeg}} {{All files} *}"))
                                                if (nchar(FileName)) {
                                                        nn <- nchar(FileName)
                                                        if (nn < 5 || substr(FileName, nn - 3, nn) != ".jpg") 
                                                                FileName <- paste(FileName, ".jpg", sep = "")
                                                        jpeg(FileName, width = 8, height = 8, units = "in", restoreConsole = FALSE, res = 96, quality = 50)
                                                        plotFunctiond(screen = FALSE)
                                                        dev.off()
                                                }#end if (nchar(FileName))
                                        }#end SaveFileJPG <- function()
                                        
                                        SaveFilePDF <- function() {
                                                FileName <- tclvalue(tkgetSaveFile(filetypes = "{{PDF files} {.pdf}} {{All files} *}"))
                                                if (nchar(FileName)) {
                                                        nn <- nchar(FileName)
                                                        if (nn < 5 || substr(FileName, nn - 3, nn) != ".pdf") 
                                                                FileName <- paste(FileName, ".pdf", sep = "")
                                                        pdf(FileName, width = 7, height = 7)
                                                        plotFunctiond(screen = FALSE)
                                                        dev.off()
                                                }#end if (nchar(FileName))
                                        }#end SaveFilePDF <- function()
                                        
                                        SaveFileBmp <- function() {
                                                FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Bitmap files} {.bmp}} {{All files} *}"))
                                                if (nchar(FileName)) {
                                                        nn <- nchar(FileName)
                                                        if (nn < 5 || substr(FileName, nn - 3, nn) != ".bmp") 
                                                                FileName <- paste(FileName, ".bmp", sep = "")
                                                        bmp(FileName, width = 8, height = 8, units = "in", restoreConsole = FALSE, res = 96)
                                                        plotFunctiond(screen = FALSE)
                                                        dev.off()
                                                }#end if (nchar(FileName))
                                        }#end SaveFileBmp <- function()
                                        
                                        SaveFilePng <- function() {
                                                FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Png files} {.png}} {{All files} *}"))
                                                if (nchar(FileName)) {
                                                        nn <- nchar(FileName)
                                                        if (nn < 5 || substr(FileName, nn - 3, nn) != ".png") 
                                                                FileName <- paste(FileName, ".png", sep = "")
                                                        png(FileName, width = 8, height = 8, units = "in", restoreConsole = FALSE, res = 96)
                                                        plotFunctiond(screen = FALSE)
                                                        dev.off()
                                                }#end if (nchar(FileName))
                                        }#end SaveFilePng <- function()
                                        
                                        SaveFileeps <- function() {
                                                FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Eps files} {.eps}} {{All files} *}"))
                                                if (nchar(FileName)) {
                                                        nn <- nchar(FileName)
                                                        if (nn < 5 || substr(FileName, nn - 3, nn) != ".eps") 
                                                                FileName <- paste(FileName, ".eps", sep = "")
                                                        postscript(FileName, width = 8, height = 8)
                                                        plotFunctiond(screen = FALSE)
                                                        dev.off()
                                                }#end if (nchar(FileName))
                                        }#end SaveFilePng <- function()
                                        
                                        projections <- function(project)
                                        {        				
                                                proj<<-project
                                                if(proj=="normal")
                                                {
                                                        tkrreplot(img)
                                                        return()
                                                }#end if(proj=="normal")
                                                if(proj=="v") 
                                                {
                                                        wproj <- tktoplevel()
                                                        tkwm.title(wproj, "Select variable")
                                                        
                                                        scrproj <- tkscrollbar(wproj, repeatinterval=5, command=function(...)tkyview(tlproj,...))
                                                        tlproj<-tklistbox(wproj,height=6,width=42,yscrollcommand=function(...)tkset(scrproj,...),background="white")
                                                        
                                                        for (i in 1:dim(Xpon)[2])
                                                        {
                                                                tkinsert(tlproj,"end",textvariables[i])
                                                        }#end for (i in 1:Numbercuant)
                                                        
                                                        tkselection.set(tlproj,0) #  Indexing starts at zero.
                                                        
                                                        OnOKproj <- function()
                                                        {
                                                                Choiceproj <<- as.numeric(tkcurselection(tlproj))+1
                                                                tkdestroy(wproj)
                                                                tkrreplot(img)
                                                        }#end OnOKproj <- function()
                                                        
                                                        OK.butp <- tkbutton(wproj, text = "   OK   ", command = OnOKproj)
                                                        tkpack(tlproj,scrproj,expand = "TRUE", side="left", fill = "both")
                                                        tkpack.configure(scrproj,side="left")
                                                        tkpack(OK.butp,expand = "TRUE", side="left", fill = "both")
                                                        tkfocus(wproj)
                                                        tkwait.window(wproj)
                                                }#end if       
                                        }#end projections <-function()
                                        
                                        clusterbip <- function(cl)
                                        {
                                                clb<<-cl
                                                if(clb=="h")
                                                {
                                                        wclust <- tktoplevel()
                                                        tkwm.title(wclust, "Hierarchical clustering")
                                                        
                                                        framehier<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                        framehier1<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                        framehier2<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                        framehier3<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                        
                                                        linkm <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
                                                        distancias <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
                                                        
                                                        comboBoxdist <- tkwidget(framehier2,"ComboBox",editable=FALSE,values=distancias,width=10, text= distancias[1])
                                                        comboBoxlink <- tkwidget(framehier2,"ComboBox",editable=FALSE,values=linkm,width=10, text= linkm[1])
                                                        
                                                        
                                                        OnOKclust <- function()
                                                        {
                                                                distchoos<<-distancias[as.numeric(tclvalue(tcl(comboBoxdist,"getvalue")))+1]
                                                                linkchoos<<-linkm[as.numeric(tclvalue(tcl(comboBoxlink,"getvalue")))+1]
                                                                nclust <<- tclvalue(n_clust)
                                                                hc<<-hclust(dist(coindividuosnam, method=distchoos), method=linkchoos)
                                                                pertb<<-cutree(hc, k=nclust)
                                                                pert<<-pertb
                                                                for(i in 1:nclust)
                                                                {
                                                                        pert[which(pertb==i)]<<-unique(pertb[hc$order])[i]
                                                                }
                                                                dcolor<<-colour_clusters(hc, nclust, col=unique(pert)+1)
                                                                plot(dcolor)
                                                                tkdestroy(wclust)
                                                                tkrreplot(img)
                                                        }#end OnOKclust <- function()
                                                        
                                                        
                                                        n_clust<-tclVar(nclust)
                                                        entry.clust <-tkentry(framehier2, width="50",textvariable=n_clust, bg="white")
                                                        tkbind(entry.clust, "<Return>",OnOKclust)
                                                        
                                                        OK.butc <- tkbutton(framehier3, text = "   OK   ", command = OnOKclust)
                                                        
                                                        tkpack(tklabel(framehier1,text="Number of cluster:    "),
                                                               tklabel(framehier1,text="Distance:     "),
                                                               tklabel(framehier1,text="Link method:   "), expand = "TRUE", side="top", fill = "both")
                                                        
                                                        tkpack(OK.butc, expand = "TRUE", side="left", fill = "both")
                                                        tkpack(entry.clust, comboBoxdist, comboBoxlink, expand = "TRUE", side="top", fill = "both")
                                                        
                                                        #tkpack(framehier21, framehier22, framehier23, side="top", expand="TRUE", fill="both")
                                                        tkpack(framehier1, framehier2, side="left", expand="TRUE", fill="both")
                                                        tkpack(framehier, framehier3, side="top", expand="TRUE", fill="both")
                                                        
                                                        tkfocus(wclust)
                                                        tkwait.window(wclust)
                                                        
                                                }else{
                                                        if(clb=="km")
                                                        {
                                                                wclust <- tktoplevel()
                                                                tkwm.title(wclust, "K-means")
                                                                
                                                                framehier<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                framehier1<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                framehier2<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                framehier3<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                
                                                                algoritmo <- c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
                                                                
                                                                comboBoxalg <- tkwidget(framehier2,"ComboBox",editable=FALSE,values=algoritmo,width=10, text= algoritmo[1])
                                                                
                                                                
                                                                OnOKclust <- function()
                                                                {
                                                                        algorchoos<<-algoritmo[as.numeric(tclvalue(tcl(comboBoxalg,"getvalue")))+1]
                                                                        nclust <<- tclvalue(n_clust)
                                                                        niterclust <<- tclvalue(n_iter)
                                                                        nstart <<- tclvalue(n_start)
                                                                        
                                                                        
                                                                        hc<<-kmeans(coindividuosnam, centers=nclust, iter.max=niterclust, nstart=nstart)
                                                                        pertb<<-hc$cluster
                                                                        pertcentr<<-pertb
                                                                        pert<<-pertb
                                                                        for(i in 1:length(unique(pertb)))
                                                                        {
                                                                                pert[which(pertb==unique(pertb)[i])]<<-i
                                                                        }
                                                                        pertb<<-pert
                                                                        tkdestroy(wclust)
                                                                        tkrreplot(img)
                                                                }#end OnOKclust <- function()
                                                                
                                                                
                                                                n_clust<-tclVar(nclust)
                                                                entry.clust <-tkentry(framehier2, width="50",textvariable=n_clust, bg="white")
                                                                tkbind(entry.clust, "<Return>",OnOKclust)
                                                                
                                                                n_iter<-tclVar(niterclust)
                                                                entry.iterclus <-tkentry(framehier2, width="50",textvariable=n_iter, bg="white")
                                                                tkbind(entry.iterclus, "<Return>",OnOKclust)
                                                                
                                                                n_start<-tclVar(nstart)
                                                                entry.start <-tkentry(framehier2, width="50",textvariable=n_start, bg="white")
                                                                tkbind(entry.start, "<Return>",OnOKclust)
                                                                
                                                                OK.butc <- tkbutton(framehier3, text = "   OK   ", command = OnOKclust)
                                                                
                                                                tkpack(tklabel(framehier1,text="Algorithm:    "),
                                                                       tklabel(framehier1,text="Number of cluster:    "),
                                                                       tklabel(framehier1,text="Iterations:     "),
                                                                       tklabel(framehier1,text="Random sets:   "), expand = "TRUE", side="top", fill = "both")
                                                                
                                                                tkpack(OK.butc, expand = "TRUE", side="left", fill = "both")
                                                                tkpack(comboBoxalg, entry.clust, entry.iterclus, entry.start, expand = "TRUE", side="top", fill = "both")
                                                                
                                                                #tkpack(framehier21, framehier22, framehier23, side="top", expand="TRUE", fill="both")
                                                                tkpack(framehier1, framehier2, side="left", expand="TRUE", fill="both")
                                                                tkpack(framehier, framehier3, side="top", expand="TRUE", fill="both")
                                                                
                                                                tkfocus(wclust)
                                                                tkwait.window(wclust)  
                                                        }else{
                                                                if(clb=="kmed")
                                                                {
                                                                        wclust <- tktoplevel()
                                                                        tkwm.title(wclust, "K-medoids")
                                                                        
                                                                        framehier<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                        framehier1<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                        framehier2<-tkframe(framehier, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                        framehier3<-tkframe(wclust, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                                        
                                                                        distancias <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
                                                                        
                                                                        comboBoxdist <- tkwidget(framehier2,"ComboBox",editable=FALSE,values=distancias,width=10, text= distancias[1])
                                                                        
                                                                        
                                                                        OnOKclust <- function()
                                                                        {
                                                                                distchoos<<-distancias[as.numeric(tclvalue(tcl(comboBoxdist,"getvalue")))+1]
                                                                                nclust <<- tclvalue(n_clust)
                                                                                
                                                                                pamx<<-pam(dist(coindividuosnam, method=distchoos), k=nclust, diss=TRUE)
                                                                                pertb<<-pamx$clustering
                                                                                tkdestroy(wclust)
                                                                                tkrreplot(img)
                                                                        }#end OnOKclust <- function()
                                                                        
                                                                        
                                                                        n_clust<-tclVar(nclust)
                                                                        entry.clust <-tkentry(framehier2, width="50",textvariable=n_clust, bg="white")
                                                                        tkbind(entry.clust, "<Return>",OnOKclust)
                                                                        
                                                                        
                                                                        OK.butc <- tkbutton(framehier3, text = "   OK   ", command = OnOKclust)
                                                                        
                                                                        tkpack(tklabel(framehier1,text="Distance:    "),
                                                                               tklabel(framehier1,text="Number of cluster:    "), expand = "TRUE", side="top", fill = "both")
                                                                        
                                                                        tkpack(OK.butc, expand = "TRUE", side="left", fill = "both")
                                                                        tkpack(comboBoxdist, entry.clust, expand = "TRUE", side="top", fill = "both")
                                                                        
                                                                        #tkpack(framehier21, framehier22, framehier23, side="top", expand="TRUE", fill="both")
                                                                        tkpack(framehier1, framehier2, side="left", expand="TRUE", fill="both")
                                                                        tkpack(framehier, framehier3, side="top", expand="TRUE", fill="both")
                                                                        
                                                                        tkfocus(wclust)
                                                                        tkwait.window(wclust)  
                                                                        
                                                                }else{
                                                                        tkrreplot(img)
                                                                        return()
                                                                }
                                                        }
                                                }
                                        }#end clusterbip
                                        
                                        showaxes <- function()
                                        {
                                                if(cbVal=="1")
                                                {
                                                        cbVal<<-"0"
                                                        tkrreplot(img)
                                                }else{
                                                        cbVal<<-"1"
                                                        tkrreplot(img)
                                                }#end if(cbVal=="1")        
                                        }#end showaxes
                                        
                                        changetit <- function()
                                        {
                                                ctwin<-tktoplevel()
                                                tkwm.title(ctwin,"Change title")
                                                OnOKchantit <- function()
                                                {
                                                        tit_graph <<- tclvalue(tit_gr)
                                                        tkrreplot(img)
                                                        tkdestroy(ctwin)
                                                        
                                                }
                                                OK.butchantit<-tkbutton(ctwin,text=" Change ", command=OnOKchantit,  bg= "lightblue", width=20, foreground = "navyblue")
                                                tkbind(OK.butchantit, "<Return>",OnOKchantit)
                                                
                                                tit_gr<-tclVar(tit_graph)
                                                entry.tit <-tkentry(ctwin, width="50",textvariable=tit_gr, bg="white")
                                                tkbind(entry.tit, "<Return>",OnOKchantit)
                                                
                                                
                                                tkpack(tklabel(ctwin,text="New title:    "),entry.tit, expand = "TRUE", side="left", fill = "both")
                                                tkpack(OK.butchantit)
                                                
                                                tkfocus(ctwin)
                                                
                                        }#end changetit
                                        
                                        
                                        
                                        topMenugr <- tkmenu(wgr)
                                        tkconfigure(wgr, menu = topMenugr)
                                        menuFile <- tkmenu(topMenugr, tearoff = FALSE)
                                        menuSaveAs <- tkmenu(topMenugr, tearoff = FALSE)
                                        menu3d <- tkmenu(topMenugr, tearoff = FALSE) 
                                        menuproj <- tkmenu(topMenugr, tearoff = FALSE)                  			
                                        menuopt <-tkmenu(topMenugr, tearoff = FALSE)
                                        menuclust <- tkmenu(topMenugr, tearoff = FALSE)
                                        
                                        tkadd(menuFile, "command", label = "Copy image", command = function() {tkrreplot(img)})
                                        tkadd(menuFile, "cascade", label = "Save image", menu = menuSaveAs)
                                        tkadd(menuSaveAs, "command", label = "PDF file", command = function() {SaveFilePDF()})
                                        tkadd(menuSaveAs, "command", label = "Eps file", command = function() {SaveFileeps()})
                                        #          tkadd(menuSaveAs, "command", label = "Bmp file", command = function() {SaveFileBmp()})
                                        tkadd(menuSaveAs, "command", label = "Png file", command = function() {SaveFilePng()})
                                        tkadd(menuSaveAs, "command", label = "Jpg/Jpeg file", command = function() {SaveFileJPG()})
                                        tkadd(menuFile, "separator")
                                        tkadd(menuFile, "command", label = "Exit", command = function() {tkdestroy(wgr)})
                                        tkadd(menu3d, "command", label = "3D", command = function() {g3d()})
                                        tkadd(topMenugr, "cascade", label = "File", menu = menuFile)
                                        tkadd(menuFile, "separator")
                                        tkadd(topMenugr, "cascade", label = "3D", menu = menu3d)
                                        tkadd(menuopt, "command", label = "Change title", command = function() {changetit()})
                                        tkadd(menuopt, "command", label = "Show/Hide axes", command = function() {showaxes()})
                                        tkadd(menuproj, "command", label = "Variables", command = function() {projections(project="v")})
                                        tkadd(menuproj, "command", label = "Back to original graph", command = function() {projections(project="normal")})						
                                        tkadd(topMenugr, "cascade", label = "Projections", menu = menuproj)
                                        tkadd(topMenugr, "cascade", label = "Options", menu = menuopt)
                                        tkadd(menuclust, "command", label = "Hierarchical cluster with biplot coordinates", command = function() {clusterbip(cl="h")})
                                        tkadd(menuclust, "command", label = "K-means with biplot coordinates", command = function() {clusterbip(cl="km")})
                                        tkadd(menuclust, "command", label = "K-medoids with biplot coordinates", command = function() {clusterbip(cl="kmed")})
                                        tkadd(menuclust, "command", label = "Back to original graph", command = function() {clusterbip(cl="normal")})
                                        tkadd(topMenugr, "cascade", label = "Cluster", menu = menuclust)
                                        
                                        img <<- tkrplot(wgr,fun=plotFunctiond, hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                                        
                                        framedim1<-tkframe(wgr, relief = "flat", borderwidth = 2, background = "whitesmoke")
                                        framecomb<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                        framescal<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                        framescala<-tkframe(framescal, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                        framescalb<-tkframe(framescal, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                        framerefr<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                        
                                        comboBoxdim1 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=10, text= dim1)
                                        comboBoxdim2 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=10, text= dim2)
                                        if (nejes>2)
                                                comboBoxdim3 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=10, text= dim3)
                                        
                                        
                                        chang.symdim1 <- function()
                                        {
                                                dim1 <<-as.numeric(tclvalue(tcl(comboBoxdim1,"getvalue")))+1
                                                dim2 <<-as.numeric(tclvalue(tcl(comboBoxdim2,"getvalue")))+1
                                                if (nejes>2)
                                                        dim3 <<-as.numeric(tclvalue(tcl(comboBoxdim3,"getvalue")))+1
                                                
                                                if((dim1==dim1ant)&(dim2==dim2ant))
                                                {
                                                        limix[1] <<- as.numeric(tclvalue(Limix1))
                                                        limix[2] <<- as.numeric(tclvalue(Limix2))
                                                        limiy[1] <<- as.numeric(tclvalue(Limiy1))
                                                        limiy[2] <<- as.numeric(tclvalue(Limiy2))
                                                }else{
                                                        dim1ant<<-dim1
                                                        dim2ant<<-dim2
                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                        
                                                }
                                                tkrreplot(img)
                                        }#end chang.symdim1 <- function()
                                        
                                        Change.symboldim1 <-tkbutton(framerefr,text="Refresh",command=chang.symdim1, bg= "lightblue", width=15, foreground = "navyblue")
                                        if(nejes>2){
                                                tkpack(tklabel(framecomb, text="Select X, Y and Z axes numbers:"), expand="FALSE", side= "top", fill ="both")
                                                tkpack(comboBoxdim1, comboBoxdim2, comboBoxdim3, side="top", fill="x")
                                                
                                                #}else{
                                                #tkpack(tklabel(framedim1, text="Select X and Y axes numbers:"), expand="FALSE", side= "left", fill ="both")
                                                #tkpack(comboBoxdim1, comboBoxdim2, Change.symboldim1, side="left", expand="FALSE")
                                        }
                                        Limix1 <- tclVar(limix[1])
                                        entry.limix1 <-tkentry(framescalb,width=10,textvariable=Limix1, bg="white")
                                        tkconfigure(entry.limix1,textvariable=Limix1)                  
                                        tkbind(entry.limix1, "<Return>",chang.symdim1)
                                        
                                        Limix2 <- tclVar(limix[2])
                                        entry.limix2 <-tkentry(framescalb,width=10,textvariable=Limix2, bg="white")
                                        tkconfigure(entry.limix2,textvariable=Limix2)                  
                                        tkbind(entry.limix2, "<Return>",chang.symdim1)
                                        
                                        Limiy1 <- tclVar(limiy[1])
                                        entry.limiy1 <-tkentry(framescalb,width=10,textvariable=Limiy1, bg="white")
                                        tkconfigure(entry.limiy1,textvariable=Limiy1)                  
                                        tkbind(entry.limiy1, "<Return>",chang.symdim1)
                                        
                                        Limiy2 <- tclVar(limiy[2])
                                        entry.limiy2 <-tkentry(framescalb,width=10,textvariable=Limiy2, bg="white")
                                        tkconfigure(entry.limiy2,textvariable=Limiy2) 
                                        tkbind(entry.limiy2, "<Return>",chang.symdim1)
                                        tkpack(tklabel(framescala, text="-X"),
                                               tklabel(framescala, text="+X"),
                                               tklabel(framescala, text="-Y"),
                                               tklabel(framescala, text="+Y"),
                                               expand = "TRUE",side="top", fill="both")
                                        tkpack(entry.limix1,entry.limix2,entry.limiy1,entry.limiy2,side="top", expand="TRUE", fill="both")
                                        tkpack(tklabel(framescal, text="Zoom:"), expand = "TRUE",side="top", fill="both")
                                        tkpack(framescala, framescalb, side="left", expand ="TRUE", fill="both")
                                        tkpack(Change.symboldim1, side="top", expand="TRUE", fill="both")
                                        tkpack(framecomb, framescal, framerefr, side="top", expand="TRUE", fill="both")
                                        tkpack(framedim1, side="left", expand="FALSE", fill="x")
                                        tkpack(img, side="top", expand="TRUE", fill="both")
                                        
                                        
                                        labelClosestPointd <- function(xClick,yClick,imgXcoords,imgYcoords)
                                        {
                                                labelsVec <- c()
                                                sizesVec <- c()
                                                colores<-c()
                                                simbolos<-c()
                                                
                                                colores<-c(colindividuos, colvariables)
                                                labelsVec<-c(textindividuos, textvariables)
                                                sizesVec<-c(cexindividuos, cexvariables)
                                                simbolos<-c(simindividuos, simvariables)
                                                
                                                squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                indexClosest <<- which.min(squared.Distance)
                                                mm<-tktoplevel() 	
                                                tkwm.title(mm, labelsVec[indexClosest])	
                                                
                                                framemm1<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm2<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm3<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm4<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                
                                                
                                                
                                                
                                                
                                                colorder <- colores[indexClosest]
                                                canvasder <- tkcanvas(framemm1,width="120",height="20",bg=colorder)
                                                
                                                ChangeColorder <- function()
                                                {
                                                        
                                                        colorder <- tclvalue(tcl("tk_chooseColor",initialcolor=colores[indexClosest],title="Choose a color"))
                                                        
                                                        if (nchar(colorder)>0)
                                                        {
                                                                tkconfigure(canvasder,bg=colorder)
                                                                colores[indexClosest]<-colorder
                                                                
                                                                colindividuos<<-colores[1:length(colindividuos)]
                                                                colvariables<<-colores[(length(colindividuos)+1):(length(colindividuos)+length(colvariables))]
                                                                
                                                        }#end if (nchar(colorder)>0)
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end ChangeColorder <- function()
                                                
                                                ChangeColor.buttonder<- tkbutton(framemm1,text="Change Color",command=ChangeColorder)
                                                tkpack(canvasder,ChangeColor.buttonder,expand = "TRUE", side="left", fill = "both")
                                                
                                                Nameder <- labelsVec[indexClosest]
                                                tclvalue(Nameder) <- labelsVec[indexClosest]
                                                entry.Nameder <-tkentry(framemm2,width="10",textvariable=Nameder)
                                                NameValder <- Nameder 
                                                
                                                OnOKlder <- function()
                                                {
                                                        NameValder <- tclvalue(Nameder)
                                                        labelsVec[indexClosest]<-NameValder
                                                        
                                                        textindividuos<<-labelsVec[1:length(textindividuos)]
                                                        textvariables<<-labelsVec[(length(textindividuos)+1):(length(textindividuos)+length(textvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end OnOKlder <- function()
                                                
                                                OK.butlder <-tkbutton(framemm2,text=" Change label",command=OnOKlder,width=2)
                                                tkbind(entry.Nameder, "<Return>",OnOKlder)
                                                tkpack(entry.Nameder,OK.butlder,expand = "TRUE", side="left", fill = "both")
                                                
                                                Cexder <- sizesVec[indexClosest]  
                                                tclvalue(Cexder) <- sizesVec[indexClosest]
                                                entry.Cexder <-tkentry(framemm3,width="10",textvariable=Cexder)
                                                NameCexder <- Cexder 
                                                
                                                OnOKcder <- function()
                                                {
                                                        NameCexder <- tclvalue(Cexder)
                                                        sizesVec[indexClosest]<-NameCexder
                                                        
                                                        cexindividuos<<-sizesVec[1:length(cexindividuos)]
                                                        cexvariables<<-sizesVec[(length(cexindividuos)+1):(length(cexindividuos)+length(cexvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end OnOKcder <- function()
                                                
                                                OK.butcder <-tkbutton(framemm3,text=" Change size",command=OnOKcder,width=2)
                                                tkbind(entry.Cexder, "<Return>",OnOKcder)
                                                tkpack(entry.Cexder,OK.butcder, expand = "TRUE", side="left", fill = "both")
                                                
                                                comboBox <- tkwidget(framemm4,"ComboBox",editable=FALSE,values=symbols,width=10, text= simbolos[indexClosest])
                                                
                                                chang.symder <- function()
                                                {
                                                        simChoice <-symbols[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
                                                        simbolos[indexClosest]<-simChoice
                                                        
                                                        simindividuos<<-simbolos[1:length(simindividuos)]
                                                        simvariables<<-simbolos[(length(simindividuos)+1):(length(simindividuos)+length(simvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end chang.symder <- function()
                                                
                                                if(indexClosest %in% c(length(simindividuos)+1):(length(simindividuos)+length(simvariables)))
                                                {}else{
                                                        Change.symbolder <-tkbutton(framemm4,text="   Change symbol   ",command=chang.symder,width=6)
                                                        tkpack(comboBox,Change.symbolder,side="left",expand="TRUE", fill="both")
                                                        
                                                }
                                                
                                                if(proj=="normal")
                                                {
                                                        tkpack(framemm1,framemm2,framemm3,framemm4, expand = "TRUE", side="top", fill = "both")
                                                }else{
                                                        tkpack(framemm1,framemm2, expand = "TRUE", side="top", fill = "both")
                                                }
                                        }#end labelClosestPointd <- function(xClick,yClick,imgXcoords,imgYcoords)
                                        
                                        OnLeftClick.up <- function(x,y)
                                        {
                                                msg <- ("-To change the label press Yes.\n-To remove it press No.\n-If you do not want to do anything press Cancel.")
                                                mbval<- tkmessageBox(title="Change of label", message=msg,type="yesnocancel",icon="question")
                                                if (tclvalue(mbval)=="yes"){  
                                                        indexLabeled <<- c(indexLabeled,indexClosest)
                                                }#end if (tclvalue(mbval)=="yes")
                                                
                                                if(tclvalue(mbval)=="no"){
                                                        indexLabeledaux<<-c()
                                                        for (i in (1:length(indexLabeled)))
                                                        {
                                                                if (indexLabeled[i]!=indexClosest)
                                                                        indexLabeledaux <<- c(indexLabeledaux,indexLabeled[i])
                                                        }#end for (i in (1:length(indexLabeled)))
                                                        indexLabeled<<-indexLabeledaux 
                                                }#end if(tclvalue(mbval)=="no")
                                                
                                                if(tclvalue(mbval)=="cancel"){
                                                        textos[indexClosest,dim1] <<- anteriorx
                                                        textos[indexClosest,dim2] <<- anteriory
                                                }#end if(tclvalue(mbval)=="cancel")
                                                tkrreplot(img)
                                        }#end OnLeftClick.up <- function(x,y)
                                        
                                        OnLeftClick.move <- function(x,y)
                                        {
                                                xClick <- x
                                                yClick <- y
                                                width  = as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height = as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin = parPlotSize[1] * width
                                                xMax = parPlotSize[2] * width
                                                yMin = parPlotSize[3] * height
                                                yMax = parPlotSize[4] * height
                                                
                                                rangeX = usrCoords[2] - usrCoords[1]
                                                rangeY = usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords = (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords = (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord = usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord = usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                
                                                textos[indexClosest,dim1]<<-xPlotCoord
                                                textos[indexClosest,dim2]<<-yPlotCoord
                                                
                                                tkrreplot(img) 
                                        }#end OnLeftClick.move <- function(x,y)
                                        
                                        OnLeftClick.down <- function(x,y)
                                        {
                                                anteriorx <- NULL
                                                anteriory <- NULL
                                                xClick <- x
                                                yClick <- y
                                                width  = as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height = as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin = parPlotSize[1] * width
                                                xMax = parPlotSize[2] * width
                                                yMin = parPlotSize[3] * height
                                                yMax = parPlotSize[4] * height
                                                
                                                rangeX = usrCoords[2] - usrCoords[1]
                                                rangeY = usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords = (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords = (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord = usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord = usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                indexClosest <<- which.min(squared.Distance) 
                                                
                                                anteriorx <<- textos[indexClosest,dim1]
                                                anteriory <<- textos[indexClosest,dim2]	  
                                        }#end OnLeftClick.down <- function(x,y)
                                        
                                        OnRightClick <- function(x,y)
                                        {
                                                xClick <- x
                                                yClick <- y
                                                width  = as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height = as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin = parPlotSize[1] * width
                                                xMax = parPlotSize[2] * width
                                                yMin = parPlotSize[3] * height
                                                yMax = parPlotSize[4] * height
                                                
                                                rangeX = usrCoords[2] - usrCoords[1]
                                                rangeY = usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords = (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords = (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord = usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord = usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                labelClosestPointd(xClick,yClick,imgXcoords,imgYcoords)
                                        }#end OnRightClick <- function(x,y)
                                        
                                        tkbind(img, "<B1-Motion>",OnLeftClick.move)
                                        tkbind(img, "<ButtonPress-1>",OnLeftClick.down)
                                        tkbind(img, "<ButtonRelease-1>",OnLeftClick.up)
                                        tkconfigure(img,cursor="pencil")
                                        
                                        tkbind(img, "<Button-3>",OnRightClick)
                                        tkconfigure(img,cursor="pencil")
                                }#end if (nejes > length(descom$d))
                                
                                
                                ###############################################################################				
                                ####	Bootstrap
                                ###############################################################################
                                varrep<-rep(list(Xpon),niter)
                                varresampletot<-lapply(varrep, resample_boot)
                                varresample<-lapply(varresampletot, function(x)x[[1]])
                                indicesample<-lapply(varresampletot, function(x)x[[2]])
                                
                                bootresult<-mapply(ex.biplot, varresample, tipo=tipo, nejes=nejes, tChoice=tChoice, tindi=lapply(varresample, rownames),
                                                   tvar=lapply(varresample, colnames))
                                bondadesalm<-bootresult[8,]
                                calidadesalm<-bootresult[9,]
                                calfilasalm<-bootresult[10,]
                                autovaloresalm<-bootresult[23,]
                                angulosalm<-bootresult[21,]
                                anguloejesalm<-bootresult[22,]
                                crtjalm<-bootresult[15,]
                                longalm<-bootresult[16,]
                                crejfqalm<-bootresult[18,]
                                crfqejalm<-bootresult[20,]
                                crtiaux<-bootresult[14,]
                                crtialm<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesample), unlist(crtiaux))
                                creifqaux<-bootresult[17,]
                                creifqalm <- vector("list",nejes)
                                for (i in 1:nejes)
                                {
                                        creifqalm[[i]]<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesample), unlist(lapply(creifqaux, function(x, ejes) x[,ejes], i)))
                                }
                                crfqeiaux<-bootresult[19,]
                                crfqeialm <- vector("list",nejes)
                                for (i in 1:nejes)
                                {
                                        crfqeialm[[i]]<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesample), unlist(lapply(crfqeiaux, function(x, ejes) x[,ejes], i)))
                                }
                                
                                
                                jackresample<-lapply(1:(dim(Xpon)[1]), function(i)Xpon[-i,]) 
                                indicesjack<-lapply(1:dim(Xpon)[1], function(i) (1:dim(Xpon)[1])[-i])
                                jackresult<-mapply(ex.biplot, jackresample, tipo=tipo, nejes=nejes, tChoice=tChoice,tindi=lapply(jackresample, rownames),
                                                   tvar=lapply(jackresample, colnames))
                                bondadesjackr<-jackresult[8,]
                                calcoljackr<-jackresult[9,]
                                calfilasjackr<-jackresult[10,]
                                descomjackr<-jackresult[23,]
                                angulosjackr<-jackresult[21,]
                                anguloejesjackr<-jackresult[22,]
                                crtjjackr<-jackresult[15,]
                                longjackr<-jackresult[16,]
                                crejfqjackr<-jackresult[18,]
                                crfqejjackr<-jackresult[20,]
                                crtiauxj<-jackresult[14,]
                                crtijackr<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesjack), unlist(crtiauxj))
                                
                                creifqauxj<-jackresult[17,]
                                creifqjackr <- vector("list",nejes)
                                for (i in 1:nejes)
                                {
                                        creifqjackr[[i]]<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesjack), unlist(lapply(creifqauxj, function(x, ejes) x[,ejes], i)))
                                }
                                crfqeiauxj<-jackresult[19,]
                                crfqeijackr <- vector("list",nejes)
                                for (i in 1:nejes)
                                {
                                        crfqeijackr[[i]]<-lapply(1:length(textindividuos), function(x, indices, datab) datab[which(indices==x)], unlist(indicesjack), unlist(lapply(crfqeiauxj, function(x, ejes) x[,ejes], i)))
                                }
                                
                                ####crear vectores con coordenadas de variables
                                coorvarrot<-bootresult[13,]
                                coorvarrot<-array(c(unlist(covariables),unlist(coorvarrot)), dim=c(dim(Xpon)[2], nejes, as.numeric(niter)+1))
                                
                                out.var<-procGPA(coorvarrot, reflect=TRUE, distances=FALSE, pcaoutput=FALSE)
                                
                                plot(out.var$rotated[,dim1,], out.var$rotated[,dim2,], type="n", main="Bootstrap Coordinates (variables)", xlab=paste("Dimension", dim1), ylab=paste("Dimension", dim2), asp=1/1)
                                
                                for(i in 1:dim(covariablesnam)[1])
                                {
                                        points(out.var$rotated[i,dim1,],out.var$rotated[i,dim2,], col=colorescoor[i], pch=20)
                                        
                                }
                                
                                abline(v=0)
                                abline(h=0)
                                
                                
                                for(i in 1:dim(covariablesnam)[1])
                                {
                                        hpts <- chull(t(out.var$rotated[i,c(dim1,dim2),]))
                                        hpts <- c(hpts, hpts[1])
                                        lines(t(out.var$rotated[i,c(dim1,dim2),hpts]), col=colorescoor[i])
                                }
                                text(out.var$rotated[,dim1,1], out.var$rotated[,dim2,1], labels=textvariables)
                                
                                
                                
                                #########################################################################
                                ### Guardar resultados
                                #########################################################################
                                
                                cat("File saved in:    ",file="Resultsbootstrap.txt")
                                cat(getwd(),file="temp.txt")					
                                file.append("Resultsbootstrap.txt","temp.txt")	
                                cat("\n",file="temp.txt")					
                                file.append("Resultsbootstrap.txt","temp.txt")
                                cat("\n",file="temp.txt")					
                                file.append("Resultsbootstrap.txt","temp.txt")		
                                
                                alphaic <<- tclvalue(Nalpha)
                                alphaic <<- as.numeric(alphaic)
                                liminf <- (1 - alphaic*0.01) / 2
                                limsup <- 1 - (1 - alphaic*0.01) / 2
                                
                                nombresvariables <-textvariables
                                
                                titulo <-c("Obs. Value","Mean","SE","Bias","IC t-boot lo","IC t-boot up","IC perc lo","IC perc up","IC BCa lo","IC BCa up")
                                
                                
                                ### Goodness of fit
                                
                                if (cgfVal=="1")
                                {
                                        if (tipo != "RCMP"){		
                                                
                                                calc.cgf <-c()
                                                cgf.mean <-c()
                                                se.cgf <-c()
                                                sesgo.cgf <-c()
                                                ic.t.cgfinf <-c()
                                                ic.t.cgfsup <-c()
                                                ic.p.cgfinf <-c()
                                                ic.p.cgfsup <-c()
                                                ic.bca.cgfinf <-c()
                                                ic.bca.cgfsup <-c()
                                                
                                                calc.cgf <-cal.ic(sapply(bondadesalm,function(x) x[1]), liminf, limsup, bonajuste, sapply(bondadesjackr,function(x) x[1]), niter)
                                                cgf.mean <- c(cgf.mean, calc.cgf[1])
                                                se.cgf <- c(se.cgf,calc.cgf[2])
                                                sesgo.cgf <- c(sesgo.cgf,calc.cgf[3])
                                                ic.t.cgfinf <- c(ic.t.cgfinf,calc.cgf[4])
                                                ic.t.cgfsup <- c(ic.t.cgfsup,calc.cgf[5])
                                                ic.p.cgfinf <- c(ic.p.cgfinf,calc.cgf[6])
                                                ic.p.cgfsup <- c(ic.p.cgfsup,calc.cgf[7])
                                                ic.bca.cgfinf <- c(ic.bca.cgfinf,calc.cgf[8])
                                                ic.bca.cgfsup <- c(ic.bca.cgfsup,calc.cgf[9])
                                                
                                                pdf(paste("Histogram of Goodness of fit", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(bondadesalm,function(x) x[1]), main="Histogram", xlab="Goodness of fit")
                                                
                                                if(cpdfVal=="Color pdf")
                                                {
                                                        abline(v=cgf.mean, lwd=2, col="blue")
                                                        abline(v=bonajuste, lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=cgf.mean, lwd=2)
                                                        abline(v=bonajuste, lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(bondadesalm,function(x) x[1]))
                                                dev.off()
                                                
                                                
                                                postscript(paste("Histogram of Goodness of fit", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(bondadesalm,function(x) x[1]), main="Histogram", xlab="Goodness of fit")
                                                
                                                if(cepsVal=="Color eps")
                                                {
                                                        abline(v=cgf.mean, lwd=2, col="blue")
                                                        abline(v=bonajuste, lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=cgf.mean, lwd=2)
                                                        abline(v=bonajuste, lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(bondadesalm,function(x) x[1]))
                                                dev.off()                                                                               
                                                
                                                
                                                calc.cgf <-array(cbind(bonajuste, cgf.mean, se.cgf, sesgo.cgf, ic.t.cgfinf, ic.t.cgfsup, ic.p.cgfinf, ic.p.cgfsup, ic.bca.cgfinf, ic.bca.cgfsup),
                                                                 dim=c(1,10))
                                                calc.cgf <- as.data.frame(calc.cgf)
                                                colnames(calc.cgf) <- titulo
                                                rownames(calc.cgf) <-c("Goodnes of Fit")
                                                
                                                
                                                cat("\n",file="temp.txt")					
                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                cat("Goodness of fit:  \n",file="temp.txt")					
                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                write.table(round(calc.cgf, digits=3),file="temp.txt", sep="\t", dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                        }#end if (tipo != "RCMP")
                                }#end if (cgfVal=="1")
                                
                                ### Quality of representation     
                                
                                if (ccrVal=="1")
                                {
                                        calc.ccr <-c()
                                        ccr.mean <-c()
                                        se.ccr <-c()
                                        sesgo.ccr <-c()
                                        ic.t.ccrinf <-c()
                                        ic.t.ccrsup <-c()
                                        ic.p.ccrinf <-c()
                                        ic.p.ccrsup <-c()
                                        ic.bca.ccrinf <-c()
                                        ic.bca.ccrsup <-c()
                                        
                                        calc.ccr <-cal.ic(sapply(calidadesalm,function(x) x[1]), liminf, limsup, calcol, sapply(calcoljackr,function(x) x[1]), niter)
                                        ccr.mean <- c(ccr.mean, calc.ccr[1])
                                        se.ccr <- c(se.ccr,calc.ccr[2])
                                        sesgo.ccr <- c(sesgo.ccr,calc.ccr[3])
                                        ic.t.ccrinf <- c(ic.t.ccrinf,calc.ccr[4])
                                        ic.t.ccrsup <- c(ic.t.ccrsup,calc.ccr[5])
                                        ic.p.ccrinf <- c(ic.p.ccrinf,calc.ccr[6])
                                        ic.p.ccrsup <- c(ic.p.ccrsup,calc.ccr[7])
                                        ic.bca.ccrinf <- c(ic.bca.ccrinf,calc.ccr[8])
                                        ic.bca.ccrsup <- c(ic.bca.ccrsup,calc.ccr[9])
                                        
                                        pdf(paste("Histogram of Quality of Approximation for columns", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                        par(mfrow=c(1,2))
                                        hist(sapply(calidadesalm,function(x) x[1]), main="Histogram", xlab="Quality of Approximation for columns")
                                        
                                        if(cpdfVal=="Color pdf")
                                        {
                                                abline(v=ccr.mean, lwd=2, col="blue")
                                                abline(v=calcol, lty =2, lwd=2, col="red")
                                        }else{
                                                abline(v=ccr.mean, lwd=2)
                                                abline(v=calcol, lty =2, lwd=2)
                                        }        
                                        qqnorm(sapply(calidadesalm,function(x) x[1]))
                                        dev.off()
                                        
                                        
                                        postscript(paste("Histogram of Quality of Approximation for columns", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                        par(mfrow=c(1,2))
                                        hist(sapply(calidadesalm,function(x) x[1]), main="Histogram", xlab="Quality of Approximation for columns")
                                        
                                        if(cepsVal=="Color eps")
                                        {
                                                abline(v=ccr.mean, lwd=2, col="blue")
                                                abline(v=calcol, lty =2, lwd=2, col="red")
                                        }else{
                                                abline(v=ccr.mean, lwd=2)
                                                abline(v=calcol, lty =2, lwd=2)
                                        }        
                                        qqnorm(sapply(calidadesalm,function(x) x[1]))
                                        dev.off()                                                                               
                                        
                                        
                                        calc.ccr <-array(cbind(calcol, ccr.mean, se.ccr, sesgo.ccr, ic.t.ccrinf, ic.t.ccrsup, ic.p.ccrinf, ic.p.ccrsup, ic.bca.ccrinf, ic.bca.ccrsup),
                                                         dim=c(1,10))
                                        calc.ccr <- as.data.frame(calc.ccr)
                                        colnames(calc.ccr) <- titulo
                                        rownames(calc.ccr) <-c("Quality of Approximation for columns")
                                        
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Quality of approximation for columns:  \n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.ccr, digits=3),file="temp.txt", sep="\t", dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (ccrVal=="1")
                                
                                ### Quality of representation rows    
                                
                                if (cciVal=="1")
                                {
                                        calc.cci <-c()
                                        cci.mean <-c()
                                        se.cci <-c()
                                        sesgo.cci <-c()
                                        ic.t.cciinf <-c()
                                        ic.t.ccisup <-c()
                                        ic.p.cciinf <-c()
                                        ic.p.ccisup <-c()
                                        ic.bca.cciinf <-c()
                                        ic.bca.ccisup <-c()
                                        
                                        calc.cci <-cal.ic(unlist(calfilasalm), liminf, limsup, calfilas, unlist(calfilasjackr), niter)
                                        cci.mean <- c(cci.mean, calc.cci[1])
                                        se.cci <- c(se.cci,calc.cci[2])
                                        sesgo.cci <- c(sesgo.cci,calc.cci[3])
                                        ic.t.cciinf <- c(ic.t.cciinf,calc.cci[4])
                                        ic.t.ccisup <- c(ic.t.ccisup,calc.cci[5])
                                        ic.p.cciinf <- c(ic.p.cciinf,calc.cci[6])
                                        ic.p.ccisup <- c(ic.p.ccisup,calc.cci[7])
                                        ic.bca.cciinf <- c(ic.bca.cciinf,calc.cci[8])
                                        ic.bca.ccisup <- c(ic.bca.ccisup,calc.cci[9])
                                        
                                        pdf(paste("Histogram of Quality of Approximation for rows", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                        par(mfrow=c(1,2))
                                        hist(unlist(calfilasalm), main="Histogram", xlab="Quality of Approximation for rows")
                                        
                                        if(cpdfVal=="Color pdf")
                                        {
                                                abline(v=cci.mean, lwd=2, col="blue")
                                                abline(v=calfilas, lty =2, lwd=2, col="red")
                                        }else{
                                                abline(v=cci.mean, lwd=2)
                                                abline(v=calfilas, lty =2, lwd=2)
                                        }        
                                        qqnorm(unlist(calfilasalm))
                                        dev.off()
                                        
                                        
                                        postscript(paste("Histogram of Quality of Approximation for rows", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                        par(mfrow=c(1,2))
                                        hist(unlist(calfilasalm), main="Histogram", xlab="Quality of Approximation for rows")
                                        
                                        if(cepsVal=="Color eps")
                                        {
                                                abline(v=cci.mean, lwd=2, col="blue")
                                                abline(v=calfilas, lty =2, lwd=2, col="red")
                                        }else{
                                                abline(v=cci.mean, lwd=2)
                                                abline(v=calfilas, lty =2, lwd=2)
                                        }        
                                        qqnorm(unlist(calfilasalm))
                                        dev.off()                                                                               
                                        
                                        
                                        calc.cci <-array(cbind(calfilas, cci.mean, se.cci, sesgo.cci, ic.t.cciinf, ic.t.ccisup, ic.p.cciinf, ic.p.ccisup, ic.bca.cciinf, ic.bca.ccisup),
                                                         dim=c(1,10))
                                        calc.cci <- as.data.frame(calc.cci)
                                        colnames(calc.cci) <- titulo
                                        rownames(calc.cci) <-c("Quality of Approximation for rows")
                                        
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Quality of approximation for rows:  \n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.cci, digits=3),file="temp.txt", sep="\t", dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (cciVal=="1")
                                
                                ### Eigenvalues
                                
                                if (ceiVal=="1")
                                {
                                        
                                        calc.cei <-c()
                                        cei.mean <-c()
                                        se.cei <-c()
                                        sesgo.cei <-c()
                                        ic.t.ceiinf <-c()
                                        ic.t.ceisup <-c()
                                        ic.p.ceiinf <-c()
                                        ic.p.ceisup <-c()
                                        ic.bca.ceiinf <-c()
                                        ic.bca.ceisup <-c()
                                        
                                        
                                        for (i in 1:length(autovaloresalm[[1]]))
                                        { 
                                                calc.cei <-cal.ic(sapply(autovaloresalm, function(x) x[i]), liminf, limsup, descom$d[i], sapply(descomjackr, function(x) x[i]), niter)
                                                cei.mean <- c(cei.mean, calc.cei[1])
                                                se.cei <- c(se.cei,calc.cei[2])
                                                sesgo.cei <- c(sesgo.cei,calc.cei[3])
                                                ic.t.ceiinf <- c(ic.t.ceiinf,calc.cei[4])
                                                ic.t.ceisup <- c(ic.t.ceisup,calc.cei[5])
                                                ic.p.ceiinf <- c(ic.p.ceiinf,calc.cei[6])
                                                ic.p.ceisup <- c(ic.p.ceisup,calc.cei[7])
                                                ic.bca.ceiinf <- c(ic.bca.ceiinf,calc.cei[8])
                                                ic.bca.ceisup <- c(ic.bca.ceisup,calc.cei[9])
                                                
                                                pdf(paste("Histogram of eigenvalue", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(autovaloresalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                
                                                if(cpdfVal=="Color pdf")
                                                {
                                                        abline(v=cei.mean[i], lwd=2, col="blue")
                                                        abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=cei.mean[i], lwd=2)
                                                        abline(v=descom$d[i], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(autovaloresalm, function(x) x[i]))
                                                dev.off()
                                                
                                                
                                                postscript(paste("Histogram of eigenvalue", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(autovaloresalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                
                                                if(cepsVal=="Color eps")
                                                {
                                                        abline(v=cei.mean[i], lwd=2, col="blue")
                                                        abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=cei.mean[i], lwd=2)
                                                        abline(v=descom$d[i], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(autovaloresalm, function(x) x[i]))
                                                dev.off()                                                                               
                                        }#end for (i in 1:length(autovaloresalm))
                                        
                                        
                                        calc.cei <-array(cbind(descom$d, cei.mean, se.cei, sesgo.cei, ic.t.ceiinf, ic.t.ceisup, ic.p.ceiinf, ic.p.ceisup, ic.bca.ceiinf, ic.bca.ceisup),
                                                         dim=c(length(descom$d),10))
                                        calc.cei <- as.data.frame(calc.cei)
                                        colnames(calc.cei) <- titulo
                                        
                                        nombreseig<-c()
                                        for (i in 1: length(descom$d))
                                        {
                                                nombreseig <-c(nombreseig, paste("Eigenvalue",i, sep=""))
                                        }#end for (i in 1: length(descom$d))
                                        
                                        rownames(calc.cei) <- nombreseig
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Eigenvalues: \n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.cei, digits=3), file="temp.txt", sep="\t", dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (ceiVal=="1")
                                
                                ### Angles between variables
                                
                                if (cavarVal=="1")
                                {
                                        cat("\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Angles between variables:",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")			
                                        
                                        for (i in 1:(dim(ang)[1]-1))
                                        { 
                                                calc.cavar <-c()
                                                cavar.mean <-c()
                                                se.cavar <-c()
                                                sesgo.cavar <-c()
                                                ic.t.cavarinf <-c()
                                                ic.t.cavarsup <-c()
                                                ic.p.cavarinf <-c()
                                                ic.p.cavarsup <-c()
                                                ic.bca.cavarinf <-c()
                                                ic.bca.cavarsup <-c()
                                                
                                                for (j in  (i+1) : dim(ang)[2])
                                                {
                                                        calc.cavar <-cal.ic(sapply(angulosalm, function(x) x[i,j]), liminf, limsup, ang[i,j], sapply(angulosjackr, function(x) x[i,j]), niter)
                                                        cavar.mean <- c(cavar.mean, calc.cavar[1])
                                                        se.cavar <- c(se.cavar,calc.cavar[2])
                                                        sesgo.cavar <- c(sesgo.cavar,calc.cavar[3])
                                                        ic.t.cavarinf <- c(ic.t.cavarinf,calc.cavar[4])
                                                        ic.t.cavarsup <- c(ic.t.cavarsup,calc.cavar[5])
                                                        ic.p.cavarinf <- c(ic.p.cavarinf,calc.cavar[6])
                                                        ic.p.cavarsup <- c(ic.p.cavarsup,calc.cavar[7])
                                                        ic.bca.cavarinf <- c(ic.bca.cavarinf,calc.cavar[8])
                                                        ic.bca.cavarsup <- c(ic.bca.cavarsup,calc.cavar[9])
                                                        
                                                        
                                                        pdf(paste("Histogram of angles between ", nombresvariables[i]," and ", nombresvariables[j], ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(angulosalm, function(x) x[i,j]), main="Histogram", xlab=paste("Angle between ", nombresvariables[i],"\n and ", nombresvariables[j]))
                                                        
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=cavar.mean[j-i], lwd=2, col="blue")
                                                                abline(v=ang[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=cavar.mean[j-i], lwd=2)
                                                                abline(v=ang[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(angulosalm, function(x) x[i,j]))
                                                        dev.off()
                                                        
                                                        
                                                        postscript(paste("Histogram of angles between ", nombresvariables[i]," and ", nombresvariables[j], ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(angulosalm, function(x) x[i,j]), main="Histogram", xlab=paste("Angle between ", nombresvariables[i],"\n and ", nombresvariables[j]))
                                                        
                                                        if(cepsVal=="Color eps")
                                                        {
                                                                abline(v=cavar.mean[j-i], lwd=2, col="blue")
                                                                abline(v=ang[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=cavar.mean[j-i], lwd=2)
                                                                abline(v=ang[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(angulosalm, function(x) x[i,j]))
                                                        dev.off()  
                                                        
                                                }#for (j in  (i+1) : dim(angulosalm)[2])
                                                
                                                calc.cavar <-array(cbind(ang[-(1:i),i], cavar.mean, se.cavar, sesgo.cavar, ic.t.cavarinf, ic.t.cavarsup, ic.p.cavarinf, ic.p.cavarsup, ic.bca.cavarinf, ic.bca.cavarsup),
                                                                   dim=c(dim(ang)[1]-i,10))
                                                calc.cavar <- as.data.frame(calc.cavar)
                                                colnames(calc.cavar) <- titulo
                                                rownames(calc.cavar) <- nombresvariables[-(1:i)]
                                                
                                                cat("\n Angles between ", nombresvariables[i]," and ", "\n",file="temp.txt")					
                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                write.table(round(calc.cavar, digits=3), file="temp.txt", sep="\t", dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")     
                                        }#end for (i in 1:(dim(angulosalm)[1]-1))
                                }#end if (cavarVal=="1")
                                
                                ### Angles between variables and axes
                                
                                if (cavarjVal=="1")
                                {
                                        cat("\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Angles between variables and axes:",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        
                                        for (i in 1:dim(covartotal)[1])
                                        {	 
                                                calc.cavarj <-c()
                                                cavarj.mean <-c()
                                                se.cavarj <-c()
                                                sesgo.cavarj <-c()
                                                ic.t.cavarjinf <-c()
                                                ic.t.cavarjsup <-c()
                                                ic.p.cavarjinf <-c()
                                                ic.p.cavarjsup <-c()
                                                ic.bca.cavarjinf <-c()
                                                ic.bca.cavarjsup <-c()
                                                
                                                for (j in  1 : 2)
                                                {
                                                        calc.cavarj <-cal.ic(sapply(anguloejesalm, function(x) x[i,j]), liminf, limsup, angax[i,j], sapply(anguloejesjackr, function(x) x[i,j]), niter)
                                                        cavarj.mean <- c(cavarj.mean, calc.cavarj[1])
                                                        se.cavarj <- c(se.cavarj,calc.cavarj[2])
                                                        sesgo.cavarj <- c(sesgo.cavarj,calc.cavarj[3])
                                                        ic.t.cavarjinf <- c(ic.t.cavarjinf,calc.cavarj[4])
                                                        ic.t.cavarjsup <- c(ic.t.cavarjsup,calc.cavarj[5])
                                                        ic.p.cavarjinf <- c(ic.p.cavarjinf,calc.cavarj[6])
                                                        ic.p.cavarjsup <- c(ic.p.cavarjsup,calc.cavarj[7])
                                                        ic.bca.cavarjinf <- c(ic.bca.cavarjinf,calc.cavarj[8])
                                                        ic.bca.cavarjsup <- c(ic.bca.cavarjsup,calc.cavarj[9])
                                                        
                                                        
                                                        pdf(paste("Histogram of angles between ", nombresvariables[i]," and axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(anguloejesalm, function(x) x[i,j]), main="Histogram", xlab=paste("Angle between ", nombresvariables[i],"\n and axis ", j))
                                                        
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=cavarj.mean[j], lwd=2, col="blue")
                                                                abline(v=angax[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=cavarj.mean[j], lwd=2)
                                                                abline(v=angax[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(anguloejesalm, function(x) x[i,j]))
                                                        dev.off()
                                                        
                                                        
                                                        postscript(paste("Histogram of angles between ", nombresvariables[i]," and axis ",j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(anguloejesalm, function(x) x[i,j]), main="Histogram", xlab=paste("Angle between ", nombresvariables[i],"\n and axis ", j))
                                                        
                                                        if(cepsVal=="Color eps")
                                                        {
                                                                abline(v=cavarj.mean[j], lwd=2, col="blue")
                                                                abline(v=angax[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=cavarj.mean[j], lwd=2)
                                                                abline(v=angax[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(anguloejesalm, function(x) x[i,j]))
                                                        dev.off()  
                                                }#end for (j in  1 : 2)
                                                
                                                calc.cavarj <-array(cbind(unlist(angax[i,]), cavarj.mean, se.cavarj, sesgo.cavarj, ic.t.cavarjinf, ic.t.cavarjsup, ic.p.cavarjinf, ic.p.cavarjsup, ic.bca.cavarjinf, ic.bca.cavarjsup),
                                                                    dim=c(2,10))
                                                calc.cavarj <- as.data.frame(calc.cavarj)
                                                colnames(calc.cavarj) <- titulo
                                                rownames(calc.cavarj) <- ejes[1:2]
                                                
                                                cat("\n Angles between ", nombresvariables[i]," and ", "\n",file="temp.txt")					
                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                write.table(round(calc.cavarj, digits=3), file="temp.txt", sep="\t", dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")     
                                        }#end for (i in 1:dim(covartotal)[1])
                                }#end if (cavarjVal=="1")
                                
                                
                                ### Relative contribution to total variability
                                
                                if (ccrtjVal=="1")
                                {
                                        calc.crtj <-c()
                                        crtj.mean <-c()
                                        se.crtj <-c()
                                        sesgo.crtj <-c()
                                        ic.t.crtjinf <-c()
                                        ic.t.crtjsup <-c()
                                        ic.p.crtjinf <-c()
                                        ic.p.crtjsup <-c()
                                        ic.bca.crtjinf <-c()
                                        ic.bca.crtjsup <-c()
                                        
                                        for (i in 1:dim(CRTj)[1])
                                        {
                                                calc.crtj <-cal.ic(sapply(crtjalm, function(x) x[i,1]), liminf, limsup, CRTj[i,1], sapply(crtjjackr, function(x) x[i,1]), niter)
                                                crtj.mean <- c(crtj.mean, calc.crtj[1])
                                                se.crtj <- c(se.crtj,calc.crtj[2])
                                                sesgo.crtj <- c(sesgo.crtj,calc.crtj[3])
                                                ic.t.crtjinf <- c(ic.t.crtjinf,calc.crtj[4])
                                                ic.t.crtjsup <- c(ic.t.crtjsup,calc.crtj[5])
                                                ic.p.crtjinf <- c(ic.p.crtjinf,calc.crtj[6])
                                                ic.p.crtjsup <- c(ic.p.crtjsup,calc.crtj[7])
                                                ic.bca.crtjinf <- c(ic.bca.crtjinf,calc.crtj[8])
                                                ic.bca.crtjsup <- c(ic.bca.crtjsup,calc.crtj[9])
                                                
                                                pdf(paste("Histogram of contribution to total variability of variable", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(crtjalm, function(x) x[i,1]), main="Histogram", xlab=paste("CRT of variable", i))
                                                
                                                if(cpdfVal=="Color pdf")
                                                {
                                                        abline(v=crtj.mean[i], lwd=2, col="blue")
                                                        abline(v=CRTj[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=crtj.mean[i], lwd=2)
                                                        abline(v=CRTj[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(crtjalm, function(x) x[i,1]))
                                                dev.off()
                                                
                                                
                                                postscript(paste("Histogram of contribution to total variability of variable", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(crtjalm, function(x) x[i,1]), main="Histogram", xlab=paste("CRT of variable", i))
                                                
                                                if(cepsVal=="Color eps")
                                                {
                                                        abline(v=crtj.mean[i], lwd=2, col="blue")
                                                        abline(v=CRTj[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=crtj.mean[i], lwd=2)
                                                        abline(v=CRTj[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(crtjalm, function(x) x[i,1]))
                                                dev.off()                                                                               
                                                
                                        }#end for (i in 1:length(crtjalm))
                                        
                                        calc.crtj <-array(cbind(CRTj[,1],crtj.mean, se.crtj, sesgo.crtj, ic.t.crtjinf, ic.t.crtjsup, ic.p.crtjinf, ic.p.crtjsup, ic.bca.crtjinf, ic.bca.crtjsup),
                                                          dim=c(dim(CRTj)[1],10))
                                        calc.crtj <- as.data.frame(calc.crtj)
                                        colnames(calc.crtj) <- titulo
                                        rownames(calc.crtj) <- nombresvariables
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Relative contribution to total variability of the column element j:\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.crtj, digits=3),file="temp.txt", sep="\t",dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (ccrtjVal=="1")
                                
                                ### Variable lengths
                                
                                if (clvVal=="1")
                                {
                                        calc.clv <-c()
                                        clv.mean <-c()
                                        se.clv <-c()
                                        sesgo.clv <-c()
                                        ic.t.clvinf <-c()
                                        ic.t.clvsup <-c()
                                        ic.p.clvinf <-c()
                                        ic.p.clvsup <-c()           
                                        ic.bca.clvinf <-c()
                                        ic.bca.clvsup <-c()           
                                        
                                        for (i in 1:dim(longitudprin)[1])
                                        {
                                                calc.clv <-cal.ic(sapply(longalm, function(x) x[i,1]), liminf, limsup, longitudprin[i,1], sapply(longjackr, function(x) x[i,1]), niter)
                                                clv.mean <- c(clv.mean, calc.clv[1])
                                                se.clv <- c(se.clv,calc.clv[2])
                                                sesgo.clv <- c(sesgo.clv,calc.clv[3])
                                                ic.t.clvinf <- c(ic.t.clvinf,calc.clv[4])
                                                ic.t.clvsup <- c(ic.t.clvsup,calc.clv[5])
                                                ic.p.clvinf <- c(ic.p.clvinf,calc.clv[6])
                                                ic.p.clvsup <- c(ic.p.clvsup,calc.clv[7])
                                                ic.bca.clvinf <- c(ic.bca.clvinf,calc.clv[8])
                                                ic.bca.clvsup <- c(ic.bca.clvsup,calc.clv[9])
                                                
                                                pdf(paste("Histogram of length of variable", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(longalm, function(x) x[i,1]), main="Histogram", xlab=paste("Length of variable", i))
                                                
                                                if(cpdfVal=="Color pdf")
                                                {
                                                        abline(v=clv.mean[i], lwd=2, col="blue")
                                                        abline(v=longitudprin[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=clv.mean[i], lwd=2)
                                                        abline(v=longitudprin[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(longalm, function(x) x[i,1]))
                                                dev.off()
                                                
                                                
                                                postscript(paste("Histogram of length of variable", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(longalm, function(x) x[i,1]), main="Histogram", xlab=paste("Length of variable", i))
                                                
                                                if(cepsVal=="Color eps")
                                                {
                                                        abline(v=clv.mean[i], lwd=2, col="blue")
                                                        abline(v=longitudprin[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=clv.mean[i], lwd=2)
                                                        abline(v=longitudprin[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(longalm, function(x) x[i,1]))
                                                dev.off()                                                                               
                                                
                                        }#end for (i in 1:length(longalm))
                                        
                                        calc.clv <-array(cbind(longitudprin[,1],clv.mean, se.clv, sesgo.clv, ic.t.clvinf, ic.t.clvsup, ic.p.clvinf, ic.p.clvsup, ic.bca.clvinf, ic.bca.clvsup),
                                                         dim=c(dim(longitudprin)[1],10))
                                        calc.clv <- as.data.frame(calc.clv)
                                        colnames(calc.clv) <- titulo
                                        rownames(calc.clv) <- nombresvariables
                                        
                                        
                                        
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Variable lengths:\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.clv, digits=3), file="temp.txt", sep="\t", dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (clvVal=="1")
                                
                                
                                
                                
                                ### Relative contribution to total variability
                                
                                if (ccrtiVal=="1")
                                {
                                        calc.crti <-c()
                                        crti.mean <-c()
                                        se.crti <-c()
                                        sesgo.crti <-c()
                                        ic.t.crtiinf <-c()
                                        ic.t.crtisup <-c()
                                        ic.p.crtiinf <-c()
                                        ic.p.crtisup <-c()
                                        ic.bca.crtiinf <-c()
                                        ic.bca.crtisup <-c()
                                        
                                        
                                        for (i in 1:dim(CRTi)[1])
                                        {
                                                
                                                calc.crti <-cal.ic(crtialm[[i]], liminf, limsup, CRTi[i,1], crtijackr[[i]], niter)
                                                crti.mean <- c(crti.mean, calc.crti[1])
                                                se.crti <- c(se.crti,calc.crti[2])
                                                sesgo.crti <- c(sesgo.crti,calc.crti[3])
                                                ic.t.crtiinf <- c(ic.t.crtiinf,calc.crti[4])
                                                ic.t.crtisup <- c(ic.t.crtisup,calc.crti[5])
                                                ic.p.crtiinf <- c(ic.p.crtiinf,calc.crti[6])
                                                ic.p.crtisup <- c(ic.p.crtisup,calc.crti[7])
                                                ic.bca.crtiinf <- c(ic.bca.crtiinf,calc.crti[8])
                                                ic.bca.crtisup <- c(ic.bca.crtisup,calc.crti[9])
                                                
                                                pdf(paste("Histogram of contribution to total variability of individual", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(crtialm, function(x) x[[i]]), main="Histogram", xlab=paste("CRT of individual", i))
                                                
                                                if(cpdfVal=="Color pdf")
                                                {
                                                        abline(v=crti.mean[i], lwd=2, col="blue")
                                                        abline(v=CRTi[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=crti.mean[i], lwd=2)
                                                        abline(v=CRTi[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(crtialm, function(x) x[[i]]))
                                                dev.off()
                                                
                                                
                                                postscript(paste("Histogram of contribution to total variability of individual", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                par(mfrow=c(1,2))
                                                hist(sapply(crtialm, function(x) x[[i]]), main="Histogram", xlab=paste("CRT of individual", i))
                                                
                                                if(cepsVal=="Color eps")
                                                {
                                                        abline(v=crti.mean[i], lwd=2, col="blue")
                                                        abline(v=CRTi[i,1], lty =2, lwd=2, col="red")
                                                }else{
                                                        abline(v=crti.mean[i], lwd=2)
                                                        abline(v=CRTi[i,1], lty =2, lwd=2)
                                                }        
                                                qqnorm(sapply(crtialm, function(x) x[[i]]))
                                                dev.off()                                                                               
                                                
                                        }#end for (i in 1:length(crtialm))
                                        
                                        calc.crti <-array(cbind(CRTi[,1],crti.mean, se.crti, sesgo.crti, ic.t.crtiinf, ic.t.crtisup, ic.p.crtiinf, ic.p.crtisup, ic.bca.crtiinf, ic.bca.crtisup),
                                                          dim=c(dim(CRTi)[1],10))
                                        calc.crti <- as.data.frame(calc.crti)
                                        colnames(calc.crti) <- titulo
                                        rownames(calc.crti) <- textindividuos
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Relative contribution to total variability of the row element i:\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                        write.table(round(calc.crti, digits=3),file="temp.txt", sep="\t",dec=",")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                }#end if (ccrtiVal=="1")
                                
                                
                                ### Relative contribution of element j to factor q
                                
                                if (ccrejfqVal=="1")
                                {
                                        cat("\n",file="temp.txt")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                        cat("Relative contribution of the column element j to the q-th factor:",file="temp.txt")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                        
                                        for(i in 1: dim(CREjFq)[1])
                                        {
                                                calc.crejfq <-c()
                                                crejfq.mean <-c()
                                                se.crejfq <-c()
                                                sesgo.crejfq <-c()
                                                ic.t.crejfqinf <-c()
                                                ic.t.crejfqsup <-c()
                                                ic.p.crejfqinf <-c()
                                                ic.p.crejfqsup <-c()
                                                ic.bca.crejfqinf <-c()
                                                ic.bca.crejfqsup <-c()
                                                
                                                for (j in  1:dim(CREjFq)[2])
                                                {
                                                        calc.crejfq <-cal.ic(sapply(crejfqalm, function(x) x[i,j]), liminf, limsup, CREjFq[i,j], sapply(crejfqjackr, function(x) x[i,j]), niter)
                                                        crejfq.mean <- c(crejfq.mean, calc.crejfq[1])
                                                        se.crejfq <- c(se.crejfq,calc.crejfq[2])
                                                        sesgo.crejfq <- c(sesgo.crejfq,calc.crejfq[3])
                                                        ic.t.crejfqinf <- c(ic.t.crejfqinf,calc.crejfq[4])
                                                        ic.t.crejfqsup <- c(ic.t.crejfqsup,calc.crejfq[5])
                                                        ic.p.crejfqinf <- c(ic.p.crejfqinf,calc.crejfq[6])
                                                        ic.p.crejfqsup <- c(ic.p.crejfqsup,calc.crejfq[7])
                                                        ic.bca.crejfqinf <- c(ic.bca.crejfqinf,calc.crejfq[8])
                                                        ic.bca.crejfqsup <- c(ic.bca.crejfqsup,calc.crejfq[9])
                                                        
                                                        
                                                        pdf(paste("Histogram of CREjFq of ", nombresvariables[i]," to axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(crejfqalm, function(x) x[i,j]), main="Histogram", xlab=paste("CREjFq of ", nombresvariables[i],"\n to axis ", j))
                                                        
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=crejfq.mean[j], lwd=2, col="blue")
                                                                abline(v=CREjFq[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=crejfq.mean[j], lwd=2)
                                                                abline(v=CREjFq[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(crejfqalm, function(x) x[i,j]))
                                                        dev.off()
                                                        
                                                        
                                                        postscript(paste("Histogram of CREjFq of ", nombresvariables[i]," to axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(crejfqalm, function(x) x[i,j]), main="Histogram", xlab=paste("CREjFq of ", nombresvariables[i],"\n to axis ", j))
                                                        
                                                        if(cepsVal=="Color eps")
                                                        {
                                                                abline(v=crejfq.mean[j], lwd=2, col="blue")
                                                                abline(v=CREjFq[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=crejfq.mean[j], lwd=2)
                                                                abline(v=CREjFq[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(crejfqalm, function(x) x[i,j]))
                                                        dev.off()
                                                        
                                                }#end for (j in  1:dim(CREjFq)[2])
                                                
                                                calc.crejfq <-array(cbind(unlist(CREjFq[i,]), crejfq.mean, se.crejfq, sesgo.crejfq, ic.t.crejfqinf, ic.t.crejfqsup, ic.p.crejfqinf, ic.p.crejfqsup, ic.bca.crejfqinf, ic.bca.crejfqsup),
                                                                    dim=c(nejes,10))
                                                calc.crejfq <- as.data.frame(calc.crejfq)
                                                colnames(calc.crejfq) <- titulo
                                                rownames(calc.crejfq) <- ejes
                                                
                                                cat("\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                cat(nombresvariables[i],"\n", file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                write.table(round(calc.crejfq, digits=3),file="temp.txt", sep="\t",dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                        }#end for(i in 1: dim(CREjFq)[1])
                                }#end if (ccrejfqVal=="1")
                                
                                ### Relative contribution of  factor q to element j
                                
                                if (ccrfqejVal=="1")
                                {
                                        cat("\n",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Relative contribution of the q-th factor to column element j:",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        
                                        for(i in 1: dim(CRFqEj)[2])
                                        {
                                                calc.crfqej <-c()
                                                crfqej.mean <-c()
                                                se.crfqej <-c()
                                                sesgo.crfqej <-c()
                                                ic.t.crfqejinf <-c()
                                                ic.t.crfqejsup <-c()
                                                ic.p.crfqejinf <-c()
                                                ic.p.crfqejsup <-c()
                                                ic.bca.crfqejinf <-c()
                                                ic.bca.crfqejsup <-c()
                                                
                                                for (j in  1:dim(CRFqEj)[1])
                                                {
                                                        
                                                        calc.crfqej <-cal.ic(sapply(crfqejalm, function(x) x[j,i]), liminf, limsup, CRFqEj[j,i], sapply(crfqejjackr, function(x) x[j,i]), niter)
                                                        crfqej.mean <- c(crfqej.mean, calc.crfqej[1])
                                                        se.crfqej <- c(se.crfqej,calc.crfqej[2])
                                                        sesgo.crfqej <- c(sesgo.crfqej,calc.crfqej[3])
                                                        ic.t.crfqejinf <- c(ic.t.crfqejinf,calc.crfqej[4])
                                                        ic.t.crfqejsup <- c(ic.t.crfqejsup,calc.crfqej[5])
                                                        ic.p.crfqejinf <- c(ic.p.crfqejinf,calc.crfqej[6])
                                                        ic.p.crfqejsup <- c(ic.p.crfqejsup,calc.crfqej[7])
                                                        ic.bca.crfqejinf <- c(ic.bca.crfqejinf,calc.crfqej[8])
                                                        ic.bca.crfqejsup <- c(ic.bca.crfqejsup,calc.crfqej[9])
                                                        
                                                        pdf(paste("Histogram of CRFqEj of axis ", i, " to variable " , nombresvariables[j], ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(crfqejalm, function(x) x[j,i]), main="Histogram", xlab=paste("CRFqEj of axis", i, "\n to variable ", nombresvariables[j]))
                                                        
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=crfqej.mean[j], lwd=2, col="blue")
                                                                abline(v=CRFqEj[j,i], lty =2, lwd=2, col="red")
                                                                
                                                        }else{
                                                                abline(v=crfqej.mean[j], lwd=2)
                                                                abline(v=CRFqEj[j,i], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(crfqejalm, function(x) x[j,i]))
                                                        dev.off()
                                                        
                                                        postscript(paste("Histogram of CRFqEj of axis ", i, " to variable " , nombresvariables[j], ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(sapply(crfqejalm, function(x) x[j,i]), main="Histogram", xlab=paste("CRFqEj of axis", i, "\n to variable ", nombresvariables[j]))
                                                        
                                                        if(cpdfVal=="Color eps")
                                                        {
                                                                abline(v=crfqej.mean[j], lwd=2, col="blue")
                                                                abline(v=CRFqEj[j,i], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=crfqej.mean[j], lwd=2)
                                                                abline(v=CRFqEj[j,i], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(sapply(crfqejalm, function(x) x[j,i]))
                                                        dev.off()
                                                        
                                                        
                                                }#end for (j in  1:dimcreqfjalm)[1])
                                                
                                                calc.crfqej <-array(cbind(unlist(CRFqEj[,i]), crfqej.mean, se.crfqej, sesgo.crfqej, ic.t.crfqejinf, ic.t.crfqejsup, ic.p.crfqejinf, ic.p.crfqejsup, ic.bca.crfqejinf, ic.bca.crfqejsup),
                                                                    dim=c(length(nombresvariables),10))
                                                calc.crfqej <- as.data.frame(calc.crfqej)
                                                colnames(calc.crfqej) <- titulo
                                                rownames(calc.crfqej) <- nombresvariables
                                                
                                                cat("\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                cat(ejes[i],"\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                write.table(round(calc.crfqej, digits=3),file="temp.txt", sep="\t",dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                        }#end for(i in 1: dim(crfqejalm)[2])
                                }#end if (ccrfqejVal=="1")
                                
                                
                                ### Relative contribution of element i to factor q
                                
                                if (ccreifqVal=="1")
                                {
                                        cat("\n",file="temp.txt")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                        cat("Relative contribution of the row element i to the q-th factor:",file="temp.txt")
                                        file.append("Resultsbootstrap.txt","temp.txt")
                                        for(i in 1: dim(CREiFq)[1])
                                        {
                                                calc.creifq <-c()
                                                creifq.mean <-c()
                                                se.creifq <-c()
                                                sesgo.creifq <-c()
                                                ic.t.creifqinf <-c()
                                                ic.t.creifqsup <-c()
                                                ic.p.creifqinf <-c()
                                                ic.p.creifqsup <-c()
                                                ic.bca.creifqinf <-c()
                                                ic.bca.creifqsup <-c()
                                                
                                                for (j in  1:dim(CREiFq)[2])
                                                {
                                                        calc.creifq <-cal.ic(creifqalm[[j]][[i]], liminf, limsup, CREiFq[i,j], creifqjackr[[j]][[i]], niter)
                                                        creifq.mean <- c(creifq.mean, calc.creifq[1])
                                                        se.creifq <- c(se.creifq,calc.creifq[2])
                                                        sesgo.creifq <- c(sesgo.creifq,calc.creifq[3])
                                                        ic.t.creifqinf <- c(ic.t.creifqinf,calc.creifq[4])
                                                        ic.t.creifqsup <- c(ic.t.creifqsup,calc.creifq[5])
                                                        ic.p.creifqinf <- c(ic.p.creifqinf,calc.creifq[6])
                                                        ic.p.creifqsup <- c(ic.p.creifqsup,calc.creifq[7])
                                                        ic.bca.creifqinf <- c(ic.bca.creifqinf,calc.creifq[8])
                                                        ic.bca.creifqsup <- c(ic.bca.creifqsup,calc.creifq[9])
                                                        
                                                        
                                                        pdf(paste("Histogram of CREiFq of ", textindividuos[i]," to axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(creifqalm[[j]][[i]], main="Histogram", xlab=paste("CREiFq of ", textindividuos[i],"\n to axis ", j))
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=creifq.mean[j], lwd=2, col="blue")
                                                                abline(v=CREiFq[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=creifq.mean[j], lwd=2)
                                                                abline(v=CREiFq[i,j], lty =2, lwd=2)
                                                        }       
                                                        
                                                        qqnorm(creifqalm[[j]][[i]])
                                                        dev.off()
                                                        
                                                        
                                                        postscript(paste("Histogram of CREiFq of ", textindividuos[i]," to axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(creifqalm[[j]][[i]], main="Histogram", xlab=paste("CREiFq of ", textindividuos[i],"\n to axis ", j))
                                                        
                                                        if(cepsVal=="Color eps")
                                                        {
                                                                abline(v=creifq.mean[j], lwd=2, col="blue")
                                                                abline(v=CREiFq[i,j], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=creifq.mean[j], lwd=2)
                                                                abline(v=CREiFq[i,j], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(creifqalm[[j]][[i]])
                                                        dev.off()
                                                        
                                                }#end for (j in  1:dim(CREiFq)[2])
                                                
                                                calc.creifq <-array(cbind(unlist(CREiFq[i,]), creifq.mean, se.creifq, sesgo.creifq, ic.t.creifqinf, ic.t.creifqsup, ic.p.creifqinf, ic.p.creifqsup, ic.bca.creifqinf, ic.bca.creifqsup),
                                                                    dim=c(nejes,10))
                                                calc.creifq <- as.data.frame(calc.creifq)
                                                colnames(calc.creifq) <- titulo
                                                rownames(calc.creifq) <- ejes
                                                
                                                cat("\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                cat(textindividuos[i],"\n", file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                write.table(round(calc.creifq, digits=3),file="temp.txt", sep="\t",dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                        }#end for(i in 1: dim(CREiFq)[1])
                                }#end if (ccreifqVal=="1")
                                
                                ### Relative contribution of  factor q to element i
                                
                                if (ccrfqeiVal=="1")
                                {
                                        cat("\n",file="temp.txt")        				
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        cat("Relative contribution of the q-th factor to row element i:",file="temp.txt")					
                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                        
                                        for(i in 1: dim(CRFqEi)[2])
                                        {
                                                calc.crfqei <-c()
                                                crfqei.mean <-c()
                                                se.crfqei <-c()
                                                sesgo.crfqei <-c()
                                                ic.t.crfqeiinf <-c()
                                                ic.t.crfqeisup <-c()
                                                ic.p.crfqeiinf <-c()
                                                ic.p.crfqeisup <-c()
                                                ic.bca.crfqeiinf <-c()
                                                ic.bca.crfqeisup <-c()
                                                
                                                for (j in  1:dim(CRFqEi)[1])
                                                {
                                                        
                                                        calc.crfqei <-cal.ic(crfqeialm[[i]][[j]], liminf, limsup, CRFqEi[j,i], crfqeijackr[[i]][[j]], niter)
                                                        crfqei.mean <- c(crfqei.mean, calc.crfqei[1])
                                                        se.crfqei <- c(se.crfqei,calc.crfqei[2])
                                                        sesgo.crfqei <- c(sesgo.crfqei,calc.crfqei[3])
                                                        ic.t.crfqeiinf <- c(ic.t.crfqeiinf,calc.crfqei[4])
                                                        ic.t.crfqeisup <- c(ic.t.crfqeisup,calc.crfqei[5])
                                                        ic.p.crfqeiinf <- c(ic.p.crfqeiinf,calc.crfqei[6])
                                                        ic.p.crfqeisup <- c(ic.p.crfqeisup,calc.crfqei[7])
                                                        ic.bca.crfqeiinf <- c(ic.bca.crfqeiinf,calc.crfqei[8])
                                                        ic.bca.crfqeisup <- c(ic.bca.crfqeisup,calc.crfqei[9])
                                                        
                                                        pdf(paste("Histogram of CRFqEi of axis ", i, " to individual " , textindividuos[j], ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(crfqeialm[[i]][[j]], main="Histogram", xlab=paste("CRFqEi of axis", i, "\n to individual ", textindividuos[j]))
                                                        
                                                        if(cpdfVal=="Color pdf")
                                                        {
                                                                abline(v=crfqei.mean[j], lwd=2, col="blue")
                                                                abline(v=CRFqEi[j,i], lty =2, lwd=2, col="red")
                                                                
                                                        }else{
                                                                abline(v=crfqei.mean[j], lwd=2)
                                                                abline(v=CRFqEi[j,i], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(crfqeialm[[i]][[j]])
                                                        dev.off()
                                                        
                                                        postscript(paste("Histogram of CRFqEi of axis ", i, " to individual " , textindividuos[j], ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                        par(mfrow=c(1,2))
                                                        hist(crfqeialm[[i]][[j]], main="Histogram", xlab=paste("CRFqEi of axis", i, "\n to individual ", textindividuos[j]))
                                                        
                                                        if(cpdfVal=="Color eps")
                                                        {
                                                                abline(v=crfqei.mean[j], lwd=2, col="blue")
                                                                abline(v=CRFqEi[j,i], lty =2, lwd=2, col="red")
                                                        }else{
                                                                abline(v=crfqei.mean[j], lwd=2)
                                                                abline(v=CRFqEi[j,i], lty =2, lwd=2)
                                                        }        
                                                        qqnorm(crfqeialm[[i]][[j]])
                                                        dev.off()
                                                        
                                                        
                                                }#end for (j in  1:dimcreqfialm)[1])
                                                
                                                calc.crfqei <-array(cbind(unlist(CRFqEi[,i]), crfqei.mean, se.crfqei, sesgo.crfqei, ic.t.crfqeiinf, ic.t.crfqeisup, ic.p.crfqeiinf, ic.p.crfqeisup, ic.bca.crfqeiinf, ic.bca.crfqeisup),
                                                                    dim=c(length(textindividuos),10))
                                                calc.crfqei <- as.data.frame(calc.crfqei)
                                                colnames(calc.crfqei) <- titulo
                                                rownames(calc.crfqei) <- textindividuos
                                                
                                                cat("\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                cat(ejes[i],"\n",file="temp.txt")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                write.table(round(calc.crfqei, digits=3),file="temp.txt", sep="\t",dec=",")
                                                file.append("Resultsbootstrap.txt","temp.txt")
                                        }#end for(i in 1: dim(crfqeialm)[2])
                                }#end if (ccrfqeiVal=="1")
                                
                                
                                file.show("Resultsbootstrap.txt")
                                file.remove("temp.txt")
                        }#end Onaxis <- function()
                        
                        numaxis <- tclVar( 1 )
                        enumaxis <-tkentry(barvp,width="50",textvariable=numaxis)
                        but.axis <-tkbutton(barvp,text="Choose",command=Onaxis, bg= "lightblue", width=10, foreground = "navyblue")
                        tkpack(imgbar, expand="TRUE", fill="both")  
                        
                        tkpack(tklabel(barvp,text="Select the number of axes:"),
                               enumaxis,but.axis,expand="FALSE", side= "left", fill ="both")
                        
                        tkfocus(barvp)
                }#end Graphics <- function()
                
                graphic.button <-tkbutton(framegraphic,text="    Graph    ",command=Graphics, bg= "lightblue", width=20, foreground = "navyblue")
                
                tkpack(tklabel(frametext3,text=""),side="right",expand = "TRUE",fill="both")
                tkpack(tklabel(framet3,text=""),side="right",expand = "TRUE",fill="both")
                tkpack(tklabel(frameok3,text=""),side="right",expand = "TRUE",fill="both")
                
                
                tkpack(tklabel(framehvtitle,text="Graph size"),side="right",expand = "TRUE",fill="both")
                tkpack(tklabel(framehnames,text="Horizontal"),side="right",expand = "TRUE",fill="both")
                tkpack(tklabel(framevnames,text="Vertical"),side="right",expand = "TRUE",fill="both")
                
                
                #####  Textbox to change the size of the graph window #####
                entryvalueh <- tclVar(hescale)
                entryh <-tkentry(framehtext,width="10",textvariable=entryvalueh, bg="white")
                tkbind(entryh, "<Return>",Graphics)
                
                entryvaluev <- tclVar(vescale)
                entryv <-tkentry(framevtext,width="10",textvariable=entryvaluev, bg="white")
                tkbind(entryv, "<Return>",Graphics)
                
                tkpack(entryh, entryv, expand = "TRUE",side="top", fill="both")
                
                tkpack(graphic.button, expand="TRUE", side= "left", fill ="both")
                
                tkpack(framecol11,framecol12, side="left", expand = "TRUE", fill="both")
                tkpack(framename11,framename12, side="left", expand = "TRUE", fill="both")
                tkpack(framecex11,framecex12, side="left", expand = "TRUE", fill="both")
                tkpack(frames11,frames12, side="left", expand = "TRUE", fill="both")
                
                tkpack(framecol21,framecol22, side="left", expand = "TRUE", fill="both")
                tkpack(framename21,framename22, side="left", expand = "TRUE", fill="both")
                tkpack(framecex21,framecex22, side="left", expand = "TRUE", fill="both")
                tkpack(tklabel(frames2, text="Show axes"),cb, expand = "TRUE", side="left",expand="TRUE", fill = "both")
                #tkpack(frames21,frames22, side="left", expand = "TRUE", fill="both")
                
                tkpack(frametext1,framet1,frameok1,framecol1,framename1,framecex1,frames1,expand = "TRUE", fill="both")
                tkpack(frametext2,framet2,frameok2,framecol2,framename2,framecex2,frames2,expand = "TRUE", fill="both")
                tkpack(framehnames, framevnames, expand = "FALSE", fill="both")
                tkpack(framehtext, framevtext, expand = "FALSE", fill="both")
                tkpack(framehvnames, framehvtext,expand = "FALSE",side="left", fill="both")
                tkpack(framehvtitle, framehv, expand = "FALSE", fill="both")
                tkpack(frametext3,framet3,frameok3,framett3auxgs, expand = "FALSE", fill="both")
                
                tkpack(framett1,framett2, framett3, expand = "TRUE",side="left", fill="both")
                tkpack(framett,framegraphic,expand = "TRUE",side="top", fill="y")
        }
        
        OK.butinf <-tkbutton(framewigr,text="   OK   ",command=OnOKinf, bg= "lightblue", width=20, foreground = "navyblue")
        
        tkpack(OK.butinf, expand="TRUE", side= "left", fill ="both")
        tkpack(framewi21,framewi22, expand = "TRUE",side="left", fill="x")
        tkpack(framewi1,framewi2, expand = "TRUE",side="top", fill="both")
        tkpack(framewi,framewigr, expand = "TRUE",side="top", fill="y")
        
        tkfocus(winfor)
}#end function
