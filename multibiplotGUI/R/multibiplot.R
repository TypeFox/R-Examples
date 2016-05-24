

multibiplot <- function(x, ni)
{
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
        
        
        multibiplotint<-function(matrices, nejes, tipo, filas)
        {
                if (filas==1)
                {
                        dimensiones<-mapply(function(x){dim(x)[1]}, matrices) 
                }else{
                        dimensiones<-mapply(function(x){dim(x)[2]}, matrices)
                }
                
                matrices<-lapply(matrices, function(x){array(unlist(x), dim=dim(x))})
                estand<-function(matriz, desvia)
                {
                        matrizst<-apply(matriz, 1, function(x,y){x/y}, desvia)
                        return(t(matrizst))
                }#end function
                
                
                normalizar<-function(matriz)
                {
                        covar <- cov(matriz)
                        vvpropios <- eigen(covar)
                        matriz <- matriz/vvpropios$values[1]
                        return(matriz)
                }#end function
                
                
                
                if (filas == 1){
                        
                        ###############################################################################                		
                        ####	We center the matrices
                        ###############################################################################
                        
                        matricesst <- lapply(matrices, function(x){transforma(x,tipo="Column centering")})
                        
                        ###############################################################################				
                        ####	We standardize the matrices
                        ###############################################################################
                        
                        desvt<-sqrt(diag(var(do.call(rbind,matrices))))
                        
                        matricesstt <- lapply(matricesst, function(x){estand(x,desvt)})
                        #matricesstt <- mapply(transforma, matrices, tipo="Standardize columns")
                        
                        ###############################################################################        			
                        ####	We normalize the matrices
                        ###############################################################################
                        
                        matricesnor <- lapply(matricesstt, normalizar)
                        Xpon<<-do.call(rbind,matricesnor)
                        
                } else{
                        
                        matricesnor <- lapply(matrices,normalizar)
                        
                        Xpon<<-do.call(cbind,matricesnor)
                        
                }#end if (filas == 1)
                
                ##############################################################################
                #####        	Coordinates
                ##############################################################################
                ejes<<-c()
                for (i in 1:nejes)
                {
                        ejes<<-c(ejes, paste("Axis",i))
                }#end for (i in 1:nejes)
                
                descom<<-La.svd(Xpon,nu=nejes,nv=nejes)
                bonajuste<<-0
                
                if (tipo == "RMP"){
                        coindividuos<<-descom$u%*%diag(descom$d[1:nejes])
                        covariables<<-t(descom$v)				
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        
                        suma2valprop<<-sum((descom$d[1:nejes])^2)
                        sumaRvalprop<<-sum((descom$d)^2)
                        inercia<<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<<-cumsum(inercia)
                        bonajuste<<-(suma2valprop/sumaRvalprop)*100
                        calcol<<-nejes/length(descom$d)*100
                        calfilas<<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo == "RMP")
                
                if (tipo == "CMP"){
                        
                        coindividuos<<-descom$u
                        covariables<<-t(descom$v)%*%diag(descom$d[1:nejes])				
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        
                        suma2valprop<<-sum((descom$d[1:nejes])^2)
                        sumaRvalprop<<-sum((descom$d)^2)
                        inercia<<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<<-cumsum(inercia)
                        bonajuste<<-(suma2valprop/sumaRvalprop)*100
                        calfilas<<-nejes/length(inerciatot)*100
                        calcol<<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo == "CMP")
                
                
                if (tipo=="RCMP")
                {
                        coindividuos<<-descom$u%*%diag(descom$d[1:nejes])
                        covariables<<-t(descom$v)%*%diag(descom$d[1:nejes])
                        
                        ##############################################################################
                        #####		Contributions, goodness of fit and qualities of representation
                        ##############################################################################
                        
                        suma2valprop<<-sum((descom$d[1:nejes])^2)
                        sumaRvalprop<<-sum((descom$d)^2)
                        inercia<<-(descom$d[1:nejes])^2/sumaRvalprop
                        cuminer<<-cumsum(inercia)
                        calcol<<-(suma2valprop/sumaRvalprop)*100
                        calfilas<<-(suma2valprop/sumaRvalprop)*100
                }#end if (tipo=="RCMP")
                
                coindividuosnam<-as.data.frame(coindividuos)
                #rownames(coindividuosnam)<-textindividuos
                colnames(coindividuosnam)<-ejes
                covariablesnam<-as.data.frame(covariables)
                #rownames(covariablesnam)<-textvariables
                colnames(covariablesnam)<-ejes
                
                coindivcuad<-coindividuos^2
                CRTi<-rowSums(coindivcuad)
                CRTi<-(CRTi*1000)/suma2valprop
                CRTi<-as.data.frame(CRTi)
                #rownames(CRTi)<-textindividuos
                
                covarcuad<-covariables^2
                CRTj<-rowSums(covarcuad)
                CRTj<-(CRTj*1000)/suma2valprop
                CRTj<-as.data.frame(CRTj)
                #rownames(CRTj)<-textvariables
                
                
                CREiFq<-array(dim=dim(coindividuos))
                CREjFq<-array(dim=dim(covariables))
                
                CRFqEi<-coindivcuad
                sumaindi<-rowSums(coindivcuad)
                
                CRFqEj<-covarcuad
                sumavar<-rowSums(covarcuad)
                
                for(i in 1:nejes)
                {
                        CREiFq[,i]<-((coindivcuad)[,i]*1000)/((descom$d[i])^2)
                        CREjFq[,i]<-((covarcuad)[,i]*1000)/((descom$d[i])^2)
                        CRFqEi[,i]<-((coindivcuad)[,i]*1000)/(sumaindi)
                        CRFqEj[,i]<-((covarcuad)[,i]*1000)/(sumavar)
                }#end for(i in 1:nejes)
                
                CREiFq<-as.data.frame(CREiFq)
                #rownames(CREiFq)<-textindividuos
                colnames(CREiFq)<-ejes
                
                CREjFq<-as.data.frame(CREjFq)
                #rownames(CREjFq)<-textvariables
                colnames(CREjFq)<-ejes
                
                CRFqEi<-as.data.frame(CRFqEi)
                #rownames(CRFqEi)<-textindividuos
                colnames(CRFqEi)<-ejes
                
                CRFqEj<-as.data.frame(CRFqEj)
                #rownames(CRFqEj)<-textvariables
                colnames(CRFqEj)<-ejes
                
                
                
                if (filas==1){
                        
                        partir<-function(x, start, end)
                        {
                                mat<-x[start:end,]
                                return(mat)
                        }
                        
                        rango<-as.list(dimensiones)
                        final<-cumsum(dimensiones)
                        princi<-as.list(final-dimensiones+1)
                        numerod<-length(dimensiones)
                        listcoindividuos<-rep(list(coindivcuad), numerod)
                        coorpart<-mapply(partir,listcoindividuos, princi, as.list(final), SIMPLIFY=FALSE)
                        
                        CRTt<-mapply(function(x){sum(x)/suma2valprop}, coorpart)
                        CRGtFq<-t(mapply(function(x,lambda){colSums(x)/lambda}, coorpart, rep(list(descom$d[1:nejes]^2),numerod))) 
                        CRFqGt<-t(mapply(function(x){colSums(x)/sum(x)}, coorpart)) 
                        
                }else{
                        
                        partir<-function(x, start, end)
                        {
                                mat<-x[start:end,]
                                return(mat)
                        }
                        
                        rango<-as.list(dimensiones)
                        final<-cumsum(dimensiones)
                        princi<-as.list(final-dimensiones+1)
                        numerod<-length(dimensiones)
                        listcovariables<-rep(list(covarcuad), numerod)
                        coorpart<-mapply(partir,listcovariables, princi, as.list(final), SIMPLIFY=FALSE)
                        
                        CRTt<-mapply(function(x){sum(x)/suma2valprop}, coorpart)
                        CRGtFq<-t(mapply(function(x,lambda){colSums(x)/lambda}, coorpart, rep(list(descom$d[1:nejes]^2),numerod))) 
                        CRFqGt<-t(mapply(function(x){colSums(x)/sum(x)}, coorpart))
                }#end if (filas==1)
                
                CRTt<-(CRTt*1000)
                CRTt<-as.data.frame(CRTt)
                CRGtFq<-CRGtFq*1000
                CRFqGt<-CRFqGt*1000
                colnames(CRGtFq)<-ejes
                colnames(CRFqGt)<-ejes
                
                resultados<- list(ejes=ejes,descom=descom, coindividuos=coindividuos, covariables=covariables, suma2valprop=suma2valprop,
                                  inercia=inercia, cuminer=cuminer, bonajuste=bonajuste, calcol=calcol, calfilas=calfilas, 
                                  coindividuosnam=coindividuosnam, covariablesnam=covariablesnam, CRTi=CRTi, CRTj=CRTj,
                                  CREiFq=CREiFq, CREjFq=CREjFq, CRFqEi=CRFqEi, CRFqEj=CRFqEj, CRTt=CRTt, CRGtFq=CRGtFq, CRFqGt=CRFqGt, valpro=descom$d)
                return(resultados)
        }#end multibiplotint
        
        
        remuestreojack<-function(cadajack,indices)
        {
                indicessep<-t(array(unlist(strsplit(unlist(indices), split="[.]")), dim=c(2,length(unlist(indices)))))
                tablas<-unique(indicessep[,1])
                tablas<-as.list(as.numeric(tablas))
                remuestra<-mapply(function(t){cadajack[[t]][as.numeric(indicessep[which(indicessep[,1]==t),2]),]},tablas, SIMPLIFY=FALSE)
                return(remuestra)
        }#end remuestreojack
        
        remuestreojackcol<-function(cadajack,indices)
        {
                remuestra<-mapply(function(x,y){x[indices,]},cadajack, SIMPLIFY=FALSE)
                return(remuestra)
        }#end remuestreojackcol
        
        resample_boot<-function(X)
        {
                dimen<-length(X)
                indicestabla<-sample(1:dimen, dimen, replace=TRUE)
                Xinter<- X[indicestabla]
                muestra<-function(x)
                {
                        indices <- sample(1:dim(x)[1], replace = T)
                        rem<-x[indices,]
                        return(list(rem,indices))
                }
                Xr <- lapply(Xinter, muestra)
                Xres<-lapply(Xr, function(x){x[[1]]})
                indices<-lapply(Xr, function(x){x[[2]]})
                indicestot<-mapply(function(x,y){paste(x,y, sep=".")}, indicestabla, indices, SIMPLIFY=FALSE)
                return(list(Xres, indicestot, indicestabla))
        }# end resample_boot<-function(X)
        
        
        resample_bootcol<-function(X)
        {
                dimen<-length(X)
                indicestabla<-1:dimen
                indicestot<-sample(1:dim(X[[1]])[1],dim(X[[1]])[1], replace=TRUE)
                
                Xres<-mapply(function(x){x[indicestot,]}, X, SIMPLIFY=FALSE)
                #indicestot<-rep(list(indices), dimen)
                return(list(Xres, indicestot, indicestabla))
        }# end resample_bootcol<-function(X)
        
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
        
        
        if(missing(ni))
        {
                msg<-("ERROR: this function requires two arguments")
                tkmessageBox(message=msg)
                stop(" this function requires two arguments")
        }#end if(missing(ni))
        
        if(sum(ni)!=dim(x)[1] & sum(ni)!=dim(x)[2])
        {
                msg<-("ERROR: dimensions do not match")
                tkmessageBox(message=msg)
                stop(" dimensions do not match")
        }#end  if(sum(ni)!=dim(x)[1] & sum(ni)!=dim(x)[2])
        
        #############################################################################
        #########	libraries
        #############################################################################
        
#         require(tcltk)
#         library(tkrplot)
#         library(tcltk2)
#         library(rgl)
#         library(shapes)
#         library(cluster)
#         library(dendroextras)
#         library(Matrix)
#         tclRequire("BWidget")
#         
        mientorno <- new.env()
        
        symbols <- c("*",".", "o","O","0","+","-","|","%","#")
        tipo<-"RCMP" 
        filas<-1
        nejes<-3  
        dim1<-1
        dim2<-2
        dim3<-3
        NameVal<-NULL
        rbVal<-NULL
        hescale <- "1.5"
        vescale <- "1.5"
        indicei<-NULL
        Namei<-NULL
        Cexi<-1
        NameCexi<-NULL
        colori<-NULL
        simChoicei<-NULL
        Namei<-NULL
        colores<-c()  
        indicev<-NULL
        Namev<-NULL
        NameValv<-NULL
        Cexv<-1
        NameCexv<-NULL  
        colorv<-NULL   
        tipobi<-"Classical Biplot"
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
        descom<-NULL
        inerciatot<-NULL
        msginertia<-NULL
        nejes<-NULL
        cbVal<-NULL
        X<-NULL
        Xpon<-NULL
        colvariables<-NULL
        colindividuos<-NULL
        textvariables<-NULL
        textindividuos<-NULL
        cexvariables<-NULL
        cexindividuos<-NULL
        simvariables<-NULL
        simindividuos<-NULL
        Choicei<-NULL
        Choicev<-NULL
        sumaRvalprop<-NULL
        suma2valprop<-NULL
        inercia<-NULL
        cuminer<-NULL
        bonajuste<-NULL
        ejes<-NULL
        coindividuos<-NULL
        covariables<-NULL
        coindividuosnam<-NULL
        covariablesnam<-NULL
        CRTi<-NULL
        CRTj<-NULL
        CRTt<-NULL
        CREiFq<-NULL
        CREjFq<-NULL
        CRGtFq<-NULL
        CRFqEi<-NULL
        CRFqEj<-NULL
        CRFqGt<-NULL
        
        calcol<-NULL
        calfilas<-NULL
        labelsVec<-NULL
        sizesVec<-NULL
        centro<-NULL
        simbolos<-NULL
        simChoice<-NULL
        proj<-"normal"
        clb <- "normal"
        Choiceproj<- 0
        colvariablesp <- c()
        tit_graph <- "Graph"
        nclust <-"3"
        niterclust <- "10"
        nstart <- "1"
        dim1ant<-1
        dim2ant<-2
        Limix1 <-tclVar("0")
        Limix2 <-tclVar("0")
        Limiy1 <-tclVar("0")
        Limiy2 <-tclVar("0")
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
        niter <- 1000
        alphaic <- 95
        cbonajusteVal <- NULL
        ccalfilasVal <- NULL
        ccalcolVal <- NULL
        ccrtiVal <- NULL
        ccrtjVal <- NULL
        ccrttVal <- NULL
        ccreifqVal <- NULL
        ccrejfqVal <- NULL
        ccrgtfqVal <- NULL
        ccrfqeiVal <- NULL
        ccrfqejVal <- NULL
        ccrfqgtVal <- NULL
        ceigenVal <- NULL
        cpdfVal <- NULL
        cepsVal <- NULL
        ccaVal <- NULL
        
        
        colorescoor <- c("skyblue","red","green","blue","yellow","pink","orange", "navyblue",
                         "violet", "brown", "grey", "navyblue", "darkgreen", "papayawhip", "paleturquoise", "purple",
                         "seagreen", "azure", "coral", "springgreen", "steelblue", "plum", "orchid", 
                         "lemonchiffon", "lavender", "honeydew", "gold", "deeppink", "darksalmon", "darkmagenta")
        
        #############################################################################
        ### Informative window
        #############################################################################
        
        winfor<-tktoplevel()
        tkwm.title(winfor,"MultibiplotGUI")
        fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
        fontFixedWidth <- tkfont.create(family="courier",size=12)
        tkgrid(tklabel(winfor, text="    "))
        tkgrid(tklabel(winfor,text="MULTIBIPLOT:",font=fontHeading, foreground = "blue"))
        tkgrid(tklabel(winfor, text="    "))
        
        rb1 <- tkradiobutton(winfor)
        rb2 <- tkradiobutton(winfor)
        rbValue <- tclVar("Some sets of individuals which have been observed in a single set of variables")
        tkconfigure(rb1,variable=rbValue,value="Some sets of individuals which have been observed in a single set of variables")
        tkconfigure(rb2,variable=rbValue,value="Some sets of variables observed on a single set of individuals")
        tkgrid(tklabel(winfor,text="Some sets of individuals which have been observed in a single set of variables"),rb1)
        tkgrid(tklabel(winfor,text="Some sets of variables observed on a single set of individuals"),rb2)
        tkgrid(tklabel(winfor, text="    "))
        OnOKinf <- function()
        {
                #############################################################################
                #########	Window to enter the number of matrices to analyze
                #############################################################################
                rbVal <<- as.character(tclvalue(rbValue))
                tipobi<-rbVal
                tkdestroy(winfor)
                if(tipobi=="Some sets of individuals which have been observed in a single set of variables")
                {
                        filas<<-1 
                }else{
                        filas<<-0
                }
                NameVal<<-length(ni)	
                
                #############################################################################
                #########	We create vectors with the names of the matrices and of the  
                #########	eigen values
                #############################################################################
                
                niac <- cumsum(ni)
                matrices<-vector("list",NameVal)
                matricesname<-vector("list",NameVal)
                
                
                
                ###############################################################################				
                ####	If the selected option is Some sets of individuals observed in a single 
                ####	 set of variables
                ###############################################################################
                X <<- array(data=unlist(x), dim=dim(x))
                
                
                
                
                if (filas == 1){
                        for (z in 1:NameVal)
                        {
                                if(z==1)
                                {
                                        matrices[[z]]<-X[1:niac[1],]   
                                        matricesname[[z]]<-x[1:niac[1],]   
                                }else{
                                        matrices[[z]]<-X[(niac[z-1]+1):niac[z],]   
                                        matricesname[[z]]<-x[(niac[z-1]+1):niac[z],]   
                                }    
                        }#end for (z in 1:NameVal)
                        
                } else{
                        for (z in 1:NameVal)
                        {
                                if(z==1)
                                {
                                        matrices[[z]]<-X[,1:niac[1]]   
                                }else{
                                        matrices[[z]]<-X[,(niac[z-1]+1):niac[z]]   
                                }    
                        }#end for (z in 1:NameVal)
                }#end if (filas == 1)
                
                ##############################################################################
                ####	We create vectors of the colors
                ##############################################################################
                
                colvariables<<-rep("blue",times = dim(x)[2])		
                colindividuos<<-rep("green",times = dim(x)[1])
                textvariables<<-colnames(x)
                textindividuos<<-rownames(x)
                
                
                ##############################################################################
                ####	We create vectors of the character size
                ##############################################################################
                
                cexvariables<<-rep(1,times = dim(x)[2])		
                cexindividuos<<-rep(1,times = dim(x)[1])
                
                ##############################################################################
                ####	We create vectors of the symbols
                ##############################################################################
                
                simvariables<<-rep(" ",times = dim(x)[2])		
                simindividuos<<-rep("+",times = dim(x)[1])
                
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
                
                
                tkadd(fileMenutt,"command",label="JK-biplot",command=function() tipo<<-"RMP")
                tkadd(fileMenutt,"command",label="HJ-biplot",command=function() tipo<<-"RCMP")
                
                tkadd(topMenutt,"cascade",label="Biplot",menu=fileMenutt)
                
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
                framehvtitle<-tkframe(framett3auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                framehv<-tkframe(framett3auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                framehvnames<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                framehnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                framevnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                
                framehvtext<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                framehtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                framevtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                
                frames2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                framegraphic<-tkframe(tt, relief = "flat", borderwidth = 2, background = "white")
                
                ##### Checkbox to show the axes or not #######
                
                cb <- tkcheckbutton(frames2)
                cbValue <- tclVar("0")
                tkconfigure(cb,variable=cbValue)
                
                ##############################################################################
                ##### 	List of individuals
                ##############################################################################
                
                indicei<-NULL
                Namei<-NULL
                NameVali<-NULL
                Cexi<-1
                NameCexi<-NULL
                scri <- tkscrollbar(framet1, repeatinterval=5, command=function(...)tkyview(tli,...))
                tli<-tklistbox(framet1,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scri,...),background="white")
                tkpack(tklabel(frametext1,text="Individuals"),side="left",expand = "TRUE",fill="both")
                
                for (i in 1:(dim(x)[1]))
                {
                        tkinsert(tli,"end",textindividuos[i])
                }#end for (i in 1:(dim(x)[1]))
                
                tkselection.set(tli,0) #  Indexing starts at zero.
                
                OnOKi <- function()
                {
                        Choicei <<- textindividuos[as.numeric(tkcurselection(tli))+1]
                        
                        ##### Color of the selected variable  #############
                        
                        indicei<<-as.numeric(tkcurselection(tli))+1
                        colori <- colindividuos[indicei[1]]
                        tkconfigure(canvasi,bg=colori)
                        
                        ##### Text of the selected variable  #############
                        
                        Namei <<- tclVar(textindividuos[indicei[1]])
                        tkconfigure(entry.Namei,textvariable=Namei)
                        
                        ##### Size of the selected variable  #############
                        
                        Cexi <<- tclVar(cexindividuos[indicei[1]])
                        tkconfigure(entry.Cexi,textvariable=Cexi)
                }#end OnOKi <- function()
                
                OK.buti <-tkbutton(frameok1,text="    OK    ",command=OnOKi)
                
                tkpack(tli,scri,expand = "TRUE", side="left", fill = "both")
                tkpack.configure(scri,side="left")
                tkpack(OK.buti,expand = "TRUE", side="left", fill = "both")
                
                
                #######Color#######################################
                
                indicei<-as.numeric(tkcurselection(tli))+1
                colori <- colindividuos[indicei[1]]
                canvasi <- tkcanvas(framecol11,width="57",height="20",bg=colori)
                
                ChangeColori <- function()
                {
                        colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colindividuos[indicei[1]],title="Choose a color"))
                        
                        if (nchar(colori)>0)
                        {
                                tkconfigure(canvasi,bg=colori)
                                colindividuos[indicei]<<-colori
                        }#end if (nchar(colori)>0)
                }#end ChangeColori <- function()
                
                ChangeColor.buttoni<- tkbutton(framecol12,text="Change Color",command=ChangeColori,width=4)
                tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")
                
                
                
                ######Labels   ###################################
                
                Namei <- textindividuos[indicei[1]]
                entry.Namei <-tkentry(framename11,width=10,textvariable=Namei)
                
                OnOKli <- function()
                {
                        NameVali <<- tclvalue(Namei)
                        textindividuos[indicei]<<-NameVali
                        
                        #####Values of listbox###############################
                        
                        for (i in 1:dim(x)[1])
                        {
                                tkdelete(tli,0)
                        }#end for (i in 1:dim(x)[1])
                        
                        for (i in 1:(dim(x)[1]))
                        {
                                tkinsert(tli,"end",textindividuos[i])
                        }#end for (i in 1:(dim(x)[1]))
                }#end OnOKli <- function()
                
                OK.butli <-tkbutton(framename12,text=" Change label",command=OnOKli,width=4)
                tkbind(entry.Namei, "<Return>",OnOKli)
                tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
                
                ###### Sizes   ###################################
                
                Cexi <- cexindividuos[indicei[1]]
                entry.Cexi <-tkentry(framecex11,width=10,textvariable=Cexi)
                
                OnOKci <- function()
                {
                        NameCexi <<- tclvalue(Cexi)
                        cexindividuos[indicei]<<-NameCexi
                }#end OnOKci <- function()
                
                OK.butci <-tkbutton(framecex12,text=" Change size",command=OnOKci,width=4)
                tkbind(entry.Cexi, "<Return>",OnOKci)
                tkpack(entry.Cexi,OK.butci,expand = "TRUE", side="left", fill = "both")
                
                
                ######Symbols  ###################################
                
                comboBoxi <- tkwidget(frames11,"ComboBox",editable=FALSE,values=symbols, width=7)
                
                chang.symi <- function()
                {
                        simChoicei <<- symbols[as.numeric(tclvalue(tcl(comboBoxi,"getvalue")))+1]
                        simindividuos[indicei]<<-simChoicei
                }#end chang.symi <- function()
                
                Change.symboli <-tkbutton(frames12,text="   Change symbol   ",command=chang.symi,width=4, height=1)
                tkpack(comboBoxi,Change.symboli,side="left",expand="TRUE", fill="both")
                
                ##### List of variables ###########################
                
                indicev<-NULL
                Namev<-NULL
                NameValv<-NULL
                Cexv<-1
                NameCexv<-NULL
                
                scrv <- tkscrollbar(framet2, repeatinterval=5, command=function(...)tkyview(tlv,...))
                tlv<-tklistbox(framet2,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scrv,...),background="white")
                tkpack(tklabel(frametext2,text="Variables"),side="left",expand = "TRUE",fill="both")
                
                for (i in 1:dim(x)[2])
                {
                        tkinsert(tlv,"end",textvariables[i])
                }#end for (i in 1:dim(x)[2])
                
                tkselection.set(tlv,0) #  Indexing starts at zero.
                
                OnOKv <- function()
                {
                        Choicev <<- textvariables[as.numeric(tkcurselection(tlv))+1]
                        
                        ##### Color of the selected variable  #############
                        
                        indicev<<-as.numeric(tkcurselection(tlv))+1
                        colorv <- colvariables[indicev[1]]
                        tkconfigure(canvasv,bg=colorv)
                        
                        
                        ##### Text of the selected variable  #############
                        
                        Namev <<- tclVar(textvariables[indicev[1]])
                        tkconfigure(entry.Namev,textvariable=Namev)
                        
                        ##### Size of the selected variable  #############
                        
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
                        NameValv <<- tclvalue(Namev)
                        textvariables[indicev]<<-NameValv
                        
                        #####Values of listbox###############################
                        
                        for (i in 1:(dim(x)[2]))
                        {
                                tkdelete(tlv,0)
                        }#end for (i in 1:(dim(x)[2]))
                        
                        for (i in 1:(dim(x)[2]))
                        {
                                tkinsert(tlv,"end",textvariables[i])
                        }#end for (i in 1:(dim(x)[2]))
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
                        cexvariables[indicev]<<-NameCexv
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
                                descom<<-multibiplotint(matrices, rankMatrix(x)[1], tipo, filas)$descom
                                sumaRvalprop<<-sum((descom$d)^2)
                                inerciatot<<-(descom$d[1:length(descom$d)])^2/sumaRvalprop
                                barplot(descom$d, col="blue", xlab="", ylab="", names.arg=round(inerciatot, digits=2))
                        }#end plotbar<-function()
                        
                        imgbar <<- tkrplot(barvp,fun=plotbar,hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                        msginertia<-"Proportion of inertia explained by each axis:"
                        for (i in 1:length(descom$d))
                        {
                                msginertia<-paste(msginertia, "\n",i, "\t", round(inerciatot[i]*100, digits=2), "%")
                        }#end for (i in 1:length(descom$d))
                        
                        tk2tip(imgbar, msginertia)
                        
                        Onaxis <- function()
                        {
                                nejes <<- tclvalue(numaxis)
                                nejes<<-as.numeric(nejes)
                                if (nejes > length(descom$d))
                                {
                                        msg <- paste("The maximum number of dimensions is ",length(descom$d))
                                        tkmessageBox(message=msg)
                                }else{
                                        tkdestroy(barvp)
                                        nejes <<- as.integer(nejes)
                                        
                                        resultados<- multibiplotint(matrices, nejes, tipo, filas)
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
                                        coindividuosnam<<-resultados$coindividuosnam
                                        covariablesnam<<- resultados$covariablesnam
                                        CRTi<<-resultados$CRTi
                                        CRTj<<-resultados$CRTj
                                        CRTt<<-resultados$CRTt
                                        CREiFq<<-resultados$CREiFq
                                        CREjFq<<-resultados$CREjFq
                                        CRGtFq<<-resultados$CRGtFq
                                        CRFqEi<<-resultados$CRFqEi
                                        CRFqEj<<-resultados$CRFqEj
                                        CRFqGt<<-resultados$CRFqGt
                                        
                                        
                                        rownames(coindividuosnam)<<-textindividuos
                                        rownames(covariablesnam)<<-textvariables
                                        rownames(CRTi)<<-textindividuos
                                        rownames(CRTj)<<-textvariables
                                        rownames(CREiFq)<<-textindividuos
                                        rownames(CREjFq)<<-textvariables
                                        rownames(CRFqEi)<<-textindividuos
                                        rownames(CRFqEj)<<-textvariables
                                        
                                        
                                        cat("File saved in:    ",file="Results.txt")
                                        cat(getwd(),file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")		
                                        cat("CONTRIBUTIONS:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        
                                        if (tipo != "RCMP"){
                                                cat("\n",file="temp.txt")        				
                                                file.append("Results.txt","temp.txt")	
                                                cat("Goodnes of Fit:  ",file="temp.txt")					
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
                                        cat("Relative contribution of the set t to total variability: \n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRTt, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the row element i to the factor q-th:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CREiFq, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the column element j to the factor q-th:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CREjFq, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")        				
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the set t to the factor q: \n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRGtFq, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the factor q-th to row element i:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRFqEi, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the factor q-th to column element j:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRFqEj, digits=2),file="temp.txt", sep="\t",dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Relative contribution of the factor q-th to the set t: \n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(CRFqGt, digits=2), file="temp.txt", sep="\t", dec=",")
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
                                        cat("Variables coordinates:\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")					
                                        write.table(round(covariablesnam, digits=2), file="temp.txt", sep="\t", dec=",")
                                        file.append("Results.txt","temp.txt")
                                        
                                        cat("\n",file="temp.txt")					
                                        file.append("Results.txt","temp.txt")	
                                        cat("Eigen values: \n",file="temp.txt")					
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
                                        
                                        
                                        indexLabeled<-c(1:length(xCoords))
                                        indexLabeledaux<-c()
                                        labeledPoints <- list()
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
                                                        colvariablesp <- rep("grey",dim(x)[2])
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
                                        }# end plotFunctiond
                                        
                                        
                                        
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
                                                        
                                                        for (i in 1:dim(x)[2])
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
                                                        
                                                        points3d(xCoords,yCoords,zCoords, color=colores)
                                                        texts3d(xCoords, yCoords, zCoords,labelsVec,color=colores, cex= as.numeric(sizesVec))
                                                        
                                                        for (i in 1:(dim(covariables)[1]))
                                                        {
                                                                linea<-rbind(covariables[i,c(dim1, dim2, dim3)],c(0,0,0))	
                                                                segments3d(linea[,1],linea[,2], linea[,3],color=colvariables[i])
                                                        }#end for (i in 1:(dim(covariables)[1]))
                                                        
                                                        rgl.bringtotop()
                                                }else{
                                                        msg <- "You have selected less than 3 dimensions. 3D-graph not available"
                                                        tkmessageBox(message=msg)
                                                }#end if (nejes>2)
                                                
                                        }#end g3d<-function()
                                        
                                        bootmulti<-function()
                                        {
                                                wboot<-tktoplevel()
                                                tkwm.title(wboot,"Bootstrap")
                                                #### Frames
                                                
                                                framewi<-tkframe(wboot, relief = "flat", borderwidth = 2, background = "white")
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
                                                framewigr<-tkframe(wboot, relief = "flat", borderwidth = 2, background = "white")
                                                
                                                fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
                                                fontFixedWidth <- tkfont.create(family="courier",size=12)
                                                tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                tkpack(tklabel(framewi1,text="BOOTSTRAP",font=fontHeading, foreground = "blue"), expand = "TRUE", side="left", fill = "both")
                                                tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                
                                                
                                                ######Iterations   ###################################
                                                Niter <- tclVar(niter)
                                                entry.Niter <-tkentry(framewi21i,width=10,textvariable=Niter, bg="white")
                                                tkconfigure(entry.Niter,textvariable=Niter)  		
                                                
                                                tkpack(tklabel(framewi21i, text="Number of resamples"),entry.Niter, expand = "TRUE", side="left", fill = "both")
                                                
                                                ######alpha confidence intervals ###################################
                                                Nalpha <- tclVar(alphaic)
                                                entry.Nalpha <-tkentry(framewi21a,width=10,textvariable=Nalpha, bg="white")
                                                tkconfigure(entry.Nalpha,textvariable=Nalpha)
                                                
                                                tkpack(tklabel(framewi21a, text="Confidence Level     "),entry.Nalpha, expand = "TRUE", side="left", fill = "both")
                                                
                                                
                                                tkpack(framewi21fl, framewi21fr, expand = "TRUE",side="left", fill="both")
                                                tkpack(framewi21ft, framewi21fb, expand = "TRUE",side="top", fill="both")
                                                tkpack(framewi21i, framewi21a, framewi21f, expand = "TRUE",side="top", fill="both")
                                                
                                                
                                                ###### Parameters to estimate   ###################################
                                                tkpack(tklabel(framewi221l, text="Calculate confidence intervals for:"), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                tkpack(tklabel(framewi221c, text=" "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                
                                                ##### Checkbox bonajuste  #######
                                                
                                                cbonajuste <- tkcheckbutton(framewi222c)
                                                cbonajusteValue <- tclVar("0")
                                                tkconfigure(cbonajuste,variable=cbonajusteValue)
                                                
                                                ##### Checkbox calfilas  #######
                                                
                                                ccalfilas <- tkcheckbutton(framewi223c)
                                                ccalfilasValue <- tclVar("0")
                                                tkconfigure(ccalfilas,variable=ccalfilasValue)
                                                
                                                ##### Checkbox calcol  #######
                                                
                                                ccalcol <- tkcheckbutton(framewi224c)
                                                ccalcolValue <- tclVar("0")
                                                tkconfigure(ccalcol,variable=ccalcolValue)
                                                
                                                ##### Checkbox CRTi  #######
                                                
                                                ccrti <- tkcheckbutton(framewi225c)
                                                ccrtiValue <- tclVar("0")
                                                tkconfigure(ccrti, variable=ccrtiValue)
                                                
                                                ##### Checkbox CRTj  #######
                                                
                                                ccrtj <- tkcheckbutton(framewi226c)
                                                ccrtjValue <- tclVar("0")
                                                tkconfigure(ccrtj,variable=ccrtjValue)
                                                
                                                ##### Checkbox CRTt #######
                                                
                                                ccrtt <- tkcheckbutton(framewi227c)
                                                ccrttValue <- tclVar("0")
                                                tkconfigure(ccrtt, variable=ccrttValue)
                                                
                                                ##### Checkbox CREiFq #######
                                                
                                                ccreifq <- tkcheckbutton(framewi228c)
                                                ccreifqValue <- tclVar("0")
                                                tkconfigure(ccreifq,variable=ccreifqValue)
                                                
                                                ##### Checkbox CREjFq #######
                                                
                                                ccrejfq <- tkcheckbutton(framewi229c)
                                                ccrejfqValue <- tclVar("0")
                                                tkconfigure(ccrejfq, variable=ccrejfqValue)
                                                
                                                ##### Checkbox CRGtFq #######
                                                
                                                ccrgtfq <- tkcheckbutton(framewi2210c)
                                                ccrgtfqValue <- tclVar("0")
                                                tkconfigure(ccrgtfq, variable=ccrgtfqValue)
                                                
                                                ##### Checkbox CRFqEi #######
                                                
                                                ccrfqei <- tkcheckbutton(framewi2211c)
                                                ccrfqeiValue <- tclVar("0")
                                                tkconfigure(ccrfqei, variable=ccrfqeiValue)
                                                
                                                ##### Checkbox CRFqEj #######
                                                
                                                ccrfqej <- tkcheckbutton(framewi2212c)
                                                ccrfqejValue <- tclVar("0")
                                                tkconfigure(ccrfqej, variable=ccrfqejValue)
                                                
                                                ##### Checkbox CRFqGt #######
                                                
                                                ccrfqgt <- tkcheckbutton(framewi2213c)
                                                ccrfqgtValue <- tclVar("0")
                                                tkconfigure(ccrfqgt, variable=ccrfqgtValue)
                                                
                                                ##### Checkbox Eigenvalues #######
                                                
                                                ceigen <- tkcheckbutton(framewi2214c)
                                                ceigenValue <- tclVar("0")
                                                tkconfigure(ceigen,variable=ceigenValue)
                                                
                                                
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
                                                
                                                
                                                
                                                tkpack( tklabel(framewi222l, text="-Goodnes of Fit", anchor="nw"), 
                                                        tklabel(framewi223l, text="-Quality of approximation for rows", anchor="nw"), 
                                                        tklabel(framewi224l, text="-Quality of approximation for columns", anchor="nw"),
                                                        tklabel(framewi225l, text="-Relative contribution to total variability of the row element i", anchor="nw"),
                                                        tklabel(framewi226l, text="-Relative contribution to total variability of the column element j", anchor="nw"),
                                                        tklabel(framewi227l, text="-Relative contribution of the set t to total variability", anchor="nw"),
                                                        tklabel(framewi228l, text="-Relative contribution of the row element i to the factor q-th", anchor="nw"),
                                                        tklabel(framewi229l, text="-Relative contribution of the column element j to the factor q-th", anchor="nw"),
                                                        tklabel(framewi2210l, text="-Relative contribution of the set t to the factor q-th", anchor="nw"),
                                                        tklabel(framewi2211l, text="-Relative contribution of the factor q-th to row element i", anchor="nw"),
                                                        tklabel(framewi2212l, text="-Relative contribution of the factor q-th to column element j", anchor="nw"),
                                                        tklabel(framewi2213l, text="-Relative contribution of the factor q-th to the set t", anchor="nw"),
                                                        tklabel(framewi2214l, text="-Eigenvalues", anchor="nw"),
                                                        expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                
                                                tkpack(cbonajuste, ccalfilas, ccalcol, ccrti,ccrtj,ccrtt,
                                                       ccreifq,ccrejfq, ccrgtfq, ccrfqei, ccrfqej, ccrfqgt,
                                                       ceigen,
                                                       expand = "TRUE",side="top", fill="both")
                                                tkpack(framewi221l,framewi222l,framewi223l,framewi224l,framewi225l,framewi226l,framewi227l,framewi228l,framewi229l,
                                                       framewi2210l,framewi2211l,framewi2212l,
                                                       framewi2213l, framewi2214l, expand = "TRUE",side="top", fill="both")
                                                tkpack(framewi221c,framewi222c,framewi223c,framewi224c,framewi225c,framewi226c,framewi227c,framewi228c,framewi229c,
                                                       framewi2210c,framewi2211c,framewi2212c,
                                                       framewi2213c, framewi2214c, expand = "TRUE",side="top", fill="both")
                                                tkpack(framewi22c, framewi22l, expand = "TRUE",side="left", fill="both")
                                                
                                                OnOKboot<-function()
                                                {
                                                        tkdestroy(wboot) 
                                                        niter <<- tclvalue(Niter)
                                                        alphaic <<- tclvalue(Nalpha)
                                                        cbonajusteVal <<- as.character(tclvalue(cbonajusteValue))
                                                        ccalfilasVal <<- as.character(tclvalue(ccalfilasValue))
                                                        ccalcolVal <<- as.character(tclvalue(ccalcolValue))
                                                        ccrtiVal <<- as.character(tclvalue(ccrtiValue))
                                                        ccrtjVal <<- as.character(tclvalue(ccrtjValue))
                                                        ccrttVal <<- as.character(tclvalue(ccrttValue))
                                                        ccreifqVal <<- as.character(tclvalue(ccreifqValue))
                                                        ccrejfqVal <<- as.character(tclvalue(ccrejfqValue))
                                                        ccrgtfqVal <<- as.character(tclvalue(ccrgtfqValue))
                                                        ccrfqeiVal <<- as.character(tclvalue(ccrfqeiValue))
                                                        ccrfqejVal <<- as.character(tclvalue(ccrfqejValue))
                                                        ccrfqgtVal <<- as.character(tclvalue(ccrfqgtValue))
                                                        ceigenVal <<- as.character(tclvalue(ceigenValue))
                                                        
                                                        cpdfVal <<- as.character(tclvalue(rbpdfValue))
                                                        cepsVal <<- as.character(tclvalue(rbepsValue))
                                                        
                                                        ## muestra bootstrap
                                                        #muestraboot<-remuestreomulti(matrices)
                                                        tablas<-as.list(1:length(matrices))
                                                        indices<-lapply(matrices, function(x){1:dim(x)[1]}) 
                                                        textostotales<-mapply(function(x,y){paste(x,y, sep=".")}, tablas, indices, SIMPLIFY=FALSE)
                                                        textosjack<-unlist(textostotales)
                                                        
                                                        if(filas==1)
                                                        {
                                                                muestrarep<-rep(list(matrices),niter)
                                                                muestrasample<-lapply(muestrarep, resample_boot)
                                                                datosresample<-lapply(muestrasample, function(x)x[[1]])
                                                                textosresample<-lapply(muestrasample, function(x)x[[2]])
                                                                tablasresample<-lapply(muestrasample, function(x)x[[3]])
                                                                bootresult<-mapply(multibiplotint, datosresample, rep(list(nejes),niter), rep(list(tipo),niter), rep(list(filas),niter))
                                                                
                                                                bondadesalm<-bootresult[8,]
                                                                calidadesalm<-bootresult[9,]
                                                                calfilasalm<-bootresult[10,]
                                                                eigenalm<-bootresult[22,]
                                                                crtjalm<-bootresult[14,]
                                                                crejfqalm<-bootresult[16,]
                                                                crfqejalm<-bootresult[18,]
                                                                
                                                                crtiaux<-bootresult[13,]
                                                                indicestab<-sort(unique(unlist(textosresample)))
                                                                tabmat<-cbind(unlist(textosresample), t(array(unlist(strsplit(unlist(textosresample),split="[.]")), dim=c(2,length(unlist(textosresample))))))
                                                                indicestablas<-sort(unique(tabmat[,2]))
                                                                
                                                                crtialm<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(crtiaux))
                                                                
                                                                crttaux<-bootresult[19,]
                                                                crttalm<-lapply(indicestablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(crttaux))
                                                                creifqaux<-bootresult[15,]
                                                                creifqalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        creifqalm[[i]]<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(lapply(creifqaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqeiaux<-bootresult[17,]
                                                                crfqeialm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqeialm[[i]]<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(lapply(crfqeiaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                crgtfqaux<-bootresult[20,]
                                                                crgtfqalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crgtfqalm[[i]]<-lapply(tablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(lapply(crgtfqaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqgtaux<-bootresult[21,]
                                                                crfqgtalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqgtalm[[i]]<-lapply(tablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(lapply(crfqgtaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                
                                                                
                                                                #  jakknife submuestra
                                                                muestratotjack<-lapply(1:length(unlist(textostotales)), function(i){unlist(textostotales)[-i]})
                                                                datosjack<-mapply(remuestreojack, rep(list(matrices), length(muestratotjack)),muestratotjack,SIMPLIFY=FALSE)
                                                                jackresult<-mapply(multibiplotint, datosjack, rep(list(nejes),length(datosjack)), rep(list(tipo),length(datosjack)), rep(list(filas),length(datosjack)))
                                                                
                                                                bondadesjackr<-jackresult[8,]
                                                                calidadesjackr<-jackresult[9,]
                                                                calfilasjackr<-jackresult[10,]
                                                                descomjackr<-jackresult[22,]
                                                                crtjjackr<-jackresult[14,]
                                                                crejfqjackr<-jackresult[16,]
                                                                crfqejjackr<-jackresult[18,]
                                                                
                                                                crtiauxjackr<-jackresult[13,]
                                                                indicestabjackr<-sort(unique(unlist(textosjack)))
                                                                tabmatjackr<-cbind(unlist(textosjack), t(array(unlist(strsplit(unlist(textosjack),split="[.]")), dim=c(2,length(unlist(textosjack))))))
                                                                indicestablasjackr<-sort(unique(tabmatjackr[,2]))
                                                                tablasresamplejackr<-rep(list(tablas), length(datosjack))
                                                                
                                                                crtijackr<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(crtiauxjackr))
                                                                crttauxjackr<-jackresult[19,]
                                                                crttjackr<-lapply(tablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(crttauxjackr))
                                                                creifqauxjackr<-jackresult[15,]
                                                                creifqjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        creifqjackr[[i]]<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(lapply(creifqauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqeiauxjackr<-jackresult[17,]
                                                                crfqeijackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqeijackr[[i]]<-lapply(as.list(textosjack), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(lapply(crfqeiauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                crgtfqauxjackr<-jackresult[20,]
                                                                crgtfqjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crgtfqjackr[[i]]<-lapply(tablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(lapply(crgtfqauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqgtauxjackr<-jackresult[21,]
                                                                crfqgtjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqgtjackr[[i]]<-lapply(tablas, function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(lapply(crfqgtauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                
                                                        }else{
                                                                
                                                                muestrarep<-rep(list(matrices),niter)
                                                                muestrasample<-lapply(muestrarep, resample_bootcol)
                                                                datosresample<-lapply(muestrasample, function(x)x[[1]])
                                                                textosresample<-lapply(muestrasample, function(x)x[[2]])
                                                                tablasresample<-lapply(muestrasample, function(x)x[[3]])
                                                                
                                                                bootresult<-mapply(multibiplotint, datosresample, rep(list(nejes),niter), rep(list(tipo),niter), rep(list(filas),niter))
                                                                
                                                                bondadesalm<-bootresult[8,]
                                                                calidadesalm<-bootresult[9,]
                                                                calfilasalm<-bootresult[10,]
                                                                eigenalm<-bootresult[22,]
                                                                crtjalm<-bootresult[14,]
                                                                crejfqalm<-bootresult[16,]
                                                                crfqejalm<-bootresult[18,]
                                                                
                                                                crtiaux<-bootresult[13,]
                                                                #                                                                 indicestab<-sort(unique(unlist(textosresample)))
                                                                #                                                                 tabmat<-cbind(unlist(textosresample), t(array(unlist(strsplit(unlist(textosresample),split="[.]")), dim=c(2,length(unlist(textosresample))))))
                                                                #                                                                 indicestablas<-sort(unique(tabmat[,2]))
                                                                #                                                                 
                                                                crtialm<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(crtiaux))
                                                                
                                                                crttaux<-bootresult[19,]
                                                                
                                                                crttalm<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(crttaux))
                                                                
                                                                creifqaux<-bootresult[15,]
                                                                creifqalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        creifqalm[[i]]<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(lapply(creifqaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqeiaux<-bootresult[17,]
                                                                crfqeialm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqeialm[[i]]<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(textosresample), unlist(lapply(crfqeiaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                crgtfqaux<-bootresult[20,]
                                                                crgtfqalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crgtfqalm[[i]]<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(lapply(crgtfqaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqgtaux<-bootresult[21,]
                                                                crfqgtalm <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqgtalm[[i]]<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresample), unlist(lapply(crfqgtaux, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                
                                                                
                                                                #  jakknife submuestra
                                                                muestratotjack<-lapply(1:dim(matrices[[1]])[1], function(i){c(1:dim(matrices[[1]])[1])[-i]})
                                                                datosjack<-mapply(remuestreojackcol, rep(list(matrices), length(muestratotjack)),muestratotjack,SIMPLIFY=FALSE)
                                                                jackresult<-mapply(multibiplotint, datosjack, rep(list(nejes),length(datosjack)), rep(list(tipo),length(datosjack)), rep(list(filas),length(datosjack)))
                                                                
                                                                bondadesjackr<-jackresult[8,]
                                                                calidadesjackr<-jackresult[9,]
                                                                calfilasjackr<-jackresult[10,]
                                                                descomjackr<-jackresult[22,]
                                                                crtjjackr<-jackresult[14,]
                                                                crejfqjackr<-jackresult[16,]
                                                                crfqejjackr<-jackresult[18,]
                                                                
                                                                crtiauxjackr<-jackresult[13,]
                                                                #                                                                 indicestabjackr<-sort(unique(unlist(textosjack)))
                                                                #                                                                 tabmatjackr<-cbind(unlist(textosjack), t(array(unlist(strsplit(unlist(textosjack),split="[.]")), dim=c(2,length(unlist(textosjack))))))
                                                                #                                                                 indicestablasjackr<-sort(unique(tabmatjackr[,2]))
                                                                tablasresamplejackr<-rep(list(tablas), length(datosjack))
                                                                
                                                                crtijackr<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(crtiauxjackr))
                                                                crttauxjackr<-jackresult[19,]
                                                                crttjackr<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(crttauxjackr))
                                                                creifqauxjackr<-jackresult[15,]
                                                                creifqjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        creifqjackr[[i]]<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(lapply(creifqauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqeiauxjackr<-jackresult[17,]
                                                                crfqeijackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqeijackr[[i]]<-lapply(as.list(1:dim(matrices[[1]])[1]), function(x, indices, datab) datab[which(indices==x)], unlist(muestratotjack), unlist(lapply(crfqeiauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                crgtfqauxjackr<-jackresult[20,]
                                                                crgtfqjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crgtfqjackr[[i]]<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(lapply(crgtfqauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                crfqgtauxjackr<-jackresult[21,]
                                                                crfqgtjackr <- vector("list",nejes)
                                                                for (i in 1:nejes)
                                                                {
                                                                        crfqgtjackr[[i]]<-lapply(tablasresample[[1]], function(x, indices, datab) datab[which(indices==x)], unlist(tablasresamplejackr), unlist(lapply(crfqgtauxjackr, function(x, ejes) x[,ejes], i)))
                                                                }
                                                                
                                                                
                                                        }
                                                        
                                                        
                                                        ####crear vectores con coordenadas de variables
                                                        coorvarrot<-bootresult[4,]
                                                        coorvarrot<-array(c(unlist(covariables),unlist(coorvarrot)), dim=c(dim(x)[2], nejes, as.numeric(niter)+1))
                                                        
                                                        out.var<-procGPA(coorvarrot, reflect=TRUE, distances=FALSE, pcaoutput=FALSE)
                                                        plot(out.var$rotated[,dim1,], out.var$rotated[,dim2,], type="n", main=paste("Bootstrap Coordinates (Variables)"), xlab=paste("Dimension", dim1), ylab=paste("Dimension", dim2), asp=1/1)
                                                        
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
                                                        
                                                        
                                                        titulo <-c("Obs. Value","Mean","SE","Bias","IC t-boot inf","IC t-boot sup","IC perc inf","IC perc sup","IC BCa inf","IC BCa sup")
                                                        
                                                        
                                                        ### Goodness of fit
                                                        
                                                        if (cbonajusteVal=="1")
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
                                                        
                                                        
                                                        ### Quality of representation rows    
                                                        
                                                        if (ccalfilasVal=="1")
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
                                                        }#end if (ccalfilasVal=="1")
                                                        
                                                        ### Quality of representation     
                                                        
                                                        if (ccalcolVal=="1")
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
                                                                
                                                                calc.ccr <-cal.ic(sapply(calidadesalm,function(x) x[1]), liminf, limsup, calcol, sapply(calidadesjackr,function(x) x[1]), niter)
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
                                                        }#end if (ccalcolVal=="1")
                                                        
                                                        
                                                        
                                                        
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
                                                                rownames(calc.crtj) <- textvariables
                                                                
                                                                cat("\n",file="temp.txt")                			
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                cat("Relative contribution to total variability of the column element j:\n",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                                write.table(round(calc.crtj, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                        }#end if (ccrtjVal=="1")
                                                        
                                                        
                                                        
                                                        ### Relative contribution to total variability of table
                                                        
                                                        if (ccrttVal=="1")
                                                        {
                                                                calc.crtt <-c()
                                                                crtt.mean <-c()
                                                                se.crtt <-c()
                                                                sesgo.crtt <-c()
                                                                ic.t.crttinf <-c()
                                                                ic.t.crttsup <-c()
                                                                ic.p.crttinf <-c()
                                                                ic.p.crttsup <-c()
                                                                ic.bca.crttinf <-c()
                                                                ic.bca.crttsup <-c()
                                                                
                                                                
                                                                for (i in 1:dim(CRTt)[1])
                                                                {
                                                                        
                                                                        calc.crtt <-cal.ic(crttalm[[i]], liminf, limsup, CRTt[i,1], crttjackr[[i]], niter)
                                                                        crtt.mean <- c(crtt.mean, calc.crtt[1])
                                                                        se.crtt <- c(se.crtt,calc.crtt[2])
                                                                        sesgo.crtt <- c(sesgo.crtt,calc.crtt[3])
                                                                        ic.t.crttinf <- c(ic.t.crttinf,calc.crtt[4])
                                                                        ic.t.crttsup <- c(ic.t.crttsup,calc.crtt[5])
                                                                        ic.p.crttinf <- c(ic.p.crttinf,calc.crtt[6])
                                                                        ic.p.crttsup <- c(ic.p.crttsup,calc.crtt[7])
                                                                        ic.bca.crttinf <- c(ic.bca.crttinf,calc.crtt[8])
                                                                        ic.bca.crttsup <- c(ic.bca.crttsup,calc.crtt[9])
                                                                        
                                                                        pdf(paste("Histogram of contribution to total variability of table", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                        par(mfrow=c(1,2))
                                                                        hist(crttalm[[i]], main="Histogram", xlab=paste("CRT of table", i))
                                                                        if(cpdfVal=="Color pdf")
                                                                        {
                                                                                abline(v=crtt.mean[i], lwd=2, col="blue")
                                                                                abline(v=CRTt[i,1], lty =2, lwd=2, col="red")
                                                                        }else{
                                                                                abline(v=crtt.mean[i], lwd=2)
                                                                                abline(v=CRTt[i,1], lty =2, lwd=2)
                                                                        }        
                                                                        qqnorm(crttalm[[i]])
                                                                        dev.off()
                                                                        
                                                                        
                                                                        postscript(paste("Histogram of contribution to total variability of table", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                        par(mfrow=c(1,2))
                                                                        hist(crttalm[[i]], main="Histogram", xlab=paste("CRT of table", i))
                                                                        
                                                                        if(cepsVal=="Color eps")
                                                                        {
                                                                                abline(v=crtt.mean[i], lwd=2, col="blue")
                                                                                abline(v=CRTt[i,1], lty =2, lwd=2, col="red")
                                                                        }else{
                                                                                abline(v=crtt.mean[i], lwd=2)
                                                                                abline(v=CRTt[i,1], lty =2, lwd=2)
                                                                        }        
                                                                        qqnorm(crttalm[[i]])
                                                                        dev.off()                                                                               
                                                                        
                                                                }#end for (i in 1:length(crttalm))
                                                                
                                                                calc.crtt <-array(cbind(CRTt[,1],crtt.mean, se.crtt, sesgo.crtt, ic.t.crttinf, ic.t.crttsup, ic.p.crttinf, ic.p.crttsup, ic.bca.crttinf, ic.bca.crttsup),
                                                                                  dim=c(dim(CRTt)[1],10))
                                                                calc.crtt <- as.data.frame(calc.crtt)
                                                                colnames(calc.crtt) <- titulo
                                                                rownames(calc.crtt) <- unlist(tablas)
                                                                
                                                                cat("\n",file="temp.txt")                			
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                cat("Relative contribution to total variability of the table t:\n",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                                write.table(round(calc.crtt, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                        }#end if (ccrttVal=="1")
                                                        
                                                        
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
                                                                                
                                                                                
                                                                                pdf(paste("Histogram of CREjFq of ", textvariables[i]," to axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(crejfqalm, function(x) x[i,j]), main="Histogram", xlab=paste("CREjFq of ", textvariables[i],"\n to axis ", j))
                                                                                
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
                                                                                
                                                                                
                                                                                postscript(paste("Histogram of CREjFq of ", textvariables[i]," to axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(crejfqalm, function(x) x[i,j]), main="Histogram", xlab=paste("CREjFq of ", textvariables[i],"\n to axis ", j))
                                                                                
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
                                                                        cat(textvariables[i],"\n", file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        write.table(round(calc.crejfq, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end for(i in 1: dim(CREjFq)[1])
                                                        }#end if (ccrejfqVal=="1")
                                                        
                                                        
                                                        
                                                        ### Relative contribution of table t to factor q
                                                        
                                                        if (ccrgtfqVal=="1")
                                                        {
                                                                cat("\n",file="temp.txt")
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                cat("Relative contribution of the table t to the q-th factor:",file="temp.txt")
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                for(i in 1: dim(CRGtFq)[1])
                                                                {
                                                                        calc.crgtfq <-c()
                                                                        crgtfq.mean <-c()
                                                                        se.crgtfq <-c()
                                                                        sesgo.crgtfq <-c()
                                                                        ic.t.crgtfqinf <-c()
                                                                        ic.t.crgtfqsup <-c()
                                                                        ic.p.crgtfqinf <-c()
                                                                        ic.p.crgtfqsup <-c()
                                                                        ic.bca.crgtfqinf <-c()
                                                                        ic.bca.crgtfqsup <-c()
                                                                        
                                                                        for (j in  1:dim(CRGtFq)[2])
                                                                        {
                                                                                calc.crgtfq <-cal.ic(crgtfqalm[[j]][[i]], liminf, limsup, CRGtFq[i,j], crgtfqjackr[[j]][[i]], niter)
                                                                                crgtfq.mean <- c(crgtfq.mean, calc.crgtfq[1])
                                                                                se.crgtfq <- c(se.crgtfq,calc.crgtfq[2])
                                                                                sesgo.crgtfq <- c(sesgo.crgtfq,calc.crgtfq[3])
                                                                                ic.t.crgtfqinf <- c(ic.t.crgtfqinf,calc.crgtfq[4])
                                                                                ic.t.crgtfqsup <- c(ic.t.crgtfqsup,calc.crgtfq[5])
                                                                                ic.p.crgtfqinf <- c(ic.p.crgtfqinf,calc.crgtfq[6])
                                                                                ic.p.crgtfqsup <- c(ic.p.crgtfqsup,calc.crgtfq[7])
                                                                                ic.bca.crgtfqinf <- c(ic.bca.crgtfqinf,calc.crgtfq[8])
                                                                                ic.bca.crgtfqsup <- c(ic.bca.crgtfqsup,calc.crgtfq[9])
                                                                                
                                                                                
                                                                                pdf(paste("Histogram of CRGtFq of table ", i," to axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(crgtfqalm[[j]][[i]], main="Histogram", xlab=paste("CRGtFq of table", i,"\n to axis ", j))
                                                                                if(cpdfVal=="Color pdf")
                                                                                {
                                                                                        abline(v=crgtfq.mean[j], lwd=2, col="blue")
                                                                                        abline(v=CRGtFq[i,j], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=crgtfq.mean[j], lwd=2)
                                                                                        abline(v=CRGtFq[i,j], lty =2, lwd=2)
                                                                                }       
                                                                                
                                                                                qqnorm(crgtfqalm[[j]][[i]])
                                                                                dev.off()
                                                                                
                                                                                
                                                                                postscript(paste("Histogram of CRGtFq of table ", i," to axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(crgtfqalm[[j]][[i]], main="Histogram", xlab=paste("CRGtFq of table", i,"\n to axis ", j))
                                                                                
                                                                                if(cepsVal=="Color eps")
                                                                                {
                                                                                        abline(v=crgtfq.mean[j], lwd=2, col="blue")
                                                                                        abline(v=CRGtFq[i,j], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=crgtfq.mean[j], lwd=2)
                                                                                        abline(v=CRGtFq[i,j], lty =2, lwd=2)
                                                                                }        
                                                                                qqnorm(crgtfqalm[[j]][[i]])
                                                                                dev.off()
                                                                                
                                                                        }#end for (j in  1:dim(CRGtFq)[2])
                                                                        
                                                                        calc.crgtfq <-array(cbind(unlist(CRGtFq[i,]), crgtfq.mean, se.crgtfq, sesgo.crgtfq, ic.t.crgtfqinf, ic.t.crgtfqsup, ic.p.crgtfqinf, ic.p.crgtfqsup, ic.bca.crgtfqinf, ic.bca.crgtfqsup),
                                                                                            dim=c(nejes,10))
                                                                        calc.crgtfq <- as.data.frame(calc.crgtfq)
                                                                        colnames(calc.crgtfq) <- titulo
                                                                        rownames(calc.crgtfq) <- ejes
                                                                        
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat(i,"\n", file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        write.table(round(calc.crgtfq, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end for(i in 1: dim(CRGtFq)[1])
                                                        }#end if (ccrgtfqVal=="1")
                                                        
                                                        
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
                                                                                
                                                                                if(cepsVal=="Color eps")
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
                                                                                
                                                                                pdf(paste("Histogram of CRFqEj of axis ", i, " to variable " , textvariables[j], ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(crfqejalm, function(x) x[j,i]), main="Histogram", xlab=paste("CRFqEj of axis", i, "\n to variable ", textvariables[j]))
                                                                                
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
                                                                                
                                                                                postscript(paste("Histogram of CRFqEj of axis ", i, " to variable " , textvariables[j], ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(crfqejalm, function(x) x[j,i]), main="Histogram", xlab=paste("CRFqEj of axis", i, "\n to variable ", textvariables[j]))
                                                                                
                                                                                if(cepsVal=="Color eps")
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
                                                                                            dim=c(length(textvariables),10))
                                                                        calc.crfqej <- as.data.frame(calc.crfqej)
                                                                        colnames(calc.crfqej) <- titulo
                                                                        rownames(calc.crfqej) <- textvariables
                                                                        
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat(ejes[i],"\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        write.table(round(calc.crfqej, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end for(i in 1: dim(crfqejalm)[2])
                                                        }#end if (ccrfqejVal=="1")
                                                        
                                                        
                                                        
                                                        ### Relative contribution of  factor q to table t
                                                        
                                                        if (ccrfqgtVal=="1")
                                                        {
                                                                cat("\n",file="temp.txt")                			
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                cat("Relative contribution of the q-th factor to table t:",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                
                                                                for(i in 1: dim(CRFqGt)[2])
                                                                {
                                                                        calc.crfqgt <-c()
                                                                        crfqgt.mean <-c()
                                                                        se.crfqgt <-c()
                                                                        sesgo.crfqgt <-c()
                                                                        ic.t.crfqgtinf <-c()
                                                                        ic.t.crfqgtsup <-c()
                                                                        ic.p.crfqgtinf <-c()
                                                                        ic.p.crfqgtsup <-c()
                                                                        ic.bca.crfqgtinf <-c()
                                                                        ic.bca.crfqgtsup <-c()
                                                                        
                                                                        for (j in  1:dim(CRFqGt)[1])
                                                                        {
                                                                                
                                                                                calc.crfqgt <-cal.ic(crfqgtalm[[i]][[j]], liminf, limsup, CRFqGt[j,i], crfqgtjackr[[i]][[j]], niter)
                                                                                crfqgt.mean <- c(crfqgt.mean, calc.crfqgt[1])
                                                                                se.crfqgt <- c(se.crfqgt,calc.crfqgt[2])
                                                                                sesgo.crfqgt <- c(sesgo.crfqgt,calc.crfqgt[3])
                                                                                ic.t.crfqgtinf <- c(ic.t.crfqgtinf,calc.crfqgt[4])
                                                                                ic.t.crfqgtsup <- c(ic.t.crfqgtsup,calc.crfqgt[5])
                                                                                ic.p.crfqgtinf <- c(ic.p.crfqgtinf,calc.crfqgt[6])
                                                                                ic.p.crfqgtsup <- c(ic.p.crfqgtsup,calc.crfqgt[7])
                                                                                ic.bca.crfqgtinf <- c(ic.bca.crfqgtinf,calc.crfqgt[8])
                                                                                ic.bca.crfqgtsup <- c(ic.bca.crfqgtsup,calc.crfqgt[9])
                                                                                
                                                                                pdf(paste("Histogram of CRFqGt of axis ", i, " to table " , j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(crfqgtalm[[i]][[j]], main="Histogram", xlab=paste("CRFqGt of axis", i, "\n to table ",j))
                                                                                
                                                                                if(cpdfVal=="Color pdf")
                                                                                {
                                                                                        abline(v=crfqgt.mean[j], lwd=2, col="blue")
                                                                                        abline(v=CRFqGt[j,i], lty =2, lwd=2, col="red")
                                                                                        
                                                                                }else{
                                                                                        abline(v=crfqgt.mean[j], lwd=2)
                                                                                        abline(v=CRFqGt[j,i], lty =2, lwd=2)
                                                                                }        
                                                                                qqnorm(crfqgtalm[[i]][[j]])
                                                                                dev.off()
                                                                                
                                                                                postscript(paste("Histogram of CRFqGt of axis ", i, " to table " , j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(crfqgtalm[[i]][[j]], main="Histogram", xlab=paste("CRFqGt of axis", i, "\n to table ", j))
                                                                                
                                                                                if(cepsVal=="Color eps")
                                                                                {
                                                                                        abline(v=crfqgt.mean[j], lwd=2, col="blue")
                                                                                        abline(v=CRFqGt[j,i], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=crfqgt.mean[j], lwd=2)
                                                                                        abline(v=CRFqGt[j,i], lty =2, lwd=2)
                                                                                }        
                                                                                qqnorm(crfqgtalm[[i]][[j]])
                                                                                dev.off()
                                                                                
                                                                                
                                                                        }#end for (j in  1:dimcreqfialm)[1])
                                                                        
                                                                        calc.crfqgt <-array(cbind(unlist(CRFqGt[,i]), crfqgt.mean, se.crfqgt, sesgo.crfqgt, ic.t.crfqgtinf, ic.t.crfqgtsup, ic.p.crfqgtinf, ic.p.crfqgtsup, ic.bca.crfqgtinf, ic.bca.crfqgtsup),
                                                                                            dim=c(length(tablas),10))
                                                                        calc.crfqgt <- as.data.frame(calc.crfqgt)
                                                                        colnames(calc.crfqgt) <- titulo
                                                                        rownames(calc.crfqgt) <- unlist(tablas)
                                                                        
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat(ejes[i],"\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        write.table(round(calc.crfqgt, digits=3),file="temp.txt", sep="\t",dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end for(i in 1: dim(crfqgtalm)[2])
                                                        }#end if (ccrfqgtVal=="1")
                                                        
                                                        
                                                        ### Eigenvalues
                                                        
                                                        if (ceigenVal=="1")
                                                        {
                                                                
                                                                calc.ceigen <-c()
                                                                ceigen.mean <-c()
                                                                se.ceigen <-c()
                                                                sesgo.ceigen <-c()
                                                                ic.t.ceigeninf <-c()
                                                                ic.t.ceigensup <-c()
                                                                ic.p.ceigeninf <-c()
                                                                ic.p.ceigensup <-c()
                                                                ic.bca.ceigeninf <-c()
                                                                ic.bca.ceigensup <-c()
                                                                
                                                                
                                                                for (i in 1:length(eigenalm[[1]]))
                                                                { 
                                                                        calc.ceigen <-cal.ic(sapply(eigenalm, function(x) x[i]), liminf, limsup, descom$d[i], sapply(descomjackr, function(x) x[i]), niter)
                                                                        ceigen.mean <- c(ceigen.mean, calc.ceigen[1])
                                                                        se.ceigen <- c(se.ceigen,calc.ceigen[2])
                                                                        sesgo.ceigen <- c(sesgo.ceigen,calc.ceigen[3])
                                                                        ic.t.ceigeninf <- c(ic.t.ceigeninf,calc.ceigen[4])
                                                                        ic.t.ceigensup <- c(ic.t.ceigensup,calc.ceigen[5])
                                                                        ic.p.ceigeninf <- c(ic.p.ceigeninf,calc.ceigen[6])
                                                                        ic.p.ceigensup <- c(ic.p.ceigensup,calc.ceigen[7])
                                                                        ic.bca.ceigeninf <- c(ic.bca.ceigeninf,calc.ceigen[8])
                                                                        ic.bca.ceigensup <- c(ic.bca.ceigensup,calc.ceigen[9])
                                                                        
                                                                        pdf(paste("Histogram of eigenvalue", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                        par(mfrow=c(1,2))
                                                                        hist(sapply(eigenalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                                        
                                                                        if(cpdfVal=="Color pdf")
                                                                        {
                                                                                abline(v=ceigen.mean[i], lwd=2, col="blue")
                                                                                abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                                        }else{
                                                                                abline(v=ceigen.mean[i], lwd=2)
                                                                                abline(v=descom$d[i], lty =2, lwd=2)
                                                                        }        
                                                                        qqnorm(sapply(eigenalm, function(x) x[i]))
                                                                        dev.off()
                                                                        
                                                                        
                                                                        postscript(paste("Histogram of eigenvalue", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                        par(mfrow=c(1,2))
                                                                        hist(sapply(eigenalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                                        
                                                                        if(cepsVal=="Color eps")
                                                                        {
                                                                                abline(v=ceigen.mean[i], lwd=2, col="blue")
                                                                                abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                                        }else{
                                                                                abline(v=ceigen.mean[i], lwd=2)
                                                                                abline(v=descom$d[i], lty =2, lwd=2)
                                                                        }        
                                                                        qqnorm(sapply(eigenalm, function(x) x[i]))
                                                                        dev.off()                                                                               
                                                                }#end for (i in 1:length(eigenalm))
                                                                
                                                                
                                                                calc.ceigen <-array(cbind(descom$d, ceigen.mean, se.ceigen, sesgo.ceigen, ic.t.ceigeninf, ic.t.ceigensup, ic.p.ceigeninf, ic.p.ceigensup, ic.bca.ceigeninf, ic.bca.ceigensup),
                                                                                    dim=c(length(descom$d),10))
                                                                calc.ceigen <- as.data.frame(calc.ceigen)
                                                                colnames(calc.ceigen) <- titulo
                                                                
                                                                nombreseig<-c()
                                                                for (i in 1: length(descom$d))
                                                                {
                                                                        nombreseig <-c(nombreseig, paste("Eigenvalue",i, sep=""))
                                                                }#end for (i in 1: length(descom$d))
                                                                
                                                                rownames(calc.ceigen) <- nombreseig
                                                                
                                                                cat("\n",file="temp.txt")        				
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                cat("Eigenvalues: \n",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")					
                                                                write.table(round(calc.ceigen, digits=3), file="temp.txt", sep="\t", dec=",")
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                        }#end if (ceigenVal=="1")
                                                        
                                                        
                                                        
                                                        file.show("Resultsbootstrap.txt")
                                                        file.remove("temp.txt")
                                                        
                                                }#end OnOKboot
                                                
                                                OK.butboot <-tkbutton(framewigr,text="   OK   ",command=OnOKboot, bg= "lightblue", width=20, foreground = "navyblue")
                                                
                                                tkpack(OK.butboot, expand="TRUE", side= "left", fill ="both")
                                                tkpack(framewi21,framewi22, expand = "TRUE",side="left", fill="x")
                                                tkpack(framewi1,framewi2, expand = "TRUE",side="top", fill="both")
                                                tkpack(framewi,framewigr, expand = "TRUE",side="top", fill="y")
                                                
                                                
                                                
                                        }#end bootmulti
                                        topMenugr <- tkmenu(wgr)
                                        tkconfigure(wgr, menu = topMenugr)
                                        menuFile <- tkmenu(topMenugr, tearoff = FALSE)
                                        menuSaveAs <- tkmenu(topMenugr, tearoff = FALSE)
                                        menu3d <- tkmenu(topMenugr, tearoff = FALSE) 
                                        menuproj <- tkmenu(topMenugr, tearoff = FALSE)                          		
                                        menuboot <- tkmenu(topMenugr, tearoff = FALSE)
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
                                        tkadd(menuboot, "command", label = "Bootstrap", command = function() {bootmulti()})
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
                                        tkadd(topMenugr, "cascade", label = "Bootstrap", menu = menuboot)
                                        
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
                                                squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                indexClosest <<- which.min(squared.Distance)
                                                mm<-tktoplevel() 	
                                                tkwm.title(mm, labelsVec[indexClosest])	
                                                
                                                framemm1<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm2<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm3<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                framemm4<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")    	
                                                
                                                colori <- colores[indexClosest]
                                                canvasi <- tkcanvas(framemm1,width="120",height="20",bg=colori)
                                                
                                                ChangeColori <- function()
                                                {
                                                        
                                                        colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colores[indexClosest],title="Choose a color"))
                                                        
                                                        if (nchar(colori)>0)
                                                        {
                                                                tkconfigure(canvasi,bg=colori)
                                                                colores[indexClosest]<<-colori
                                                                colindividuos<<-colores[1:length(colindividuos)]
                                                                colvariables<<-colores[(length(colindividuos)+1):(length(colindividuos)+length(colvariables))]
                                                        }#end if (nchar(colori)>0)
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end ChangeColori <- function()
                                                
                                                ChangeColor.buttoni<- tkbutton(framemm1,text="Change Color",command=ChangeColori)
                                                tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")
                                                
                                                tclvalue(Namei) <<- labelsVec[indexClosest]
                                                entry.Namei <<-tkentry(framemm2,width="10",textvariable=Namei)
                                                NameVali <<- Namei 
                                                
                                                OnOKli <- function()
                                                {
                                                        NameVali <<- tclvalue(Namei)
                                                        labelsVec[indexClosest]<<-NameVali
                                                        
                                                        textindividuos<<-labelsVec[1:length(textindividuos)]
                                                        textvariables<<-labelsVec[(length(textindividuos)+1):(length(textindividuos)+length(textvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                        
                                                }#end OnOKli <- function()
                                                
                                                OK.butli <-tkbutton(framemm2,text=" Change label",command=OnOKli,width=2)
                                                tkbind(entry.Namei, "<Return>",OnOKli)
                                                tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
                                                
                                                tclvalue(Cexi) <<- sizesVec[indexClosest]
                                                entry.Cexi <<-tkentry(framemm3,width="10",textvariable=Cexi)
                                                NameCexi <<- Cexi 
                                                
                                                OnOKci <- function()
                                                {
                                                        NameCexi <<- tclvalue(Cexi)
                                                        sizesVec[indexClosest]<<-NameCexi
                                                        
                                                        cexindividuos<<-sizesVec[1:length(cexindividuos)]
                                                        cexvariables<<-sizesVec[(length(cexindividuos)+1):(length(cexindividuos)+length(cexvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end OnOKci <- function()
                                                
                                                OK.butci <-tkbutton(framemm3,text=" Change size",command=OnOKci,width=2)
                                                tkbind(entry.Cexi, "<Return>",OnOKci)
                                                tkpack(entry.Cexi,OK.butci,expand = "TRUE", side="left", fill = "both")
                                                
                                                comboBox <- tkwidget(framemm4,"ComboBox",editable=FALSE,values=symbols,width=10, text= simbolos[indexClosest])
                                                
                                                chang.sym <- function()
                                                {
                                                        simChoice <<-symbols[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
                                                        simbolos[indexClosest]<<-simChoice
                                                        
                                                        simindividuos<<-simbolos[1:length(simindividuos)]
                                                        simvariables<<-simbolos[(length(simindividuos)+1):(length(simindividuos)+length(simvariables))]
                                                        
                                                        tkrreplot(img)
                                                        tkdestroy(mm)
                                                }#end chang.sym <- function()
                                                if(indexClosest %in% c(length(simindividuos)+1):(length(simindividuos)+length(simvariables)))
                                                {}else{
                                                        Change.symbol <-tkbutton(framemm4,text="   Change symbol   ",command=chang.sym,width=6)
                                                        tkpack(comboBox,Change.symbol,side="left",expand="TRUE", fill="both")
                                                        
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
                                                mbval<- tkmessageBox(title="Change of label",
                                                                     message=msg,type="yesnocancel",icon="question")
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
                                                width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin <- parPlotSize[1] * width
                                                xMax <- parPlotSize[2] * width
                                                yMin <- parPlotSize[3] * height
                                                yMax <- parPlotSize[4] * height
                                                
                                                rangeX <- usrCoords[2] - usrCoords[1]
                                                rangeY <- usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                
                                                textos[indexClosest,dim1]<<-xPlotCoord
                                                textos[indexClosest,dim2]<<-yPlotCoord
                                                
                                                tkrreplot(img) 
                                        }#end OnLeftClick.move <- function(x,y)
                                        
                                        
                                        OnLeftClick.down <- function(x,y)
                                        {
                                                xClick <- x
                                                yClick <- y
                                                width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin <- parPlotSize[1] * width
                                                xMax <- parPlotSize[2] * width
                                                yMin <- parPlotSize[3] * height
                                                yMax <- parPlotSize[4] * height
                                                
                                                rangeX <- usrCoords[2] - usrCoords[1]
                                                rangeY <- usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                indexClosest <<- which.min(squared.Distance) 
                                                
                                                anteriorx <<- textos[indexClosest,dim1]
                                                anteriory <<- textos[indexClosest,dim2]
                                                
                                        }#end OnLeftClick.down <- function(x,y)
                                        
                                        
                                        
                                        OnRightClick <- function(x,y)
                                        {
                                                xClick <- x
                                                yClick <- y
                                                width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                
                                                xMin <- parPlotSize[1] * width
                                                xMax <- parPlotSize[2] * width
                                                yMin <- parPlotSize[3] * height
                                                yMax <- parPlotSize[4] * height
                                                
                                                rangeX <- usrCoords[2] - usrCoords[1]
                                                rangeY <- usrCoords[4] - usrCoords[3]
                                                
                                                imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                
                                                xClick <- as.numeric(xClick)+0.5
                                                yClick <- as.numeric(yClick)+0.5
                                                yClick <- height - yClick
                                                
                                                xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                
                                                
                                                labelClosestPointd(xClick,yClick,imgXcoords,imgYcoords)
                                                
                                        }#end OnRightClick <- function(x,y)
                                        
                                        tkbind(img, "<B1-Motion>",OnLeftClick.move)
                                        tkbind(img, "<ButtonPress-1>",OnLeftClick.down)
                                        tkbind(img, "<ButtonRelease-1>",OnLeftClick.up)
                                        tkconfigure(img,cursor="pencil")
                                        tkbind(img, "<Button-3>",OnRightClick)
                                        tkconfigure(img,cursor="pencil")
                                }#end if (nejes > length(descom$d))
                        }# end Onaxis <- function()  
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
                
                
                
                
        }#end OnOKinf <- function()
        
        OK.butinf <-tkbutton(winfor,text="   OK   ",command=OnOKinf, bg= "lightblue", width=20, foreground = "navyblue")
        
        tkgrid(OK.butinf)
        tkfocus(winfor)
        
}#end function

