ea1 <-
function (data, design=1, alpha=0.05, list=FALSE, p.adjust=1, plot=2) 
{
    list=ifelse(list==FALSE,1,2)

sk=function(means, df1, QME, nrep, alpha=0.05){
sk1=function(means, df1, QME, nrep, alpha=alpha) {
means=sort(means,decreasing=TRUE)
n=1:(length(means)-1)
n=as.list(n)
f=function(n){list(means[c(1:n)],means[-c(1:n)])}
g=lapply(n, f)
b1=function(x){(sum(g[[x]][[1]])^2)/length(g[[x]][[1]]) + 
(sum(g[[x]][[2]])^2)/length(g[[x]][[2]])-
(sum(c(g[[x]][[1]],g[[x]][[2]]))^2)/length(c(g[[x]][[1]],g[[x]][[2]]))}
p=1:length(g)
values=sapply(p,b1)
minimo=min(values); maximo=max(values)
alfa=(1/(length(means)+df1))*(sum((means-mean(means))^2)+(df1*QME/nrep))
lambda=(pi/(2*(pi-2)))*(maximo/alfa)
vq=qchisq((alpha),lower.tail=FALSE, df=length(means)/(pi-2))
ll=1:length(values); da=data.frame(ll,values); da=da[order(-values),]
ran=da$ll[1]
r=g[[ran]]; r=as.list(r)
i=ifelse(vq>lambda|length(means)==1, 1,2)
means=list(means)
res=list(means, r)
return(res[[i]]) 
}
u=sk1(means, df1, QME, nrep, alpha=alpha)
u=lapply(u, sk1, df1=df1, QME=QME, nrep=nrep, alpha=alpha)
sk2=function(u){
v1=function(...){c(u[[1]])};v2=function(...){c(u[[1]],u[[2]])};v3=function(...){c(u[[1]],u[[2]],u[[3]])} 
v4=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]])}; v5=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]])}
v6=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]])}
v7=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]])}
v8=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]])}
v9=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]],u[[9]])}
v10=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]],u[[9]],u[[10]])}
lv=list(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
l=length(u)
ti=lv[[l]]
u=ti()	
u=lapply(u, sk1, df1=df1, QME=QME, nrep=nrep, alpha=alpha)
return(u)
}
u=sk2(u);u=sk2(u);u=sk2(u);u=sk2(u);u=sk2(u)
u=sk2(u);u=sk2(u);u=sk2(u);u=sk2(u);u=sk2(u)
v1=function(...){c(u[[1]])};v2=function(...){c(u[[1]],u[[2]])};v3=function(...){c(u[[1]],u[[2]],u[[3]])} 
v4=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]])}; v5=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]])}
v6=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]])}
 v7=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]])}
v8=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]])}
v9=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]],u[[9]])}
v10=function(...){c(u[[1]],u[[2]],u[[3]],u[[4]],u[[5]],u[[6]],u[[7]],u[[8]],u[[9]],u[[10]])}
lv=list(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
l=length(u)
ti=lv[[l]]
u=ti()	
rp=u
l2=lapply(rp, length)
l2=unlist(l2)
rp2=rep(letters[1:length(rp)], l2)
return(rp2)
}

    cv <- function(x) {
        sd = (deviance(x)/df.residual(x))^0.5
        mm = mean(fitted(x))
        r = 100 * sd/mm
        return(round(r, 4))
    }
        fr=function(m,data){
        r=resid(m)
	i=ifelse(length(r)>5000, 2,1)
	jr=function(r,aa)r+aa-aa
	jsample=function(r,aa)sample(r,aa)
	rr=list(jr,jsample)
	rr=rr[[i]](r,5000)
        s <- shapiro.test(rr)
        b <- bartlett.test(response ~ treatments, na.action = na.omit, 
                           data = data)
        cvf = cv(m)
        rd=as.data.frame((sort(sqrt(r^2),decreasing=TRUE)))
        rl=as.list(rownames(rd))
        r1=rl[[1]];r2=rl[[2]];r3=rl[[3]]
        d=data.frame(round(s$"p.value",4),round(b$"p.value",4), round(cvf,2),as.numeric(r1),as.numeric(r2),as.numeric(r3)); d=t(d)
        rownames(d)=c("p.value Shapiro-Wilk test","p.value Bartlett test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
        colnames(d)="values"
        return(d)}
    
    fa1=function(a){
        res=a;d=data.frame(res); d=round(d,4); d1=d[,5]; d2=ifelse(d1<0.001, "<0.001", d1); 
        d2=d2[-length(d2)];d2=c(d2,"-"); d=d[,-5];d=data.frame(d,d2);d[is.na(d)] <- "-"
        names(d)=c("df", "type I SS", "mean square", "F value", "p>F")
        return(d)
    }
    
    fa2=function(a){
        res=a; d=data.frame(res); d=data.frame(d[,2],d[,1],d[,1]/d[,2],d[,3],d[,4])
        d=round(d,4); d1=d[,5]; d2=ifelse(d1<0.001, "<0.001", d1); 
        d2=d2[-length(d2)];d2=c(d2,"-"); d=d[,-5];d=data.frame(d,d2);d[is.na(d)] <- "-"
        names(d)=c("df", "type III SS", "mean square", "F value", "p>F"); rownames(d)=rownames(res)
        return(d)
    }
    
    fm=function(ma,dff){
        ma=data.frame(ma,co=ma[,1])
        ma=ma[order(ma[,2], decreasing=TRUE),]
        j=ma[,1];j=as.character(j)
        aux <- combn(j, 2)
        w <- apply(aux, 2, paste, collapse = " - ")
        jj=ma[,2]
        auxj <- combn(jj, 2)
        yi=auxj[1,]-auxj[2,]
        jjj=ma$standard.error^2
        auxjj <- combn(jjj, 2)
        si=sqrt((auxjj[1,]+auxjj[2,])/2)
        yx=yi/si; yx=yx^2; yx=sqrt(yx)
        nmeans=length(ma[,1])
        ft=function(yx, nmeans){1-ptukey(yx,nmeans, dff)}
        st=ft(yx,nmeans)
        st=round(st,4)
        fs=function(ns){s=2:ns;return(s)}
        ns=nmeans:2
        se=sapply(ns, fs)
        ns=unlist(se)
        ssnk=ft(yx, ns)
        ssnk=round(ssnk,4)
        sd=1-ptukey(yx,ns, dff)^(1/(ns-1))
        sd=round(sd,4)
        yxx=yi/(si*sqrt(2))
        vt=1-pt(yxx,dff); vt=vt*2
        vt=round(vt,4)
	lp=list("none","holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")
	pf=p.adjust(vt, lp[[p.adjust]])
        ggs=data.frame(w,round(yi,4),st, ssnk, sd, round(pf,4))
	nam=list("p(t)","p(t)adjust.holm", "p(t)adjust.hochberg", "p(t)adjust.hommel", "p(t)adjust.bonferroni", "p(t)adjust.BH", "p(t)adjust.BY","p(t)adjust.fdr")
        colnames(ggs)=c("pair", "contrast","p(tukey)", "p(snk)", "p(duncan)", nam[[p.adjust]])
        return(ggs)
    }

    ft=function(test, alpha=0.05){
        level=alpha
        tes1=test[,3]
        tes2=test[,4]
        tes3=test[,5]
        tes4=test[,6]
        names(tes1)=test$pair
        names(tes2)=test$pair
        names(tes3)=test$pair
        names(tes4)=test$pair
        tes1=ifelse(tes1<=level,TRUE,FALSE)
        tes2=ifelse(tes2<=level,TRUE,FALSE)
        tes3=ifelse(tes3<=level,TRUE,FALSE)
        tes4=ifelse(tes4<=level,TRUE,FALSE)
        x1=tes1;x2=tes2;x3=tes3;x4=tes4
        inab <- function(x, Letters=c(letters, LETTERS), separator=".", decreasing = decreasing){
            obj_x <- deparse(substitute(x))
            namx <- names(x)
            namx <- gsub(" ", "", names(x))
            if(length(namx) != length(x))
                stop("Names required for ", obj_x)
            split_names <- strsplit(namx, "-")
            stopifnot( sapply(split_names, length) == 2 )
            comps <- t(as.matrix(as.data.frame(split_names)))
            rownames(comps) <- names(x)
            lvls <- unique(as.vector(comps))
            n <- length(lvls)
            lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )
            if( sum(x) == 0 ){                                                       
                ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
                names(ltrs) <- lvls
                colnames(lmat) <- ltrs[1]
                msl <- ltrs
                ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
                return(ret)
            }
            else{
                signifs <- comps[x,,drop=FALSE]
                absorb <- function(m){
                    for(j in 1:(ncol(m)-1)){
                        for(k in (j+1):ncol(m)){
                            if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){               
                                m <- m[,-k, drop=FALSE]
                                return(absorb(m))
                            }
                            else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){           
                                m <- m[,-j, drop=FALSE]
                                return(absorb(m))
                            }
                        }
                    }
                    return(m)
                }
                for( i in 1:nrow(signifs) ){                                           
                    tmpcomp <- signifs[i,]
                    wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               
                    if(any(wassert)){
                        tmpcols <- lmat[,wassert,drop=FALSE]
                        tmpcols[tmpcomp[2],] <- FALSE
                        lmat[tmpcomp[1],wassert] <- FALSE
                        lmat <- cbind(lmat, tmpcols)
                        colnames(lmat) <- get_letters( ncol(lmat), Letters=Letters,
                                                       separator=separator)
                        if(ncol(lmat) > 1){                                               
                            lmat <- absorb(lmat)
                            colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                                           separator=separator )
                        }
                    }
                }
            }
            lmat <- lmat[,order(apply(lmat, 2, sum))]
            lmat <- sweepLetters(lmat)                                                                  
            lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                
            colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                           separator=separator)
            lmat <- lmat[,order(apply(lmat, 2, sum))]                                                  
            lmat <- sweepLetters(lmat)
            lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x)))), 
                                     decreasing = decreasing))]                
            colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                           separator=separator)
            ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
            msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                           
            for( i in 1:nrow(lmat) ){
                msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
                absent <- which(!lmat[i,])
                if( length(absent) < 2 ){
                    if( length(absent) == 0 )
                        next
                    else{
                        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
                    }
                }
                else{
                    msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
                                                             function(x) return(rep( " ",x)) ),
                                                     paste, collapse="") )
                }
            }
            msl <- apply(msl, 1, paste, collapse="")
            names(msl) <- rownames(lmat)
            ret <- list(Letters=ltrs)
            return(ret)
        }
        
        sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){
            stopifnot( all(start.col %in% 1:ncol(mat)) )
            locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))         
            cols <- 1:ncol(mat)
            cols <- cols[c( start.col, cols[-start.col] )]
            if( any(is.na(cols) ) )
                cols <- cols[-which(is.na(cols))]
            for( i in cols){
                tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
                tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        
                one <- which(tmp[,i]==1)
                if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){    
                }
                for( j in one ){                                                    
                    if( locked[j,i] == 1 ){                                           
                        next
                    }
                    chck <- 0
                    lck <- list()
                    for( k in one ){
                        if( j==k ){
                            next
                        }
                        else{                                                           
                            rows <- tmp[c(j,k),]
                            dbl <- rows[1,] & rows[2,]
                            hit <- which(dbl)
                            hit <- hit[-which(hit==i)]
                            dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
                            if( any(dbl) ){
                                chck <- chck + 1
                                lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))     
                            }
                        }
                    }
                    if( (chck == (length(one)-1)) && chck != 0 ){                    
                        for( k in 1:length(lck) ){                                     
                            locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
                            locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
                        }
                        mat[j,i] <- FALSE                                             
                    }
                }
                if(all(mat[,i]==FALSE)){                                          
                    mat <- mat[,-i,drop=FALSE]
                    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
                    return(sweepLetters(mat, Letters=Letters, separator=separator))
                }
            }
            onlyF <- apply(mat, 2, function(x) return(all(!x)))
            if( any(onlyF) ){                                                     
                mat <- mat[,-which(onlyF),drop=FALSE]
                colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
            }
            return( mat )
        }
        
        get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){
            n.complete <- floor(n / length(Letters))       
            n.partial <- n %% length(Letters)              
            lett <- character()
            separ=""
            if( n.complete > 0 ){
                for( i in 1:n.complete ){	
                    lett <- c(lett, paste(separ, Letters, sep="") )
                    separ <- paste( separ, separator, sep="" )
                }
            }
            if(n.partial > 0 )
                lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
            return(lett)
        }
        decreasing=FALSE;jjj1=inab(x1, decreasing = decreasing,); jjj1=jjj1[[1]]
        jjj2=inab(x2, decreasing = decreasing,); jjj2=jjj2[[1]]
        jjj3=inab(x3, decreasing = decreasing,); jjj3=jjj3[[1]]
        jjj4=inab(x4, decreasing = decreasing,); jjj4=jjj4[[1]]
	nam=list("t","t.adjust.holm", "t.adjust.hochberg", "t.adjust.hommel", "t.adjust.bonferroni", "t.adjust.BH", "t.adjust.BY","t.adjust.fdr")
        hgy=data.frame(jjj1,jjj2,jjj3,jjj4); names(hgy)=c("tukey","snk","duncan",nam[[p.adjust]])
        return(hgy)
    }

	pres=function(m){
	r=resid(m)
	r=scale(r)
	t=1:length(r)
	g1=function(r){boxplot(r,col="grey80",ylab="Standardized residuals", main="Box plot for residuals")}
	g2=function(r){plot(r~t, pch="",ylim=c(-4,4), ylab="Standardized residuals", xlab="Sequence data", 		main="Standardized residuals vs Sequence data", axes=FALSE);axis(2,c(-4,-3.5,-3,-2.5,-2,-1,0,1,2,2.5,3,3.5,4));abline(h=2.5, lty=2);abline(h=-2.5,lty=2);abline(h=3.5, lty=2, col=2);abline(h=-3.5,lty=2, col=2); text(2.5,2.7, "2.5 z-score");text(2.5,-2.7, "-2.5 z-score");text(2.5,3.7, "3.5 z-score");text(2.5,-3.7, "-3.5 z-score");text(t,r,labels=1:length(r))}
	a=qqnorm(r,plot.it = FALSE)
	d=data.frame(a1=a$x,a2=a$y,a3=sqrt((a$y)^2),a4=1:length(r))	
	do=d[order(d[,3], decreasing=TRUE),]
d1=do[1,c(1,2)]
d2=do[2,c(1,2)]
d3=do[3,c(1,2)]
n1=as.character(do[1,4])
n2=as.character(do[2,4])
n3=as.character(do[3,4])
	g3=function(r){qqnorm(r, ylab="Standardized residuals", xlab="Theoretical quantiles", main="Standardized residuals vs Theoretical quantiles");qqline(r, col = "grey50");text(d1,n1,adj=-0.5,col=2, cex=0.8);text(d2,n2,adj=-0.5,col=2, cex=0.8);text(d3,n3,adj=-0.5,col=2,cex=0.8)}
	g=list(g1,g2,g3)
	g[[plot]](r)
	}

	    
    f1 = function(data) {
        names(data) = c("treatments", "response")
        data <- data.frame(treatments = factor(data$treatments), 
                           response = data$response)
        m <- aov(response ~ treatments, data = data, contrasts = list(treatments = contr.sum))
        m1 <- aov(response ~ -1 + treatments, data = data, contrasts = list(treatments = contr.sum))
        a <- anova(m)
        a<-fa1(a)
	data2 <- na.omit(data)
	res=fr(m,data)
        mean <- round(coef(m1),4)
        standard.error <- round(sqrt(diag(vcov(m1))),4)
        treatment <- levels(data$treatments)
        ma = data.frame(treatment, mean, standard.error)
        rownames(ma) = NULL;dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
	nrep=length(data2[,1])/nlevels(data2[,1])
	means=mean; names(means)=treatment
	QME=deviance(m)/dff
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l <- list(a, mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Means", "Multiple comparison test", "Residual analysis")
      	pres(m)  
	return(l)
    }	
    
    f2 <- function(data) {
        names(data) = c("treatments", "blocks", "response")
        data <- data.frame(treatments = factor(data$treatments), 
                           blocks = factor(data$blocks), response = data$response)
        m <- aov(response ~ treatments + blocks, data = data, 
                 contrasts = list(treatments = contr.sum, blocks = contr.sum))
        m1 <- aov(response ~ -1 + treatments + blocks, data = data, 
                  contrasts = list(treatments = contr.sum, blocks = contr.sum))
        a <- anova(m)
        data2 <- na.omit(data)
        res=fr(m,data2)
        a2 <- Anova(m, type = 3)
        a3 <- a2[-1, ]
        a3<-fa2(a3)
        adjusted.mean <- round(coef(m1)[c(1:nlevels(data$treatments))],4)
        Standart.Error <- round(sqrt(diag(vcov(m1))),4)
        standard.error <- Standart.Error[1:nlevels(data$treatments)]
        treatment <- levels(data$treatments)
        ma = data.frame(treatment, adjusted.mean, standard.error)
        rownames(ma) = NULL;dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
	nrep=length(data2[,1])/nlevels(data2[,1])
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l <- list(a3,mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)}
    f3 <- function(data) {
        names(data) = c("treatments", "rows", "columns", "response")
        data <- data.frame(treatments = factor(data$treatments), 
                           rows = factor(data$rows), columns = factor(data$columns), 
                           response = data$response)
        m <- aov(response ~ treatments + rows + columns, data = data, 
                 contrasts = list(treatments = contr.sum, rows = contr.sum, 
                                  columns = contr.sum))
        m1 <- aov(response ~ -1 + treatments + rows + columns, data = data, 
                  contrasts = list(treatments = contr.sum, rows = contr.sum, 
                                   columns = contr.sum))
        a <- anova(m)
        data2 <- na.omit(data)
        res=fr(m,data2)
        a2 <- Anova(m, type = 3)
        a3 <- a2[-1, ]
        a3<-fa2(a3)
        adjusted.mean <- round(coef(m1)[c(1:nlevels(data$treatments))],4)
        Standart.Error <- round(sqrt(diag(vcov(m1))),4)
        standard.error <- Standart.Error[1:nlevels(data$treatments)]
        treatment <- levels(data$treatments)
        ma = data.frame(treatment, adjusted.mean, standard.error)
        rownames(ma) = NULL;dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l <- list(a3,mf, test,res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)
    }
    
    f4 = function(data) {
        names(data) = c("treatments", "squares", "rows", "columns", 
                        "response")
        data <- data.frame(treatments = factor(data$treatments), 
                           squares = factor(data$squares), rows = factor(data$rows), 
                           columns = factor(data$columns), response = data$response)
        m <- aov(response ~ treatments + squares + rows + columns, 
                 data = data, contrasts = list(treatments = contr.sum, 
                                               squares = contr.sum, rows = contr.sum, columns = contr.sum))
        m1 <- aov(response ~ -1 + treatments + squares + rows + 
                      columns, data = data, contrasts = list(treatments = contr.sum, 
                                                          squares = contr.sum, rows = contr.sum, columns = contr.sum))
        a <- anova(m)
        a<-fa1(a)
        data2 <- na.omit(data)
        res=fr(m,data2)
        adjusted.mean <- round(coef(m1)[1:nlevels(data$treatments)],4)
        standard.error <- round(sqrt(diag(vcov(m1)))[1:nlevels(data$treatments)],4)
        treatment <- levels(data$treatments)
        ma = data.frame(treatment, adjusted.mean, standard.error)
        rownames(ma) = NULL;dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l <- list(a, mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)
    }
    
    # covariate
    f5<-function(data){ 
        names(data)=c("treatments", "covariate", "response")
        data<-data.frame(treatments=factor(data$treatments), covariate=as.numeric(data$covariate), response=data$response)
        m<-lm(response~ covariate+treatments,data=data, contrasts=list(treatments=contr.sum))
        m1<-lm(response~-1+treatments+covariate,data=data, contrasts=list(treatments=contr.sum))
        a<-anova(m)
        a<-fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        b1=coef(m)[[2]]
        b1
        ra=data$response-(b1*(data$covariate-mean(data$covariate)))
        data=data.frame(data,ra)
        m2=lm(ra~-1+treatments, data=data)
        a2<-Anova(m, type=2) 
        a3<-a2[-1,]
        adjusted.mean<-round(coef(m2)[c(1:nlevels(data$treatments))],4)
        aaa=aggregate(.~treatments,data,FUN=length)
        dff=df.residual(m)
        standard.error<-round(sqrt((deviance(m)/dff)/aaa$ra),4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)}
    
    f6<-function(data){ 
        names(data)=c("treatments", "covariate", "blocks", "response")
        data<-data.frame(treatments=factor(data$treatments), covariate=as.numeric(data$covariate), blocks=factor(data$blocks),response=data$response)
        m<-lm(response~ covariate+treatments+blocks,data=data, contrasts=list(treatments=contr.sum, blocks=contr.sum))
        m1<-lm(response~-1+treatments+covariate+blocks,data=data, contrasts=list(treatments=contr.sum, blocks=contr.sum))
        a<-anova(m)
        a<-fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        b1=coef(m)[[2]]
        b1
        ra=data$response-(b1*(data$covariate-mean(data$covariate)))
        data=data.frame(data,ra)
        m2=lm(ra~-1+treatments+blocks, data=data)
        a2<-Anova(m, type=2) 
        a3<-a2[-1,]
        adjusted.mean<-round(coef(m2)[c(1:nlevels(data$treatments))],4)
        aaa=aggregate(.~treatments,data,FUN=length)
        dff=df.residual(m)
        standard.error<-round(sqrt((deviance(m)/dff)/aaa$ra),4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
	pres(m)  
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
        return(l)}
    
    # incomplete blocks
    f7<-function(data){
        names(data)=c("treatments", "repetition","blocks","response")
        data<-data.frame(treatments=factor(data$treatments), repetition=factor(data$repetition),  blocks=factor(data$blocks), response=data$response)
        m<-lm(terms(response~ repetition/blocks +treatments, keep.order=TRUE),data=data, contrasts=list(treatments=contr.sum, repetition=contr.sum, blocks=contr.sum))
        m1<-lm(terms(response~-1+treatments+repetition/blocks, keep.order=TRUE ),data=data, contrasts=list(treatments=contr.sum, repetition=contr.sum, blocks=contr.sum))
        a<-anova(m)
        a=fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        adjusted.mean<-round(coef(m1)[c(1:nlevels(data$treatments))],4)
        standard.error<-round(sqrt(diag(vcov(m1)))[c(1:nlevels(data$treatments))],4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)}
    
    f8<-function(data){ 
        names(data)=c("treatments", "blocks", "response")
        data<-data.frame(treatments=factor(data$treatments), blocks=as.factor(data$blocks), response=data$response)
        m<-lm(response~ blocks+treatments,data=data, contrasts=list(treatments=contr.sum, blocks=contr.sum))
        m1<-lm(response~-1+treatments+blocks,data=data, contrasts=list(treatments=contr.sum, blocks=contr.sum))
        a<-anova(m)
        data2<-na.omit(data)
        a=fa1(a)
        r=resid(m)
        s=shapiro.test(r)
        cvf=cv(m)
        adjusted.mean<-round(coef(m1)[c(1:nlevels(data$treatments))],4)
        standard.error<-round(sqrt(diag(vcov(m1)))[c(1:nlevels(data$treatments))],4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,s,cvf,mf, test)
        names(l) = list("Analysis of variance",  "Normality test","Coefficient of variation (%)",
                        "Adjusted means", "Multiple comparison test")
		pres(m)  
        return(l)}
    
    f9<-function(data){ 
        names(data)=c("treatments", "subject", "period", "response")
        data<-data.frame(treatments=factor(data$treatments), subject=factor(data$subject), period=as.factor(data$period), response=data$response)
        m<-lm(response~ treatments +subject+period,data=data, contrasts=list(treatments=contr.sum, subject=contr.sum, period=contr.sum))
        m1<-lm(response~-1+treatments+ subject+period, data=data, contrasts=list(treatments=contr.sum, subject=contr.sum, period=contr.sum))
        a<-anova(m)
        a=fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        adjusted.mean<-round(coef(m1)[c(1:nlevels(data$treatments))],4)
        standard.error<-round(sqrt(diag(vcov(m1)))[c(1:nlevels(data$treatments))],4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
        names(l) = list("Analysis of variance",  
                        "Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(m)  
        return(l)}
    
    f10<-function(data){ 
        names(data)=c("treatments","repetition", "blocks", "response")
        data<-data.frame(treatments=factor(data$treatments), repetition=factor(data$repetition), blocks=factor(data$blocks), response=data$response)
        m<-aov(response~repetition/blocks+treatments,data=data, contrasts=list(repetition=contr.sum, blocks=contr.sum, treatments=contr.sum))
	mrr=m
        m1<-aov(response~-1+treatments+repetition/blocks,data=data, contrasts=list(repetition=contr.sum, blocks=contr.sum, treatments=contr.sum))
        a2<-Anova(m, type=3) 
        a3<-a2[-1,]
        a3=fa2(a3)
        data2<-na.omit(data)
        res=fr(m,data2)
        adjusted.mean<-round(coef(m1)[c(1:nlevels(data$treatments))],4)
        standard.error<-round(sqrt(diag(vcov(m1)))[c(1:nlevels(data$treatments))],4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
	means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(ma)=NULL
        QME=deviance(m)/df.residual(m)
        m=nlevels(data$repetition) 
        k=nlevels(data$blocks) 
        vef=QME*(1+(m/((m-1)*(k+1))))
        mb=lm(response~repetition+treatments, data=data)
        QMEB=deviance(mb)/df.residual(mb)
        Ef=100*QMEB/vef
        l<-list(a3,round(vef,4), round(Ef,4), mf,test, res)
        names(l)= list("Analysis of variance", "Effective variance", "Efficiency of the design (%)","Adjusted means", "Multiple comparison test", "Residual analysis")
	pres(mrr)  
        return(l)}
    
    f11<-function(data){ 
        names(data)=c("treatments","repetition", "blocks", "response")
        data<-data.frame(treatments=factor(data$treatments), repetition=factor(data$repetition), blocks=factor(data$blocks), response=data$response)
        block=interaction(data$repetition,data$blocks)
        data=data.frame(data,block)
        m<-lme(response~ repetition +treatments,random=~1|block, data=data, contrasts=list(repetition=contr.sum, treatments=contr.sum), na.action=na.omit)
	mrr=m
        m1<-lme(response~-1+treatments+repetition, random=~1|repetition/blocks,data=data, contrasts=list(repetition=contr.sum, blocks=contr.sum, treatments=contr.sum), na.action=na.omit)
        a3<-anova(m, type="marginal")
        a3<-a3[-1,]
        data2<-na.omit(data)
        r=resid(m)
        s=shapiro.test(r)
        b1=bartlett.test(r~treatments, data=data2)
        adjusted.mean<-round(fixef(m1)[c(1:nlevels(data$treatments))],4)
        standard.error<-round(sqrt(diag(vcov(m1)))[c(1:nlevels(data$treatments))],4)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean,standard.error)
        rownames(ma)=NULL
        dff=a3[[2]][[2]]
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME= m$sigma^2
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(ma)=NULL
        QME= m$sigma^2
        m=nlevels(data$repetition) 
        k=nlevels(data$blocks) 
        mb=lm(response~repetition/blocks+treatments, data=data)
        QMB= anova(mb)[[3]][3]
        vef=QME*(1+(m/((m-1)*(k+1)))*(QMB-QME)/QMB)
        mbb=lm(response~repetition+treatments, data=data)
        Ef=100*(deviance(mbb)/df.residual(mbb))/vef
        l<-list(a3,s,b1, round(vef,4), round(Ef,4), mf,test)
        names(l)= list("Analysis of variance (marginal anova = type III SS)", "Normality test", "Homogeneity of variances", "Effective variance", "Efficiency of the design (%)","Adjusted means", "Multiple comparison test")
		pres(mrr)  
        return(l)}
    
    f12<-function(data){
        names(data)=c("treatments", "period","animal","response")
        data<-data.frame(treatments=factor(data$treatments), period=factor(data$period),  animal=factor(data$animal), response=data$response, p=as.numeric(data$period))
        m=aov(response~p*animal-p+period+animal+treatments, contrasts=list(treatments=contr.sum, animal=contr.sum, period=contr.sum), data=data)
        a<-anova(m)
        a<-fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        c=coef(m)
        names(c)
        class(c)
        c=as.list(c)
        t=c[1:(nlevels(data$animal)+nlevels(data$period)+nlevels(data$treatments)-2)]
        t1=t[-c(1:(nlevels(data$animal)+nlevels(data$period)-1))]
        tf=-(sum(unlist(t1)))
        ef=c(unlist(t1),tf)
        adjusted.mean<- mean(data$response, na.rm = T)+as.numeric(ef)
        Standart.Error<-sqrt(diag(vcov(m))) [1:(nlevels(data$animal)+nlevels(data$period)+nlevels(data$treatments)-2)]
        Standart.Error= Standart.Error[-c(1:(nlevels(data$animal)+nlevels(data$period)-1))]
        standart.error<-c(as.numeric(Standart.Error), as.numeric(Standart.Error)[1])*sqrt(1.5)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean=round(adjusted.mean,4),standard.error=round(standart.error,4))
        rownames(ma)=NULL; dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
        names(l)= list("Analysis of variance", "Adjusted means", "Multiple comparison test","Residual analysis")
	pres(m)  
        return(l)}
    
    f13<-function(data){
        names(data)=c("treatments", "blocks", "period","animal","response")
        data<-data.frame(treatments=factor(data$treatments), blocks=factor(data$blocks), period=factor(data$period),  animal=factor(data$animal), response=data$response, p=as.numeric(data$period))
        m=aov(response~blocks+p*animal-p+animal+treatments+blocks/period, contrasts=list(treatments=contr.sum, animal=contr.sum, period=contr.sum, blocks=contr.sum), data=data)
        a<-anova(m)
        a<-fa1(a)
        data2<-na.omit(data)
        res=fr(m,data2)
        c=coef(m)
        names(c)
        class(c)
        c=as.list(c)
        t=c[1:(a[[1]][1]+ a[[1]][2]+a[[1]][3])+1]
        t1=t[-c(1:(a[[1]][1]+a[[1]][2]))]
        tf=-(sum(unlist(t1)))
        ef=c(unlist(t1),tf)
        adjusted.mean<- mean(data$response, na.rm = T)+as.numeric(ef)
        Standart.Error<-sqrt(diag(vcov(m)))[1:(a[[1]][1]+ a[[1]][2]+a[[1]][3])+1]
        Standart.Error= Standart.Error[-c(1:(a[[1]][1]+a[[1]][2]))]
        standart.error<-c(as.numeric(Standart.Error), as.numeric(Standart.Error)[1])*sqrt(1.5)
        treatment<-levels(data$treatments)
        ma=data.frame(treatment,adjusted.mean=round(adjusted.mean,4),standard.error=round(standart.error,4))
        rownames(ma)=NULL
        dff=df.residual(m)
        test=fm(ma,dff)
        groups=ft(test, alpha); ma=ma[order(ma[,2], decreasing=TRUE),]
        means=adjusted.mean; names(means)=treatment
	QME=deviance(m)/dff
	nrep=length(data2[,1])/nlevels(data2[,1])
	scott_knott=sk(means, dff, QME, nrep, alpha)
        mf=data.frame(ma,groups, scott_knott)
        rownames(mf) = NULL
        l<-list(a,mf, test, res)
        names(l)= list("Analysis of variance", "Adjusted means", "Multiple comparison test","Residual analysis")
	pres(m)  
        return(l)}

f14 <-
function(data, alpha=0.05){
names(data)=c("treatments", "response")
data=data.frame(treatments=factor(data$treatments), response=data$response)
ran=rank(data$response, ties.method = c("average"))
d=data.frame(data,ran)
s = split(d[, -1], d[,1])
fx=function(u){s[[u]][,2]}
u=1:length(s)
ss=lapply(u,fx)
soma=lapply(ss, sum, na.rm=TRUE)
f1=function(x){length(s[x][[1]][[1]])}
x=1:nlevels(data$treatments)
n=lapply(x, f1); n=lapply(n, as.numeric)
f2=function(x){soma[[x]]/n[[x]]}
ra=lapply(x, f2)
names(ra)=names(soma) #ranks
f3=function(x){(soma[[x]]^2)/n[[x]]}
nn=length(data[,1])
S=(1/(nn-1))*((sum(d[,3]^2))-((nn*(nn+1)^2)/4)) #var
ra2=(sum(as.numeric(lapply(x, f3)))-((nn*(nn+1)^2)/4))*(1/S) #T
gl=nlevels(data$treatments)-1
pq=pchisq(ra2,lower.tail=FALSE, df=gl) #chisq
fm=function(x){mean(s[x][[1]][[1]], na.rm=TRUE)}
fmed=function(x){median(s[x][[1]][[1]], na.rm=TRUE)}
means=as.numeric(lapply(x, fm))
md=as.numeric(lapply(x, fmed))
da=data.frame(names=names(s), ra=as.numeric(ra), n=as.numeric(n), means,md, rank=as.numeric(ra))
ma=da[order(da[,2], decreasing=TRUE),]
j=ma[,1];j=as.character(j)
aux <- combn(j, 2)
w <- apply(aux, 2, paste, collapse = " - ")
jj=ma[,2]
auxj <- combn(jj, 2)
yi=auxj[1,]-auxj[2,]#names
jjj=ma[,3]
auxjj <- combn(jjj, 2)
yii=(1/auxjj[1,])+(1/auxjj[2,])#means for ranks
ct=S*((nn-1-ra2)/(nn-nlevels(data$treatments)))
x=1:length(yi)
f4=function(x){yi[x]/sqrt(ct*yii[x])}
tcal=lapply(x, f4)
p=function(x){(1-pt(tcal[[x]],df=nn-nlevels(data$treatments)))*2}
p=lapply(x,p)
pa1=p.adjust(as.numeric(p), method=c("holm"))
pa2=p.adjust(as.numeric(p), method=c("bonferroni"))
pa3=p.adjust(as.numeric(p), method=c("fdr"))
resp=data.frame(w,round(yi,4),round(as.numeric(tcal),4), round(as.numeric(p),4), round(pa1,4), round(pa2,4), round(pa3,4))
colnames(resp)=c("pair", "contrast", "tcal", "p(t)" ,"p.adj(Holm)", "p.adj(Bonferroni)", "p.adj(fdr)")
med=ft(resp[,-2], alpha=alpha)
colnames(med)=c("t", "adjust.Holm", "adjust.Bonferroni", "adjust.fdr" ) 
fxx=function(u){s[[u]][,1]}
sr=lapply(u,fxx)
med=data.frame(treatment=rownames(med),rank=round(ma$rank,4),mean=round(ma$means,4), median=round(ma$md,4),med)
rownames(med)=NULL
med
pri=data.frame(round(c(ra2,pq),4));colnames(pri)="Estimates"
rownames(pri)=c("Kruskal-Wallis chi-squared = ","p.value = ")
l=list(pri,med,resp) ; names(l)=c("Kruskal-Wallis Rank Sum Test", "Ranks, Means and Medians", "Multiple comparison test for ranks")
n=2:(nlevels(data$treatments)+1)
plot(response~treatments, data=data, col=n)
return(l)
}

f15 <-
function(data, alpha=0.05){
names(data)=c("treatments", "blocks", "response")
data=data.frame(treatments=factor(data$treatments),blocks=as.factor(data$blocks), response=data$response)
d = split(data, data[,2])
f1=function(data){rank(data$response, ties.method = c("average"))}
ran=lapply(d, f1)
da=data.frame(ran)
rownames(da)=d[[1]]$treatments
su=rowSums(da, na.rm=TRUE); sua=su^2; suu=sum(sua)
b=nlevels(data[,2])
t=nlevels(data[,1])
gl=t-1
cal=((12/(b*t*(t+1)))*(suu))-(3*b*(t+1))
pq=pchisq(cal,lower.tail=FALSE, df=gl) #chisq
f11=function(data){mean(data$response, na.rm=TRUE)}
ran2=lapply(d, f11)
so=data.frame(su, na=names(su))
so=so[order(so[,2], decreasing=FALSE),]
a1=aggregate(.~treatments, data, FUN=mean)
a2=aggregate(.~treatments, data, FUN=median)
da=data.frame(names=levels(a1$treatments), ra=so[,1], means=a1[,3],md=a2[,3])
ma=da[order(da[,2], decreasing=TRUE),]
j=ma[,1];j=as.character(j)
aux <- combn(j, 2)
w <- apply(aux, 2, paste, collapse = " - ")
jj=ma[,2]
auxj <- combn(jj, 2)
yi=auxj[1,]-auxj[2,]
dm=sqrt(b*t*(t+1)/6)
valores=yi/dm
prob=1-pnorm(valores)
prob1=prob*2
prob2=p.adjust(prob1, method="holm")
prob3=p.adjust(prob1, method="bonferroni")
prob4=p.adjust(prob1, method="fdr")
resp=data.frame(w,round(yi,4),round(prob1,4), round(prob2,4), round(prob3,4), round(prob4,4))
resp
colnames(resp)=c("pair", "contrast", "p(non adjusted)" ,"p.adj(Holm)", "p.adj(Bonferroni)", "p.adj(fdr)")
med=ft(resp, alpha=alpha)
colnames(med)=c("non adjusted", "adjust.Holm", "adjust.Bonferroni", "adjust.fdr" ) 
med=data.frame(treatment=rownames(med),rank=round(ma$ra,4), mean=round(ma$means,4), median=round(ma$md,4),med)
rownames(med)=NULL
med
pri=data.frame(round(c(cal,pq),4));colnames(pri)="Estimates"
rownames(pri)=c("Friedman chi-squared = ","p.value = ")
l=list(pri,med,resp) ; names(l)=c("Friedman Rank Sum Test", "Ranks, Means and Medians", "Multiple comparison test for ranks")
n=2:(nlevels(data$treatments)+1)
plot(response~treatments, data=data, col=n)
return(l)
}

   
    de1=c(1); de2=c(1,2);de3=c(1,2,3);de4=c(1,2,3,4);de5=c(1,2);de6=c(1,2,3)
    de7=c(1,2,3); de8=c(1,2); de9=c(1,2,3); de10=c(1,2,3);de11=c(1,2,3);de12=c(1,2,3);de13=c(1,2,3,4)
    de14=c(1);de15=c(1,2)
    de=list(de1,de2,de3,de4,de5,de6,de7,de8,de9,de10,de11,de12,de13,de14,de15)
    de=de[[design]]
    d=as.list(data)
    d1=d[de]
    d2=d[-de]
    f=function(h){data.frame(d1,d2[h])}
    h=length(d2)
    h=1:h
    l=lapply(h, f)
    l2=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15)
    fun=l2[[design]]
    li1=lapply(l, fun)
    names(li1)=names(d2)
    li=list(fun(data),li1)
    li=li[[list]]
    return(li)
}



