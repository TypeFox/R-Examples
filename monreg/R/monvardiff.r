# Erläuterung:
# ============
# x,y       : Stichprobe
# a,b       : Grenzen des Definitionsbereichs
# N         : Länge des Vektors oder Vektor des unbedingten Schätzers m^dach
# t         : Länge des Vektors oder Vektor der Auswertungsstellen des monotonen Schätzers m^dachi_I (invers)
# hr,hd      : Bandbreiten
# Kr,Kd      : Kerne
# degree    : Wahl des unbedingten Schätzers (0 = Nadaraya Watson, 1 = Local Linear)
# inverse   : Flag ob originäre Regressionsfunktion (inverse = 0) oder inverse (inverse = 1) geschätzt werden soll


monvardiff <- function(x,y,a=min(x),b=max(x),N=length(x),t=length(x),r=2,hr,Kr="epanech",hd,Kd="epanech",degree=1,inverse=0,monotonie="isoton"){

    xs <- (x-a)/(b-a)

    if(length(N) == 1 && length(degree) == 1)
    z <- seq(min(xs),max(xs),length=N)
    if(length(N) > 1)
    z <- N

    if(length(t) == 1){
        t = rep(0,length=t)
        lt = as.integer(length(t))
        tflag = as.integer(1)
    } else{
        lt = as.integer(length(t))
        tflag = as.integer(0)
    }

    Kern0 = switch(as.character(Kr), epanech=1, rectangle=2, biweight=3, triweight=4, triangle=5, cosine=6)
    Kern0 = as.integer(Kern0)
    Kern1 = switch(as.character(Kd), epanech=1, rectangle=2, biweight=3, triweight=4, triangle=5, cosine=6)
    Kern1 = as.integer(Kern1)

    lx <- as.integer(length(x))
    lz <- as.integer(length(z))
    ldeg <- as.integer(length(degree))
    erg <- seq(1,1,length=length(t))
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(xs) <- "double"
    storage.mode(z) <- "double"
    storage.mode(t) <- "double"
    storage.mode(tflag) <- "integer"
    storage.mode(r)<-"integer"
    storage.mode(hr) <-"double"
    storage.mode(hd) <- "double"
    storage.mode(erg) <- "double"
    storage.mode(degree) <- "double"
    storage.mode(inverse) <- "integer"

    if(length(x) != length(y)) print("x and y must have the same length")
    else{
    if (Kern0 == 0 || Kern1 == 0 ) print("no valid kernel selected")
        else{
        if (inverse!=0 & inverse!=1) print ("inverse must have value 0 or 1")
        else{
        if(length(degree) > 1 && length(degree) != length(N)) print("lengths of N and sdegree differ")
        else{
            if(monotonie == "antiton") temp <- .C("sdach2_a_inv",xs,y,z,t,lx,lz,lt,r,tflag,hr,hd,Kern0,Kern1,degree,ldeg,inverse,erg)
            else temp <- .C("sdach2_i_inv",xs,y,z,t,lx,lz,lt,r,tflag,hr,hd,Kern0,Kern1,degree,ldeg,inverse,erg)
            }
        }
    }
    }

    K = switch(as.integer(temp[[12]]),"epanech","rectangle","biweight","triweight","triangle","cosine")
    Kd = switch(as.integer(temp[[13]]),"epanech","rectangle","biweight","triweight","triangle","cosine")

    if(inverse == 0) ergebnis <- list(xs=temp[[1]],y=temp[[2]],z=temp[[3]],t=temp[[4]]*(b-a)+a,length.x=temp[[5]],length.z=temp[[6]],length.t=temp[[7]],r=temp[[8]],hr=temp[[10]],hd=temp[[11]],Kr=Kr,Kd=Kd,degree=temp[[14]],ldeg.vektor=temp[[15]],inverse=temp[[16]],estimation=temp[[17]])
    else ergebnis <- list(xs=temp[[1]],y=temp[[2]],z=temp[[3]],t=temp[[4]],length.x=temp[[5]],length.z=temp[[6]],length.t=temp[[7]],r=temp[[8]],hr=temp[[10]],hd=temp[[11]],Kr=Kr,Kd=Kd,degree=temp[[14]],ldeg.vektor=temp[[15]],inverse=temp[[16]],estimation=temp[[17]]*(b-a)+a)
}

