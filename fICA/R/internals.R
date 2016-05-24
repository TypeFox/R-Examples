mat.sqrt <- function(A)
{
 eig <- eigen(A, symmetric=TRUE)
 eig$vectors%*%(diag(eig$values^(1/2)))%*%t(eig$vectors)
}

mat.norm <- function(A)
{
 sqrt(sum(A^2))
}

G_pow3 <- function(x){.Call("Gpow3",as.matrix(x),PACKAGE="fICA")$gx-0.75}
G_tanh <- function(x){log(cosh(x))-0.3745672}
G_gaus <- function(x){-exp(-x^2/2)+0.7071068}
G_rt0.6 <- function(x){.Call("Grn",as.matrix(x),0.6,PACKAGE="fICA")$gx-0.07783844}
G_lt0.6 <- function(x){.Call("Gln",as.matrix(x),0.6,PACKAGE="fICA")$gx+0.07783844}

G_bt <- function(x){.Call("Gbn",as.matrix(x),0,PACKAGE="fICA")$gx-0.531923}
G_bt0.2 <- function(x){.Call("Gbn",as.matrix(x),0.2,PACKAGE="fICA")$gx-0.361278}
G_bt0.4 <- function(x){.Call("Gbn",as.matrix(x),0.4,PACKAGE="fICA")$gx-0.2399444}
G_bt0.6 <- function(x){.Call("Gbn",as.matrix(x),0.6,PACKAGE="fICA")$gx-0.1556769}
G_bt0.8 <- function(x){.Call("Gbn",as.matrix(x),0.8,PACKAGE="fICA")$gx-0.09857521}
G_bt1.0 <- function(x){.Call("Gbn",as.matrix(x),1,PACKAGE="fICA")$gx-0.06086077}
G_bt1.2 <- function(x){.Call("Gbn",as.matrix(x),1.2,PACKAGE="fICA")$gx-0.03660588}
G_bt1.4 <- function(x){.Call("Gbn",as.matrix(x),1.4,PACKAGE="fICA")$gx-0.02143076}
G_bt1.6 <- function(x){.Call("Gbn",as.matrix(x),1.6,PACKAGE="fICA")$gx-0.012203}

Gf <- c(G_pow3,G_tanh,G_gaus,G_lt0.6,G_rt0.6,G_bt,G_bt0.2,G_bt0.4,G_bt0.6,G_bt0.8,G_bt1.0, G_bt1.2,G_bt1.4,G_bt1.6)


g_pow3 <- function(x){.Call("gpow3",as.matrix(x),PACKAGE="fICA")$gx}
g_tanh <- function(x){tanh(x)}
g_gaus <- function(x){x*exp(-x^2/2)}

g_rt0.6 <- function(x){.Call("grn",as.matrix(x),0.6,PACKAGE="fICA")$gx}
g_lt0.6 <- function(x){.Call("gln",as.matrix(x),0.6,PACKAGE="fICA")$gx}

g_bt <- function(x){.Call("gbn",as.matrix(x),0,PACKAGE="fICA")$gx}
g_bt0.2 <- function(x){.Call("gbn",as.matrix(x),0.2,PACKAGE="fICA")$gx}
g_bt0.4 <- function(x){.Call("gbn",as.matrix(x),0.4,PACKAGE="fICA")$gx}
g_bt0.6 <- function(x){.Call("gbn",as.matrix(x),0.6,PACKAGE="fICA")$gx}
g_bt0.8 <- function(x){.Call("gbn",as.matrix(x),0.8,PACKAGE="fICA")$gx}
g_bt1.0 <- function(x){.Call("gbn",as.matrix(x),1,PACKAGE="fICA")$gx}
g_bt1.2 <- function(x){.Call("gbn",as.matrix(x),1.2,PACKAGE="fICA")$gx}
g_bt1.4 <- function(x){.Call("gbn",as.matrix(x),1.4,PACKAGE="fICA")$gx}
g_bt1.6 <- function(x){.Call("gbn",as.matrix(x),1.6,PACKAGE="fICA")$gx}

gf <- c(g_pow3,g_tanh,g_gaus,g_lt0.6,g_rt0.6,g_bt,g_bt0.2,g_bt0.4,g_bt0.6,g_bt0.8, g_bt1.0,g_bt1.2,g_bt1.4,g_bt1.6)



dg_pow3 <- function(x){.Call("dgpow3",as.matrix(x),PACKAGE="fICA")$gx}
dg_tanh <- function(x){1-tanh(x)^2}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}

dg_rt0.6 <- function(x){.Call("dgrn",as.matrix(x),0.6,PACKAGE="fICA")$gx}
dg_lt0.6 <- function(x){.Call("dgln",as.matrix(x),0.6,PACKAGE="fICA")$gx}

dg_bt <- function(x){.Call("dgbn",as.matrix(x),0,PACKAGE="fICA")$gx}
dg_bt0.2 <- function(x){.Call("dgbn",as.matrix(x),0.2,PACKAGE="fICA")$gx}
dg_bt0.4 <- function(x){.Call("dgbn",as.matrix(x),0.4,PACKAGE="fICA")$gx}
dg_bt0.6 <- function(x){.Call("dgbn",as.matrix(x),0.6,PACKAGE="fICA")$gx}
dg_bt0.8 <- function(x){.Call("dgbn",as.matrix(x),0.8,PACKAGE="fICA")$gx}
dg_bt1.0 <- function(x){.Call("dgbn",as.matrix(x),1,PACKAGE="fICA")$gx}
dg_bt1.2 <- function(x){.Call("dgbn",as.matrix(x),1.2,PACKAGE="fICA")$gx}
dg_bt1.4 <- function(x){.Call("dgbn",as.matrix(x),1.4,PACKAGE="fICA")$gx}
dg_bt1.6 <- function(x){.Call("dgbn",as.matrix(x),1.6,PACKAGE="fICA")$gx}


dgf <- c(dg_pow3,dg_tanh,dg_gaus,dg_lt0.6,dg_rt0.6,dg_bt,dg_bt0.2,dg_bt0.4,dg_bt0.6, dg_bt0.8,dg_bt1.0,dg_bt1.2,dg_bt1.4,dg_bt1.6)

gnames <- c("pow3","tanh","gaus","lt0.6","rt0.6","bt","bt0.2","bt0.4",
"bt0.6","bt0.8","bt1.0","bt1.2","bt1.4","bt1.6")

