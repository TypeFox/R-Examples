pdf("nomographs.pdf")

fig2.1t3 <- read.table("TWRI3(A5)fig2_1t3.dat", header=TRUE)
fig2.2t3 <- read.table("TWRI3(A5)fig2_2t3.dat", header=TRUE)
fig2.3t3 <- read.table("TWRI3(A5)fig2_3t3.dat", header=TRUE)
fig2.vert <- read.table("TWRI3(A5)fig2_vert.dat", header=TRUE)
with(fig2.1t3, plot(h.over.P, C, type="l",
					xlim=c(0,5), ylim=c(3.0,5.0)))
with(fig2.2t3, lines(h.over.P, C))
with(fig2.3t3, lines(h.over.P, C))
with(fig2.vert, lines(h.over.P, C))
mtext("Discharge coefficients for full width, vertical, and
inclined sharp-crested rectangular weirs")
tmp <- lowess(fig2.1t3$C~fig2.1t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig2.2t3$C~fig2.2t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig2.3t3$C~fig2.3t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig2.vert$C~fig2.vert$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)



fig3.0.90 <- read.table("TWRI3(A5)fig3_0.90.dat", header=TRUE)
fig3.0.20 <- read.table("TWRI3(A5)fig3_0.20.dat", header=TRUE)
fig3.0.40 <- read.table("TWRI3(A5)fig3_0.40.dat", header=TRUE)
fig3.0.60 <- read.table("TWRI3(A5)fig3_0.60.dat", header=TRUE)
fig3.0.70 <- read.table("TWRI3(A5)fig3_0.70.dat", header=TRUE)
fig3.0.80 <- read.table("TWRI3(A5)fig3_0.80.dat", header=TRUE)
with(fig3.0.90, plot(h.over.P, kc, type="l",
					xlim=c(0,2.8), ylim=c(0.70,1.00)))
with(fig3.0.20, lines(h.over.P, kc))
with(fig3.0.40, lines(h.over.P, kc))
with(fig3.0.60, lines(h.over.P, kc))
with(fig3.0.70, lines(h.over.P, kc))
with(fig3.0.80, lines(h.over.P, kc))
mtext("Definition fo adjustment factor, Kc, for contracted sharp-crested weirs")
tmp <- lowess(fig3.0.20$kc~fig3.0.20$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig3.0.40$kc~fig3.0.40$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig3.0.60$kc~fig3.0.60$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig3.0.70$kc~fig3.0.70$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig3.0.80$kc~fig3.0.80$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig3.0.90$kc~fig3.0.90$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)




fig4.0.20 <- read.table("TWRI3(A5)fig4_0.20.dat", header=TRUE)
fig4.0.25 <- read.table("TWRI3(A5)fig4_0.25.dat", header=TRUE)
fig4.0.50 <- read.table("TWRI3(A5)fig4_0.50.dat", header=TRUE)
fig4.1.00 <- read.table("TWRI3(A5)fig4_1.00.dat", header=TRUE)
fig4.2.00 <- read.table("TWRI3(A5)fig4_2.00.dat", header=TRUE)
with(fig4.0.20, plot(ht.over.H, kt, type="l",
					xlim=c(0,1), ylim=c(0.3,1.1)))
with(fig4.0.25, lines(ht.over.H, kt))
with(fig4.0.50, lines(ht.over.H, kt))
with(fig4.1.00, lines(ht.over.H, kt))
with(fig4.2.00, lines(ht.over.H, kt))
mtext("Definition of adjustment factor, kt, for submerged sharp-crested rectangular weirs")
tmp <- lowess(fig4.0.20$kt~fig4.0.20$ht.over.H, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig4.0.25$kt~fig4.0.25$ht.over.H, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig4.0.50$kt~fig4.0.50$ht.over.H, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig4.1.00$kt~fig4.1.00$ht.over.H, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig4.2.00$kt~fig4.2.00$ht.over.H, f=1/10); lines(tmp$x, tmp$y, col=2)



fig6.1t3 <- read.table("TWRI3(A5)fig6_1t3.dat", header=TRUE)
fig6.2t3 <- read.table("TWRI3(A5)fig6_2t3.dat", header=TRUE)
fig6.3t3 <- read.table("TWRI3(A5)fig6_3t3.dat", header=TRUE)
fig6.vert <- read.table("TWRI3(A5)fig6_vert.dat", header=TRUE)
with(fig6.1t3, plot(h.over.P, h.over.L, type="l",
					xlim=c(0,5), ylim=c(1.2,3.0)))
with(fig6.2t3, lines(h.over.P, h.over.L))
with(fig6.3t3, lines(h.over.P, h.over.L))
with(fig6.vert, lines(h.over.P, h.over.L))
mtext("Minimum h/L ratio for which the lower nappe will
clear the downstream corner of a broad-crested weir")
tmp <- lowess(fig6.1t3$h.over.L~fig6.1t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig6.2t3$h.over.L~fig6.2t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig6.3t3$h.over.L~fig6.3t3$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig6.vert$h.over.L~fig6.vert$h.over.P, f=1/10); lines(tmp$x, tmp$y, col=2)



fig7.1p5t1 <- read.table("TWRI3(A5)fig7_1p5t1.dat", header=TRUE)
fig7.1t1 <- read.table("TWRI3(A5)fig7_1t1.dat", header=TRUE)
fig7.2t1 <- read.table("TWRI3(A5)fig7_2t1.dat", header=TRUE)
fig7.vert <- read.table("TWRI3(A5)fig7_vert.dat", header=TRUE)
with(fig7.1p5t1, plot(h.over.L, C, type="l",
                      xlim=c(0,2), ylim=c(2.4,3.8)))
with(fig7.1t1, lines(h.over.L, C))
with(fig7.2t1, lines(h.over.L, C))
with(fig7.vert, lines(h.over.L, C))
mtext("Coefficients of discharge for full width, broad-crested weirs
with downstream slope 1:1 and various upstream slopes")
tmp <- lowess(fig7.1p5t1$C~fig7.1p5t1$h.over.L, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig7.1t1$C~fig7.1t1$h.over.L, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig7.2t1$C~fig7.2t1$h.over.L, f=1/10); lines(tmp$x, tmp$y, col=2)
tmp <- lowess(fig7.vert$C~fig7.vert$h.over.L, f=1/10); lines(tmp$x, tmp$y, col=2)


weir.nomographs <- new.env(hash=TRUE)
assign("fig2",new.env(hash=TRUE), envir=weir.nomographs)
assign("fig3",new.env(hash=TRUE), envir=weir.nomographs)
assign("fig4",new.env(hash=TRUE), envir=weir.nomographs)
assign("fig6",new.env(hash=TRUE), envir=weir.nomographs)
assign("fig7",new.env(hash=TRUE), envir=weir.nomographs)


tmp <- get("fig2", envir=weir.nomographs)
assign("0.0000", fig2.vert, envir=tmp)
assign(as.character(sprintf("%.4f",1/3)), fig2.1t3, envir=tmp)
assign(as.character(sprintf("%.4f",2/3)), fig2.2t3, envir=tmp)
assign(as.character(sprintf("%.4f",3/3)), fig2.3t3, envir=tmp)




tmp <- get("fig3", envir=weir.nomographs)
assign("0.20", fig3.0.20, envir=tmp)
assign("0.40", fig3.0.40, envir=tmp)
assign("0.60", fig3.0.60, envir=tmp)
assign("0.70", fig3.0.70, envir=tmp)
assign("0.80", fig3.0.80, envir=tmp)
assign("0.90", fig3.0.90, envir=tmp)




tmp <- get("fig4", envir=weir.nomographs)
assign("0.20", fig4.0.20, envir=tmp)
assign("0.25", fig4.0.25, envir=tmp)
assign("0.50", fig4.0.50, envir=tmp)
assign("1.00", fig4.1.00, envir=tmp)
assign("2.00", fig4.2.00, envir=tmp)




tmp <- get("fig6", envir=weir.nomographs)
assign(as.character(sprintf("%.4f",1/3)), fig6.1t3, envir=tmp)
assign(as.character(sprintf("%.4f",2/3)), fig6.2t3, envir=tmp)
assign(as.character(sprintf("%.4f",3/3)), fig6.3t3, envir=tmp)
assign("0.0000", fig6.vert, envir=tmp)




tmp <- get("fig7", envir=weir.nomographs)
assign(as.character(sprintf("%.4f",0.5/1)), fig7.1p5t1, envir=tmp)
assign(as.character(sprintf("%.4f",1/1)), fig7.1t1, envir=tmp)
assign(as.character(sprintf("%.4f",2/1)), fig7.2t1, envir=tmp)
assign("0.0000", fig7.vert, envir=tmp)


X <- c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14);
Y <- c(1, 1.01, 1.03, 1.04, 1.05, 1.06, 1.08, 1.09);
tmp <- data.frame(R.over.h=X, kr=Y);
assign("broadcrest.roundingtable", tmp, envir=weir.nomographs)


assign("broadcrest.downstreamtable", new.env(hash=TRUE), envir=weir.nomographs)
tmp <- get("broadcrest.downstreamtable", envir=weir.nomographs)
assign("h.over.L", c(0.1, 0.4, 1.0, 2.0), envir=tmp);
assign(as.character(sprintf("%.4f",2/1)), c(1,1,.98,.98), envir=tmp);
assign(as.character(sprintf("%.4f",3/1)), c(1,1,.96,.94), envir=tmp);
assign(as.character(sprintf("%.4f",4/1)), c(1,1,.95,.91), envir=tmp);
assign(as.character(sprintf("%.4f",5/1)), c(1,1,.94,.90), envir=tmp);
slopes <- c(as.character(sprintf("%.4f",2/1)),
            as.character(sprintf("%.4f",3/1)),
            as.character(sprintf("%.4f",4/1)),
            as.character(sprintf("%.4f",5/1)));
assign("slopes", slopes, envir=tmp);




fig23.gravel <-
as.data.frame(matrix(c(
0,  2.5,
0.0459917998543,  2.52068965517,
0.103577674013,  2.54062066747,
0.153419376942,  2.56131254090,
0.203290656591,  2.58047381906,
0.256997050955,  2.60116791074,
0.295318621380,  2.61802664123,
0.356739610253,  2.63949046705,
0.406625678262,  2.65788644758,
0.479596376359,  2.68165282106,
0.544867268306,  2.70388416278,
0.602512295905,  2.72075398454,
0.687077068306,  2.74452701279,
0.794682105715,  2.7759663269,
0.894646490414,  2.80280941871,
0.963870803957,  2.82045119287,
1.02923042606,  2.83809074877,
1.10239337284,  2.85190825301,
1.17555631962,  2.86572575726,
1.25640428419,  2.88184359091,
1.35271101449,  2.89797029758,
1.44903253316,  2.91333170661,
1.52223984502,  2.92485331795,
1.59546194524,  2.93560963166,
1.68412802283,  2.94714011601,
1.75355937342,  2.95406772329,
1.82681105036,  2.96329344173,
1.96182384845,  2.9763811404,
2.03513467883,  2.9825456683,
2.10840114413,  2.99100608911,
2.19330604881,  2.99717727177,
2.30136952538,  3.00489235922,
2.40168883072,  3.0133683078,
2.50206728950,  3.01878306585,
2.60243095992,  3.02496312153,
2.69123013276,  3.02960592717,
2.79163816826,  3.03349008995,
2.89204620377,  3.03737425273,
2.98860433619,  3.04049089961,
3.09286227477,  3.04514257828,
3.28215821327,  3.04907776089,
3.44832515660,  3.04993844345,
3.57970987702,  3.05077916172,
3.73428274605,  3.05163318952,
3.85799723704,  3.04940828074,
3.95073504311,  3.05022681647,
4.0,  3.05178181254), ncol=2, byrow=TRUE));
names(fig23.gravel) <- c("H", "C");

fig23.paved <-
as.data.frame(matrix(c(
0,  2.84132496312,
0.0158974870879,  2.87807034083,
0.0425658297687,  2.8979836071,
0.0731136522443,  2.91713379399,
0.115225972304,  2.93782123091,
0.161232560518,  2.95698029081,
0.203463187458,  2.97154534666,
0.257302677063,  2.98535175963,
0.315051223182,  2.99686449795,
0.376694037456,  3.00684885926,
0.442245908246,  3.01453954592,
0.496248069811,  3.01992768492,
0.550265019736,  3.02455052628,
0.612011352531,  3.02917780415,
0.666087455897,  3.03073945497,
0.739487016436,  3.03231219707,
0.812842211895,  3.03618083207,
0.913324189197,  3.03623850667,
1.006047206910,  3.03782234004,
1.098799801340,  3.03787557813,
1.195417087210,  3.03793103448,
1.404095636320,  3.03881611783,
1.605059590930,  3.03893146704,
1.883287797510,  3.04062177660,
2.188583632490,  3.04156231630,
2.474556010300,  3.04249176473,
2.609775845430,  3.04486529653,
2.799160514090,  3.04420869334,
3.015524081000,  3.04739410610,
3.185555715760,  3.04825700691,
3.409693030620,  3.04915096328,
3.575845185590,  3.05077694347,
3.718838768670,  3.05085901887,
3.850223489090,  3.05169973714,
4.0,  3.05178181254), ncol=2, byrow=TRUE));
names(fig23.paved) <- c("H", "C");

assign("fig23.gravel", new.env(hash=TRUE), envir=weir.nomographs)
assign("fig23.paved", new.env(hash=TRUE), envir=weir.nomographs)
tmp <- get("fig23.gravel", envir=weir.nomographs)
assign("H", fig23.gravel$H, envir=tmp);
assign("C", fig23.gravel$C, envir=tmp);
tmp <- get("fig23.paved", envir=weir.nomographs)
assign("H", fig23.paved$H, envir=tmp);
assign("C", fig23.paved$C, envir=tmp);


############################################# END TWRI

############################################# BEGIN BOS

Bos1989C  <- read.table("Bos1989C.dat",  header=TRUE)
Bos1989Cv <- read.table("Bos1989Cv.dat", header=TRUE)
with(Bos1989C, plot(h.over.L, C, type="l"))
mtext("Bos(1989) C")
with(Bos1989Cv, plot(CAstar.over.A, Cv, type="l"))
mtext("Bos(1989) CAstar.over.A")
assign("Bos1989C",  Bos1989C,  envir=weir.nomographs)
assign("Bos1989Cv", Bos1989Cv, envir=weir.nomographs)



dev.off()

.weir.nomographs <- weir.nomographs;
save(.weir.nomographs, file="WEIRSnomographs.RData");


