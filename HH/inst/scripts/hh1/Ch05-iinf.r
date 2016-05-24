library(HH)

#### iinf/code/vocab.s
data(vocab)
stem(vocab$score)
t.test(vocab$score, mu=10)


#### iinf/code/teachers.s
data(teachers)
teachers[c(1,2,32),]
if.R(s=
     t(bwplot(teachers))
     ,r={
       tmp <- data.frame(x=unlist(teachers),
                         language=rep(names(teachers), c(32,32)))
       bwplot(x ~ language, data=tmp)
   }
)
if.R(r=parallelplot( ~ teachers),
     s=parallel( ~ teachers))

teachers.diff <- teachers$English - teachers$Greek
teachers.diff
if.R(r=
     parallelplot( ~ teachers,
                  lty=c(1,3)[1+(teachers.diff > 0)],  ## control line type
                  scales=list(cex=1.5))               ## control axis labels
     , s=
     parallel( ~ teachers,
              lty=c(1,3)[1+(teachers.diff > 0)],  ## control line type
              scales=list(cex=1.5))               ## control axis labels
     )
stem(teachers.diff)                  ## first try
if.R(s=
     stem(teachers.diff, nl=5, scale=-1)  ## better scaling
     ,r=
     stem(teachers.diff, scale=2)  ## different scaling
     )
stem(sqrt(teachers.diff + 17))
t.test(sqrt(teachers.diff + 17), mu=sqrt(17))
t.test(teachers.diff)


#### iinf/code/teachers2.s
## teachers is defined in teachers.s (above)

tmp <- data.frame(errors=as.vector(unlist(teachers)),
                  sent.lang=rep(names(teachers), c(32,32)),
                  sentence=rep(1:32,2))
tmp[c(1,2,32,33,34,64),]


new.order <- order(sign(teachers$English - teachers$Greek),
                   apply(teachers, 1, min))
if.R(s=
     ordered(tmp$sentence) <- new.order
     ,r=
     tmp$sentence <- ordered(tmp$sentence, levels=new.order)
     )

tpg <- trellis.par.get("superpose.symbol")
tpg$pch[3] <- 32 ## ascii for blank

tmp.trellis <-
dotplot(sentence ~ errors, data=tmp,
        panel=function(x,y,...) {
          superpose.symbol <- trellis.par.get("superpose.symbol")
          dot.line <- trellis.par.get("dot.line")
          panel.abline(h = y, lwd = dot.line$lwd, lty = dot.line$lty,
                 col = dot.line$col)
          panel.superpose(x,y,...)
          panel.abline(h=21.5, lty=2)
          if.R(r={
            panel.segments(-12,21.5,44,21.5, lty=3, xpd=TRUE)
          },
               s={
                 segments(-12,21.5,44,21.5, lty=3, xpd=TRUE)
                 mtext(side=2, at=c(28.5,13.5), line=6,
                       c("Greek\nsmaller\nerror\ncount",
                         "English\nsmaller\nerror\ncount"))
               })
        },
        groups=sent.lang,
        key=list(text=
          list(c("",
            "English for", "Greek Speakers","",
            "Greek for", "English Speakers","")),
          points=Rows(tpg, c(3,1,3,3,2,3,3)),
          space="right",
          border=TRUE),
        ylab=if.R(
          r=list(
            rev(c("Greek\nsmaller\nerror\ncount",
                  "English\nsmaller\nerror\ncount")),
            rot=0),
          s=NULL)
)
print(position=c(.15,0,.85,1),
      tmp.trellis
      )
## export.eps(hh("iinf/figure/teachers-dot.eps"))

## based on trellis library example.overplot


#### iinf/code/ex.qqnorm.s

x1 <- rnorm(100,0,1)
x3 <- runif(100,-3,3)
x2 <- x1/x3
x4 <- x1*exp(-abs(x1))
x5 <- x1^2
x6 <- -x5

dist.names <- c("normal",
                "heavy-tailed",
                "uniform",
                "thin-tailed",
                "positively skewed",
                "negatively skewed")

qqnorm.dat <-
data.frame(dist=factor(rep(dist.names, rep(100,6)), levels=dist.names),
           x=c(sort(x1), sort(x2), sort(x3), sort(x4), sort(x5), sort(x6)),
           qn=rep(qnorm(ppoints(length(x1))),6)
)

xyplot(x ~ qn | dist, data=qqnorm.dat,
       par.strip.text=list(cex=1.5),
       layout=c(2,3), as.table=TRUE,
       scales=list(cex=1, y=list(relation="free"), alternating=FALSE),
       between=list(x=1,y=1),
       xlab="Quantiles of Standard Normal",
       ylab="distribution")
## export.eps(hh("iinf/figure/iinf.f.qqnorm.ps"))


#### iinf/code/ex.tquant.s
old.par <- par(cex=1.2, pch=16, pty="s",
               mar=par("mar")+if.R(s=c(0,1,4,0), r=c(0,7,4,3)))
par(plt=par()$plt+c(-.21,-.21,0,0))

## d1 <- rt(100,5)
## ## the d1 we used is in file ex.tquant.d1.s
## source(hh.old("scripts/ex.tquant.d1.s")) ## dump of the d1 object
## ## ex.tquant.d1.s  lives in the same script directory as Ch05-iinf.r
## ## source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/ex.tquant.d1.s")

## this is the d1 used in file iind.f1.ps.gz
d1 <-
c(0.193109809054234, -0.526543795724567, -0.296739977034572, 0.710778693915572,
	-0.0629813594467444, -0.736529190800633, -0.699735812697692,
	0.171675764341657, -0.0641886971261035, 9.42364991049266e-005,
	-0.0343373973117536, 1.60189771844038, 1.48622846237989,
	2.12724167273368, 1.46847478555014, 0.413877241989935, 2.30799285820845,
	-1.27487619533415, 0.82769703851613, 0.241492057718752,
	-1.04986928969301, -1.78187014856777, 1.37792650954825, 1.01605624380149,
	-1.51233891092346, 0.357773786351872, -0.36562569561658,
	0.191330659513942, 0.460337380085628, -1.81233854096677,
	-1.99487322592861, -1.21726353848997, -0.353299333611969,
	-0.576351526287996, -1.43023555053075, -0.859289205958409,
	1.10032220523717, 0.177163843112181, 0.857271241430201,
	-0.926170364073415, 0.179639843675247, 0.229129068673923,
	-0.814701472508782, 0.65542043438568, -1.10451979096874,
	-0.522508906468807, -0.906101410767418, 0.448010406027382,
	0.0886217586096557, 0.195298848183735, -0.982665528368598,
	-0.213506984965018, 0.244095808280918, 1.5510214000569,
	0.358774799743031, 0.56109745560833, 2.22559931541281, -1.39999413496839,
	0.605403184840874, 1.55420055828709, -0.130139565877701,
	-0.157874208968318, -0.580899935654221, 0.114980781465297,
	0.0658414724467252, 0.802104022917183, -0.532722466312124,
	-6.27799675313213, 0.654897330007934, -0.279461025710469,
	2.11694277017957, -0.118076779874461, 0.139101531208278,
	-1.69225497386311, 1.81941339507338, 0.0909917468413372,
	0.624684066432928, 1.66675829939121, 2.50983875236152, -0.12387663806769,
	-0.237335394129467, -1.18366243931167, 0.854060836324612,
	0.288833945493047, 0.479631374808358, -0.36800123495364,
	-0.640393391388562, -1.36056793817533, -0.875377991107211,
	0.900716863288967, 0.77760323740068, -3.16539175731614,
	-0.760007741408089, -2.4064647620079, -0.816193659240182,
	-1.34868267174933, -0.975141669495992, -0.364302703839452,
	-0.439874237226035, -0.831565509836057)


q3 <- qt(ppoints(d1), 3)
q5 <- qt(ppoints(d1), 5)
q7 <- qt(ppoints(d1), 7)
qn <- qnorm(ppoints(d1))


qqplot(d1, q3,
       type="l", lty=1,
       main="Selected Q-Q Plots",
       xlab="sample from t with 5 df",
       ylab="comparison values",
       xlim=c(-6,6), ylim=c(-6,6))

lines(sort(d1), q5, lty=2)

lines(sort(d1), q7, lty=3)

lines(sort(d1), qn, lty=4)

abline(a=0, b=1, lty=5)

if.R(r=
     old.xpd <- par(xpd=TRUE)
     ,s=
     {}
   )
legend(7, 2,
       c("quantiles of t, 3 df",
         "quantiles of t, 5 df",
         "quantiles of t, 7 df",
         "quantiles of normal",
         "slope=1"),
       lty=c(1,2,3,4,5))
if.R(r=
     par(old.xpd)
     ,s=
     {}
     )
par(old.par)
## export.eps(hh("iinf/figure/iinf.f1.ps"))


#### iinf/code/ex.tquant.d1.s
## read above


#### iinf/code/ex.ttail.s
## continues with data from ex.tquant.s
ppd <- ppoints(d1)
    d3 <- dt(q3, 3)
    d5 <- dt(q5, 5)
    d7 <- dt(q7, 7)
    dn <- dnorm(qn)

    plot( x=q3, y=d3, lty=1, type="l", ylim=c(0,.4),
         main="density functions",
         xlab="quantile q", ylab="density")
    lines(x=q5, y=d5, lty=2)
    lines(x=q7, y=d7, lty=3)
    lines(x=qn, y=dn, lty=4)

    legend(-5.5, .3,
           c("t, 3 df",
             "t, 5 df",
             "t, 7 df",
             "normal"),
           lty=c(1,2,3,4))

    abline(h=0)
## export.eps(hh("iinf/figure/iinf.f2.ps"))


#### iinf/code/goffit.s
par(mfcol=c(2,3))
old.par <- par(mar=par("mar")+c(0,1,0,0))

## t9 <- rt(100,9)
## ## the t9 we used is in file goffit-t9-t3-t5.s
## source(hh.old("scripts/goffit-t9-t3-t5.s")) ## dump of the t9 object
## ## goffit-t9-t3-t5.s lives in the same script directory as Ch05-iinf.r
## ## source("c:/HOME/rmh/HH-R.package/HH/inst/scripts/goffit-t9-t3-t5.s")
## ## to reproduce the data in hh("iinf/figure/iinf.f3.ps.gz")

## dump(c("t9","t3","t5"), fileout=hh("iinf/code/goffit-t9-t3-t5.s"))
"t9" <-
c(-0.85967772221988004, -0.57192563874817859, 0.22448081393808378,
	0.9103287131836052, 0.13956168888555937, -0.52813077887697146,
	0.3433207222357153, -0.12410310794500161, 0.077707893120165641,
	1.4244997372322139, -0.27453471933480078, 0.65600728413868503,
	-1.7523778488185169, 0.463447546394335, 1.217357370681688,
	-0.54666775559422864, -0.63167259759492866, 1.3334307196712061,
	-0.13016624422995282, 0.65151727948733396, 0.463824551191167,
	-2.4502795055155056, 0.56777445717250685, -0.76989008060760933,
	-0.17753676844646382, -0.4930398790337352, 0.02474359708132072,
	-0.36939076435974677, -0.2117184183074467, -2.2272560043825478,
	-0.16624943700175271, -1.3904860745365339, 0.13966322523135341,
	-0.81593594015069093, -1.6218115861330897, -0.0080633824465966997,
	1.6650626670996551, 1.4032508603986984, -1.4013814582270485,
	-1.312748051317548, -0.090368515603434826, 0.37468319734206923,
	-1.145916733973221, -0.025756950401791798, 0.046356151066127012,
	-0.77557111639345266
	, 0.14974347820626671, -0.40712895827746964, 0.59859494249520162,
	0.25815291138896385, -0.96058425229139466, -1.1853584990942738,
	0.93224019681193016, -1.2420820468170841, -0.073626703106132174,
	0.67574163980361224, -1.0315527642599782, 0.66732149270820007,
	-0.42773165242262562, 0.73001207283437097, -0.68343711866788892,
	-0.98926402485740439
	, 0.5387171519416577, 0.41534388902569852, -1.0772539618913597,
	1.243178978993343, -0.20812170333719371, -1.2282395079434336,
	0.55218500580155638, 1.3732290883022595, 1.3336601016658585,
	-0.00078745556507992074, -0.24980086929099499, 0.41930574502793866,
	0.50459775498651582, 1.135592968330712, -1.0310937460495082,
	-1.2854182264736749, -1.6075694815150505, 0.46314413596402526,
	-1.600096912369231, 1.0956752771848868, 0.95719572945661058,
	-0.30235104369106264, 0.53049014821666141, -1.9111258560606157,
	-0.73815821835477657, -0.1579025234470007, 1.1825072137845929,
	-1.7452538498485046, -0.2111533814457561, -0.52588811634389854,
	1.4383051706926433, 0.84850085961802335, 0.33421002583981113,
	0.22709755974686791, -0.20674276948560583, 3.7334376945022125,
	1.9611678294594663, -0.37692692697311869)

"t3" <-
c(0.57599255621037659, -0.55630535060784692, -1.3022604847252943,
	0.61480652994671736, -0.67262434171894037, -3.0970497704103175,
	0.58093803148922918, -0.5944239977023118, -13.322589657205517,
	-0.75577238396596125, -0.50916646155825973, -0.41043304164324523,
	-0.60066023268342539
	, -4.3653426599370659, -0.51444065685286611, -1.4177220616698756,
	0.18293426956092879, -1.7300101589489107, 0.34178247524087024,
	-1.0551366727128002, -1.3181065525462337, -0.25328414596899573,
	0.49009289178610299, 1.5237714226179329, -1.4877123363761851,
	-0.43954137801569609, -0.73381747906824213, -1.3212774403296652,
	0.93315785956004504, -3.0871119310137249, 0.75499775539323899,
	-0.25212813911678233, 0.40526634681640916, 0.20545841981898308,
	-0.20902930724953708, -0.28134826619538583, 2.1153618146708428,
	1.189770661663667, 0.072637444427467834, -0.51774844329499647,
	-0.75211336609129675, 0.40184584960270359, -0.20981632178267054,
	-1.6444752292240008
	, 0.28975445426996244, -1.6541167157490826, 0.16783772313014933,
	0.57211175188103458, 1.1290603706079052, 0.50009077961525406,
	0.54383472506476527, -0.13953381529677111, -0.18967367068109497,
	-0.46863438493568232
	, -1.4802349130080592, -0.32007469623369361, -0.23431027063261689,
	-0.018705522516544525, 0.23814109510942971, -0.066649708592885007,
	-1.3890928446020436, -3.7403771686736276, 0.52215965379912166,
	-0.99192239667475035, 1.2570880691973261, 0.42419064756741093,
	0.24336077649528215, -0.090207541228995533, -0.1140770811584083,
	-2.3905162859816826
	, 1.3967233251728994, 1.4139374598188128, 1.4841301305798686,
	0.65369629751248648, -0.49753024009825186, 0.1909927751762876,
	2.3409508270524695, 0.4100769354407085, -0.0018660057280440445,
	0.19978118501088679, -1.1928377554298497, -1.5674006895694292,
	1.1825384915585826, 0.028265746928356438, -1.2210923859691278,
	0.82810283190372258, 1.2868715083798068, 0.29821282684944594,
	-0.16049779311640869, -0.026288393653549534, 0.54939273732913152,
	-3.164593723611544
	, -0.16654064090003901, 0.3208301705364941, -0.0089215623454682782,
	-2.7921913463967454, 0.99338318090918887, -0.13546227870728586,
	-0.79393354931631466, 1.8522002951187633)

"t5" <-
c(-2.1371965377220454, 0.25206598032195426, -3.6418647380419249,
	-2.2866067103111027, -0.52913480318952244, -0.91227663260401126,
	0.93754923126854162, 2.1227543290044788, -0.12458941304792399,
	-1.4799526595446653, 2.4342018019633862, -0.78288289775408959,
	2.0727445198888046, 1.9875546492245184, -0.11090554362299379,
	-2.3863435976884277, 0.021511933160897799, 0.26472400133103707,
	0.18488291227110537, 3.3672749181515171, 0.3319315308185517,
	-0.69868915210150084, 0.38967659609368088, -1.4099499630331633,
	1.3007936890606187, -0.59734403641677614, -0.67338750978424056,
	1.8920588841870651, -1.3298425544921169, 0.80485321544944066,
	-0.53534736185235698, -1.0658262234403584, 3.4814935933211633,
	0.12367331653115093, 0.45327120676573018, 0.82944949084145714,
	-0.42523878913929741, 2.2507747111456804, -0.018286227810069555,
	3.1729395467633825, -1.20793386812646, -0.96238163505895913,
	-4.1052058008126107, 0.78235565479096658, 1.2003870560549568,
	0.26751924607976729, 1.1724315406435233, -0.49734535427049686,
	0.39913258717321648, -2.2067328396520107, -0.32343537925015736,
	1.8141637898667986, 0.33898868010510014, -1.0207813793098146,
	0.64587901696870686, -0.97277369503049327, 0.26632850587123325,
	1.9898463355904088, -0.89815146433051296, 0.40723338278218002,
	0.46354350640026559, -0.36368467032831137, -0.26109002028189504,
	-0.93279034042074704
	, -0.53257375892533221, 1.363144370338716, -1.9407622093616841,
	1.1795314154730565, 0.58873054176395256, 0.026616725793430038,
	3.3093101302579795, 4.3426732225285196, -2.0691881650684847,
	0.65069867614812282, -0.76046072287552446, -1.5188184866472569,
	0.84812588367117303, -0.030685977590295361, -0.58813032367258689,
	0.65884204975922112, -0.33623944689042418, 3.9335297719896891,
	0.58821266747465484, -0.70200329637462144, 0.17108815538899855,
	1.3419841516002873, 2.5090903406077021, -0.6317267339390582,
	-0.060160998450137337, -0.56280627836437458, 2.050657400734111,
	-1.0894470875620508, -1.5936176820376864, 0.74911665591102417,
	-1.1867292962345946, 0.86373098144543237, -0.74844491172293481,
	1.233453177774954, -0.88427213552444683, -0.5085656557929864)



t9s <- sort(t9)
p.fractions <- ppoints(t9s)
length(p.fractions)
p.fractions[1:5]
p.fractions[96:100]



## one sample when H0 is true
if.R(s=
     ks.gof(t9, distribution="t", df=9)
     ,r=
     ks.test(t9, y="pt", df=9)
     )
t9.quantiles <- qt(p.fractions, 9)

plot(x=range(t9, t9.quantiles), y=c(0,1), type="n", cex=1)
title(main="KS for H0: t9 ~ t with 9 df", cex=.9)
lines(x=t9.quantiles, y=p.fractions)    #reference distribution
points(x=t9s, y=p.fractions, pch=16, cex=.3)

vert.obs.99 <- pt(t9s, 9)               #additional reference distribution
lines(x=t9s, y=vert.obs.99, lty=2)      #to be used for finding ks statistic

if.R(s=
     segments(x1=t9s, y1=p.fractions, x2=t9s, y2=vert.obs.99)
     ,r=
     segments(x0=t9s, y0=p.fractions, x1=t9s, y1=vert.obs.99)
     )
max(abs(vert.obs.99-p.fractions))       #approx ks statistic

plot(y=p.fractions-vert.obs.99, x=t9s, type="h", ylim=c(-.15,.15), cex=1)
title(main="deviations: H0 is true", cex=.9)



## one sample when H0 is false

if.R(s=
     ks.gof(t9, distribution="t", df=2)
     ,r=
     ks.test(t9, y="pt", df=2)
     )
t2.quantiles <- qt(p.fractions, 2)

plot(x=range(t9, t2.quantiles), y=c(0,1), type="n", cex=1)
title(main="KS for H0: t9 ~ t with 2 df", cex=.9)
lines(x=t2.quantiles, y=p.fractions)    #reference distribution
points(x=t9s, y=p.fractions, pch=16, cex=.3)

vert.obs.92 <- pt(t9s, 2)               #additional reference distribution
lines(x=t9s, y=vert.obs.92, lty=2)      #to be used for finding ks statistic

if.R(s=
     segments(x1=t9s, y1=p.fractions, x2=t9s, y2=vert.obs.92)
     ,r=
     segments(x0=t9s, y0=p.fractions, x1=t9s, y1=vert.obs.92)
     )
max(abs(vert.obs.92-p.fractions))       #approx ks statistic

plot(y=p.fractions-vert.obs.92, x=t9s, type="h", ylim=c(-.15,.15), cex=1)
title(main="deviations: H0 is false", cex=.9)


## two sample

## t3 <- rt(100,3)
## t5 <- rt(100,5)
if.R(s=
     ks.gof(t3, t5)
     ,r=
     ks.test(t3, t5)
     )

t3s <- sort(t3)
t5s <- sort(t5)

tpg <- trellis.par.get("superpose.line")

matplot(x=cbind(t3s, t5s), y=(1:100)/101,
        type="n",  xaxt="n", yaxt="n", xlab="", ylab="")
matpoints(x=cbind(t3s, t5s), y=(1:100)/101,
        pch=c(16,17), cex=.3, col=tpg$col)
title(main="two-sample KS", cex=.9)
title(xlab="t3 and t5",
      ylab="empirical distribution",
      cex=1)
axis(1, cex=1)
axis(2, cex=1)

matplot(x=cbind(t3s, t5s), y=(1:100)/101,
        type="n",  xaxt="n", yaxt="n", xlab="", ylab="")
matlines(x=cbind(t3s, t5s), y=(1:100)/101,
         lty=c(1,6), col=tpg$col)
title(main="two-sample KS", cex=.9)
title(xlab="t3 and t5",
      ylab="empirical distribution",
      cex=1)
axis(1, cex=1)
axis(2, cex=1)

par(old.par)
par(mfcol=c(1,1))
## export.eps(hh("iinf/figure/iinf.f3.ps"))


#### iinf/code/goffit-t9-t3-t5.s
## used by iinf/code/goffit.s, and not directly by user.
