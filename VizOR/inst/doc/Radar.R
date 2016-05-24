### R code from vignette source 'Radar.Rnw'

###################################################
### code chunk number 1: Init
###################################################
library(rms)
library(VizOR)


###################################################
### code chunk number 2: PrepareData
###################################################
df <- upData(mtcars,
             cyl=factor(cyl,levels=2*(2:4),labels=paste(2*(2:4),"cyl", sep="-")),
             am=factor(am,levels=0:1,labels=c("automatic","manual")),
             gear=factor(gear,levels=3:5,labels=paste(3:5,"speed", sep="-")),
             ## TODO: Add a region factor?
             labels=c(
               mpg="Miles per gallon"
               ,cyl="Number of cylinders"
               ,disp="Displacement"
               ,hp="Gross horsepower"
               ,drat="Rear axle ratio"
               ,wt="Weight"
               ,qsec="1/4 mile time"
               ,am="Transmission type"
               ,gear="Number of forward gears"
               ,carb="Number of carburetors"
               ),
             units=c(
               wt="lb/1000"
               ,disp="in^3"
               ,qsec="sec"
               ),
             drop='vs' # I have no idea what this poorly documented variable means!
             )


###################################################
### code chunk number 3: Describe
###################################################
latex(describe(df, descript="Built-in dataset `mtcars'"), file="")


###################################################
### code chunk number 4: ColorScheme
###################################################
## There may be an excellent opportunity here to develop and demonstrate
## an approach to colorizing factors.  Each factor employed in the 3-way
## experimental cross can be given its own color scheme, to be applied
## automatically when plotting.


###################################################
### code chunk number 5: Table1
###################################################
## Here, demonstrate a basic "Table 1" which subsequently will appear
## in the individual panels of the trellised radar plot.
s <- summary(cyl ~ mpg + disp + hp + drat + wt, method='reverse', data=df)
options(digits=3)
latex(s, npct='both', nptc.size='normalsize', file="", label="tbl:Table-1")


###################################################
### code chunk number 6: CrossPlot
###################################################
## Let's keep the overall=TRUE option (which is the default for method='cross'),
## so that a single summary can be generated and used for many plots.  The marginal
## factor level should probably be ignored automatically for the overlaid polygons.
s <- summary(cbind(mpg, disp, hp, drat, wt) ~ cyl + gear + am,
             method='cross', overall=TRUE, data=df)
dd <- datadist(df)
print(radarplot(S ~ cyl | gear*am, data=s, datadist=dd, rescale="range"))


