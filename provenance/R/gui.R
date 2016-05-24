#' Menu-based interface for \code{provenance}
#'
#' For those less familiar with the syntax of the \code{R}
#' programming language, the \code{provenance()} function provides
#' a user-friendly way to access the most important functionality
#' in the form of a menu-based query interface. Further details
#' and examples are provided on \url{http://provenance.london-geochron.com}
#' @author Pieter Vermeesch
#' @references Vermeesch, P., Resentini, A. and Garzanti, E., an R
#' package for statistical provenance analysis, Sedimentary Geology,
#' doi:10.1016/j.sedgeo.2016.01.009.
#' @seealso \url{http://provenance.london-geochron.com}
#' @export
provenance <- function(){
    version <- as.character(utils::packageVersion('provenance'))
    message("This is provenance version ", version)
    while (TRUE){
        message("Pick an option:\n",
                "1 - sample size calculation\n",
                "2 - plot a single dataset\n",
                "3 - plot multiple datasets\n",
                "4 - Minsorting\n",
                "5 - MDS/PCA\n",
                "6 - Procrustes analysis\n",
                "7 - 3-way MDS\n",
                "8 - save plots (.pdf)\n",
                "9 - help\n",
                "q - quit")
        response <- readline()
        if (response == '1'){
            gui.samplesize()
        } else if (response == '2'){
            gui.plot.single()
        } else if (response == '3'){
            gui.plot.multiple()
        } else if (response == '4'){
            gui.minsorting()
        } else if (response == '5'){
            gui.mds()
        } else if (response == '6'){
            gui.3way(TRUE)
        } else if (response == '7'){
            gui.3way(FALSE)
        } else if (response == '8'){
            gui.save.plots()
        } else if (response == '9'){
            gui.help()
        } else if (response == 'q'){
            break
        }
    }
}

# sigdig = number of significant digits
gui.samplesize <- function(sigdig=3){
    message("Sample size calculations:\n",
            "1 - Probability of missing a fraction of given size\n",
            "2 - Largest fraction not to be missed with a given confidence\n",
            "3 - Number of grains needed to sample all fractions of a given size")
    response <- readline()
    if (response == '1'){
        n <- readline("enter the number of grains [e.g., 117]: ")
        f <- readline(paste0("enter a number between 0 and 100% indicating\n",
                             "the size of the smallest fraction of interest [e.g., 5]: "))
        p <- get.p(as.numeric(n),as.numeric(f)/100)
        message("\nResult: the likelihood that all fractions greater than ", f , "%\n",
                "of the population are missed in a sample with ", n ," grains is ",
                "p = ",signif(100*p,sigdig) , "%\n")
    } else if (response == '2'){
        n <- readline("enter the number of grains [e.g., 117]: ")
        p <- readline(paste0("enter a number between 0 and 100% indicating\n",
                             "the desired level of confidence [e.g., 95]: "))
        f <- get.f(as.numeric(n),as.numeric(p)/100)
        message("\nResult: the largest fraction of which a ", n , "-grain sample\n",
                "has not missed with ", p , "% confidence is f = ", signif(100*f,3),"%\n")
    } else if (response == '3'){
        f <- readline(paste0("enter a number between 0 and 100% indicating\n",
                             "the size of the smallest fraction of interest [e.g., 5]: "))
        p <- readline(paste0("enter a number between 0 and 100% indicating\n",
                             "the desired level of confidence [e.g., 95]: "))
        n <- get.n((100-as.numeric(p))/100,as.numeric(f)/100)
        message("\nResult: the minimum sample size which guarantees with ",p,"% certainty not\n",
                "to have missed any fraction greater than ",f,"% of the population is n = ",n,"\n")
    } else {
        message('Incorrect input.')
    }
}

gui.plot.single <- function(){
    grDevices::graphics.off()
    message("Plot a single dataset:\n",
            "1 - Ternary diagram\n",
            "2 - Pie charts\n",
            "3 - Cumulative Age Distributions\n",
            "4 - Kernel Density Estimates")
    response1 <- readline()
    if (response1 == '1'){
        showpath = FALSE
        dat <- gui.open.compositional()
        message("Plot background lines?\n",
                "1 - descriptive QFL diagram\n",
                "2 - Folk's classification\n",
                "3 - Dickinson's QFL diagram\n",
                "4 - no lines")
        response2 <- readline()
        if (response2 == '1'){
            type <- "QFL.descriptive"
        } else if (response2 == '2'){
            type <- "QFL.folk"
        } else if (response2 == '3'){
            type <- "QFL.dickinson"
        } else {
            type <- "empty"
        }
        if (methods::is(dat,"SRDcorrected")){
            response3 <- readline("Show SRD correction [Y or n]? ")
            if (response3 %in% c('n','N')) showpath <- FALSE
            else showpath <- TRUE
        }
        graphics::plot(ternary(dat),type=type,showpath=showpath)
    } else if (response1 == '2'){
        dat <- gui.open.compositional()
        numcol <- as.numeric(readline("Number of columns in the plot? "))
        if (is.na(numcol)) numcol <- 1
        summaryplot(dat,ncol=numcol)
    } else if (response1 == '3'){
        dat <- gui.open.distributional()
        graphics::plot(dat,CAD=TRUE)
    } else if (response1 == '4'){
        dat <- gui.open.distributional()
        kdes <- gui.get.kdes(dat)
        if (methods::is(kdes,"KDE")){
            graphics::plot(kdes)
        } else { # is(kdes,'KDEs')
            numcol <- as.numeric(readline("Number of columns in the plot? "))
            if (is.na(numcol)) numcol <- 1
            summaryplot(kdes,ncol=numcol)
        }
    } else {
        message('Incorrect input.')
    }
}

gui.plot.multiple <- function(){
    command <- "summaryplot("
    datasets <- gui.get.datasets(kdes=TRUE)
    dnames <- names(datasets)
    for (dname in dnames){
        command <- paste0(command,"datasets[[\'",dname,"\']],")
    }
    numcol <- as.numeric(readline("Number of columns? "))
    if (is.na(numcol)) numcol <- 1
    command <- paste0(command,"ncol=",numcol,")")
    eval(parse(text=command))
}

gui.minsorting <- function(){
    phi <- 2
    sigmaphi <- 1
    medium <- "seawater"
    fromphi <- -2.25
    tophi <- 5.5
    byphi <- 0.05
    endcomp <- provenance::endmembers
    mydens <- provenance::densities
    while (TRUE){
        message("Change default:\n",
                "1 - bulk composition [default=tectonic endmembers]\n",
                "2 - average grain size in phi units [default=",phi,"]\n",
                "3 - grain size standard deviation [defaults=",sigmaphi,"]\n",
                "4 - transport medium [default=",medium,"]\n",
                "5 - plot resolution [default: from ",fromphi," to ", tophi," by ",byphi,"]\n",
                "6 - selected minerals [default: all]\n",
                "c - continue")
        response1 <- readline()
        if (response1 == "1") {
            message('Open a compositional dataset:')
            endcomp <- read.compositional(file.choose())
        } else if (response1 == "2"){
            mydens <- gui.load.densities()
        } else if (response1 == "3"){
            phi <- as.numeric(readline("Enter the average grain size (in phi units): "))
        } else if (response1 == "3"){
            message("Enter the standard deviation of the grain size (in phi units): ")
            sigmaphi <- as.numeric(readline())
        } else if (response1 == "4"){
            message("Choose one of the following options:\n",
                    "1 - seawater\n",
                    "2 - freshwater\n",
                    "3 - air")
            response2 <- readline()
            if (response2 == "2"){
                medium <- "freshwater"
            } else if (response2 == "3"){
                medium <- "air"
            } else {
                medium <- "seawater"
            }
        } else if (response1 == "5"){
            message("Change:\n",
                    "1 - minimum [default = ",minphi,"]\n",
                    "2 - maximum [default = ",maxphi,"]\n",
                    "3 - step size [default = ",byphi,"]")
            response2 <- readline()
            if (response2 == "1"){
                minphi <- as.numeric(response2)
            } else if (response2 == "2"){
                maxphi <- as.numeric(response2)
            } else if (response2 == "3"){
                byphi <- as.numeric(response2)
            } else {
                message("Invalid input")
            }
        }  else if (response1 == "6"){
            endcomp <- gui.subset.components(endcomp)
        } else {
            break
        }
    }
    while (TRUE){
        sname <- NULL
        if (length(names(endcomp)) > 1){
            samplist <- paste(names(endcomp),collapse=',')
            message("Select one of the following samples to plot:\n",
                    samplist,
                    "\nor press [Return] to exit")
            smp <- readline()
            if (smp != "") sname <- smp
        }
        grDevices::graphics.off()
        graphics::plot(minsorting(endcomp,mydens,sname=sname,
                                  phi,sigmaphi,medium,fromphi,tophi,byphi))
        if (length(names(endcomp)) == 1 | smp == "") break
    }
}

gui.mds <- function(){
    grDevices::graphics.off()
    dat <- gui.get.dataset()
    redo <- FALSE
    classical <- FALSE
    method <- dat$method
    if (method=="aitchison"){
        message("1 - Use Aitchison distance (PCA, default)")
        message("2 - Use Aitchison distance (MDS)")
        message("3 - Use Bray-Curtis dissimilarity (MDS)")
        response <- readline()
        if (response == "2"){
            # do nothing
        } else if (response == "3"){
            method <- "bray"
        } else {
            graphics::plot(PCA(dat))
            return()
        }
    }
    if (methods::is(dat,"distributional") & length(dat$x$err)>0){
        message("Choose:\n",
                "1 - Kolmogorov-Smirnov distance\n",
                "2 - Sircombe-Hazelton distance")
        response <- readline()
        if (response == "1") method <- "KS"
        else method <- "SH"
        if (response %in% c("y","Y")) method <- "SH"
    }
    mymds <- MDS(dat,classical,method=method)
    if (mymds$stress < 0.05){
        message("Warning: extremely low stress (=",mymds$stress,").\n",
                "Looks like you\'re overfitting your data!")
    }
    if (!dat$method %in% c("bray","SH")){
        if (mymds$stress < 0.05){
            message("Do you want to switch to classical scaling?")
            response <- readline("[Y or n]: ")
            if (response %in% c("n","N")) classical <- FALSE
            else classical <- TRUE
        }
        if (classical){
            mymds <- MDS(dat,classical,method=method)
        }
    }
    thennlines=FALSE
    thepch=NA
    thecex=1
    thepos=NULL
    thexlab=""
    theylab=""
    thexaxt='n'
    theyaxt='n'
    while (TRUE){
        message("Options:\n",
                "1 - Add nearest neighbour lines\n",
                "2 - Change plot character\n",
                "3 - Change size of plot character\n",
                "4 - Change position of text label relative to plot character\n",
                "5 - Add X and Y axis ticks\n",
                "c - Continue")
        response1 <- readline()
        if (response1 == "1"){
            thennlines <- TRUE
        } else if (response1 == "2"){
            thepch <- readline(
                "Specify a plot character [e.g. *, ., o, O, +, or 1-25]: ")
            if (grepl("[[:digit:]]",thepch)) thepch <- as.numeric(thepch)
        } else if (response1 == "3"){
            thecex <- as.numeric(readline(
                "Magnification of the default plot character [1 = normal]: "))
        } else if (response1 == "4"){
            thepos <- as.numeric(readline(
                "Position of the text label [1 = below, 2 = left, 3 = above, 4 = right]"))
        } else if (response1 == "5"){
            thexlab <- "X"
            theylab <- "Y"
            thexaxt <- 's'
            theyaxt <- 's'
        } else {
            break
        }
    }
    graphics::plot(mymds,nnlines=thennlines,pch=thepch,pos=thepos,cex=thecex,
                   xlab=thexlab,ylab=theylab,xaxt=thexaxt,yaxt=theyaxt)
}

gui.3way <- function(doprocrustes=FALSE){
    grDevices::graphics.off()
    if (doprocrustes) command <- "procrustes("
    else command <- "indscal("
    datasets <- gui.get.datasets()
    dnames <- names(datasets)
    for (dname in dnames){
        command <- paste0(command,dname,"=","datasets[[\'",dname,"\']],")
    }
    substr(command,start=nchar(command),stop=nchar(command)) <- ")"
    proc <- eval(parse(text=command))
    graphics::plot(proc)    
}

gui.save.plots <- function(){
    devices <- grDevices::dev.list()
    for (device in devices){
        fname <- readline(paste0("Enter a name for plot ",device, ": "))
        grDevices::dev.set(device)
        grDevices::dev.copy2pdf(file=paste0(fname,".pdf"))
    }
}

gui.help <- function(){
    while (TRUE){
    message("Choose a help topic:\n",
            "1 - compositional data\n",
            "2 - distributional data\n",
            "3 - mineral densities\n",
            "4 - sample size calculation\n",
            "5 - cumulative age distributions\n",
            "6 - kernel density estimation\n",
            "7 - Minsorting and SRD correction\n",
            "8 - multidimensional scaling\n",
            "9 - Procrustes analysis and 3-way MDS\n",
            "c - continue")
        response <- readline()
        if (response == "1"){
            message("Formatting requirements of compositional data files:\n",
                    "---------------------------------------------------\n",
                    "A compositional data file consists of a comma separated table\n",
                    "in which the rows correspond to samples and the columns to categories.\n",
                    "For example:\n",
                    "sample,zr,tm,rt,TiOx,sph,ap,ep,gt,st,and,ky,sil,amp,cpx,opx\n",
                    "N1,1,1,0,0,0,1,10,24,2,0,0,0,1,163,1\n",
                    "N2,0,0,1,0,3,1,11,18,1,0,0,1,3,170,1\n",
                    "N3,1,0,0,1,1,0,9,21,1,0,0,0,10,157,1\n",
                    "N4,0,2,4,1,4,0,59,54,0,1,0,0,19,54,2\n",
                    "N5,1,2,2,5,12,2,91,43,5,2,0,0,32,10,1\n",
                    "N6,3,1,2,1,2,0,11,12,1,0,0,0,31,139,3\n",
                    "N7,4,0,1,0,0,0,12,19,0,1,0,0,10,155,1\n",
                    "N8,5,0,0,0,2,0,23,46,1,0,0,1,11,113,0\n",
                    "N9,10,2,0,0,9,0,26,80,2,0,0,0,18,56,0\n",
                    "N10,2,0,0,0,4,0,16,28,2,0,1,0,34,126,0\n",
                    "If you want to use the Minsorting function and SRD correction,\n",
                    "it is important to make sure that the component labels are consistent with\n",
                    "the densities table.\n")
        } else if (response == "2"){
            message("Formatting requirements for distributional data files:\n",
                    "-----------------------------------------------------\n",
                    "A distributional data file consists of a comma separated table\n",
                    "in which the columns correspond to different samples, each of which\n",
                    "may contain a different number of single grain analyses.\n",
                    "For example:\n",
                    "N1,N2,N3,N4,N5,N6,N7,N8,N9,N10\n",
                    "645.4,2764.6,418.6,1228.4,1227,605.1,931.2,2440.3,2117.6,2727.2\n",
                    "496.9,1998.6,1036.1,1445.7,1269.7,1462.6,980,1246.4,1919.1,768\n",
                    "1000.3,997.7,1118.4,1400.5,2127.1,1240.1,1383.9,1131.9,640.9,1294.2\n",
                    "1168.5,1105.5,1221.4,1053.9,1185.1,316.5,2017.4,1116.2,679.7,1462.1\n",
                    "2263.5,620.6,944.7,2073.4,362.6,1524.5,1353.8,1985.7,732.3,1924\n",
                    "1878.6,1133.7,937.8,1107.1,1114.6,1141.3,948.2,1055.6,1190.3,1154.3\n",
                    "769.8,583.6,982,1386.1,701.2,1521.5,618.5,296.2,2005.4,1350.8\n",
                    "1161.7,1216.6,752.6,1221.8,1064.5,1166.7,2104.7,556.1,529.7,538.5\n",
                    "519.4,1181.8,1240.3,871,1133.7,660.1,307.2,815.9,2088.8,1412.2\n",
                    "1213.3,1203.7,1278.2,1418,938.7,675.5,2052.9,2266.4,1584.3,1577.6\n",
                    "271.3,1251.2,1088.2,1853.9,991.8,1211.4,1977.3,1094.6,995.4,2809.5\n",
                    "1065.4,1129.1,2041,2905.9,2039.6,1391.2,1064.2,634.9,1069,659.3\n",
                    "1114.6,,577.8,1928.9,1064.5,2076.7,703,981.3,1858.4,1302.3\n",
                    "998.7,,976.9,728.5,698.1,,1137.6,1088,1157.4,1135.8\n",
                    "603.3,,1411.5,1371.7,1507.9,,1325,670.6,215.9,973.3\n",
                    "465.9,,547.5,814.7,668,,583.8,663.5,,1194.1\n",
                    "489.4,,791.5,,1279.3,,1120.2,1767.5,,1157.1\n",
                    "744.5,,2160.2,,1386.9,,1989.2,,,811.5\n",
                    "1287.4,,1043.6,,1210.8,,1163.4,,,1554.9\n",
                    "It is also possible (though optional) to specify the analytical\n",
                    "uncertainties of these measurements. They are stored in a separate\n",
                    ".csv file with exactly the same format as the measurement as shown above.\n")
        } else if (response == "3"){
            message("Formatting requirements for mineral density files:\n",
                    "-------------------------------------------------\n",
                    "The Minsorting function and SRD correction require a table of\n",
                    "mineral densities. The provenance package comes with a default table\n",
                    "but it is also possible to specify a different set of densities. It\n",
                    "is important that the labels used for the different components in this\n",
                    "table are consistent with the compositional data files used elsewhere.\n ",
                    "For example:\n",
                    "Q,F,Lv,Ls,Lm,mica,opaques,zr,tm,rt,sph,ap,\n",
                    "mon,oth,ep,gt,ctd,st,and,ky,sil,amp,px,ol\n",
                    "2.65,2.61,2.6,2.65,2.75,2.4,5,4.65,3.15,4.25,3.5,3.2,\n",
                    "5.15,3.5,3.45,4,3.6,3.75,3.15,3.6,3,3.2,3.3,3.35\n")
        } else if (response == "4"){
            message("Sample size calculations:\n",
                    "------------------------\n",
                    "On the most basic level, provenance analysis requires the geologist to\n",
                    "identify certain properties in a representative number of grains from\n",
                    "each sample. The question then arises how many grains constitute a\n",
                    "'representative' number of grains. The answer to this question depends\n",
                    "on the geological problem of interest. The sample size calculators implemented\n",
                    "in provenance() assume that we are interested in collecting at least one grain\n",
                    "from every fraction greater than f=1/N, where N is an integer denoting the\n",
                    "number of fractions. Using simple combinatorics, it is possible to calculate:\n",
                    "(a) the number of grains n required to ensure that all N fractions were sampled\n",
                    "at least once with p% certainty.\n",
                    "(b) the largest fraction f of which it is p% certain that all were sampled\n",
                    "at least one out of n grains.\n",
                    "(c) the likelihood that a sample of size n has missed at least one\n",
                    "fraction of size f.\n",
                    "\nCiteable references:\n",
                    "Vermeesch, P., 2004. How many grains are needed for a provenance study?\n",
                    "Earth and Planetary Science Letters 224, 441-451.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n")
        } else if (response == "5"){
            message("Cumulative Age Distributions (CADs):\n",
                    "-----------------------------------\n",
                    "The probability distribution of detrital zircon U-Pb ages or any other\n",
                    "distributional dataset may be estimated using histograms or kernel density\n",
                    "estimates (KDEs). For a sample of limited size, these estimates never exactly\n",
                    "agree with the true age distribution, but are smooth approximation thereof.\n",
                    "In contrast, the Cumulative Age Distribution (also known as the 'Empirical\n",
                    "Cumulative Distribution Function) is a method to visualise distributional\n",
                    "datasets without the need for any smoothing.\n",
                    "\nCiteable references:\n",
                    "Vermeesch, P., 2007. Quantitative geomorphology of the White Mountains\n",
                    "(California) using detrital apatite fission track thermochronology.\n",
                    "Journal of Geophysical Research (Earth Surface) 112, 3004. doi:10.1029/2006JF000671.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n")
        } else if (response == "6"){
            message("Kernel Density Estimation (KDE):\n",
                    "-------------------------------\n",
                    "A large body of statistical literature has been written on the subject\n",
                    "of probability density estimation. The most widely used tools are the\n",
                    "traditional histogram, which is discrete and discontinuous, and the Kernel\n",
                    "Density Estimator (KDE), which is a smooth and continuous alternative.\n",
                    "A KDE is produced by arranging the measurements along a line and stacking\n",
                    "a so-called 'kernel' (typically a Gaussian bell curve) of a certain\n",
                    "width (the 'bandwidth') on top of each of them. Multiple samples can be\n",
                    "plotted together, either using automatically selected bandwidths using the\n",
                    "algorithm of Botev et al. (2010), or any other value\n",
                    "\nCiteable references:\n",
                    "Botev, Z. I., Grotowski, J. F., Kroese, D. P., 2010. Kernel density\n",
                    "estimation via diffusion. Annals of Statistics 38, 2916-2957.\n",
                    "Vermeesch, P., 2012. On the visualisation of detrital age distributions.\n",
                    "Chemical Geology 312-313, 190-194. doi:10.1016/j.chemgeo.2012.04.021.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n")
        } else if (response == "7"){
            message("Minsorting and SRD correction:\n",
                    "-----------------------------\n",
                    "To facilitate the comparison of detrital modes for provenance analysis\n",
                    "or stratigraphic correlation, we need to first remove the often significant\n",
                    "compositional differences among sediment samples that are caused by\n",
                    "hydrodynamic processes in the depositional environment. Intersample modal\n",
                    "variability can be corrected for by a simple principle. In the absence of\n",
                    "provenance changes and environmental bias, the weighted average\n",
                    "Source Rock Density (SRD) of terrigenous grains should be equal,\n",
                    "for each sample and each grain-size class of each sample, to the\n",
                    "weighted average density of source rocks. By correcting relative\n",
                    "abundances of detrital minerals in proportion to their densities,\n",
                    "we can restore the appropriate SRD index for any provenance and\n",
                    "subprovenance type in each sample or grain-size class (Garzanti et al., 2009).\n",
                    "Within a sample of a certain average grain size, the denser minerals will\n",
                    "tend to be enriched in the finer fraction, and the lighter minerals in the\n",
                    "coarser fraction. This so-called principle of 'settling equivalence'\n",
                    "(Garzanti et al., 2008), is implemented in the Minsorting function, which is\n",
                    "based on the Excel spreadsheet of Resentini et al. (2013).\n",
                    "\nCiteable references:\n",
                    "Garzanti, E., Ando, S., Vezzoli, G., 2009. Grain-size dependence of\n",
                    "sediment composition andenvironmental bias in provenance studies.\n",
                    "Earth and Planetary Science Letters 277, 422-432.\n",
                    "Garzanti, E, Ando, S and Vezzoli, G. Settling equivalence of detrital\n",
                    "minerals and grain-size dependence of sediment composition.\n",
                    "Earth and Planetary Science Letters 273.1 (2008): 138-151.\n",
                    "Resentini, A., Malusa, M.G., Garzanti, E., 2013. MinSORTING: An Excel\n",
                    "worksheet for modelling mineral grain-size distribution in sediments,\n",
                    "with application to detrital geochronology and provenance studies.\n",
                    "Computers & Geosciences 59, 90-97.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n")
        } else if (response == "8"){
            message("Multidimensional Scaling:\n",
                    "------------------------\n",
                    "An increasing number of provenance studies are based on not just a few\n",
                    "but many samples. The large datasets resulting from such studies call\n",
                    "for a dimension-reducing technique such as Multi-Dimensional Scaling (MDS).\n",
                    "Given a dissimilarity matrix (i.e., a table of pairwise distances),\n",
                    "MDS constructs a 'map' on which 'similar' samples cluster closely\n",
                    "together and 'dissimilar' samples plot far apart. provenance() implements\n",
                    "the following dissimilarity measures:\n",
                    "Compositional data:\n",
                    "(a) the Aitchison distance (all components are strictly positive)\n",
                    "(b) the Bray-Curtis distance (some zero components)\n",
                    "Distributional data:\n",
                    "(c) the Kolmogorov-Smirnov distance\n",
                    "(d) the Sircombe-Hazelton distance\n",
                    "All these dissimilarities can be fitted using a 'non-metric' MDS\n",
                    "algorithm, which aims to preserve the rank order of the relative distances.\n",
                    "For dissimilarities (a) and (c), it is also possible to perform\n",
                    "'classical MDS', which aims to represent the actual values of the dissimilarities\n",
                    "Finally, for the Aitchison distance (a) which turns out to be a Euclidean\n",
                    "distance, it can be shown that MDS is equivalent to Principal Component Analysis\n",
                    "(PCA), in which the configuration can be plotted as a 'compositional biplot',\n",
                    "providing insights into the genetic origin of the observed intersample variability\n",
                    "\nCiteable references:\n",
                    "Vermeesch, P., 2013. Multi-sample comparison of detrital age distributions.\n",
                    "Chemical Geology 341, 140-146.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n")
        } else if (response == "9"){
            message("Procrustes analysis and 3-way MDS:\n",
                    "---------------------------------\n",
                    "Procrustes analysis and 3-way MDS are simple yet powerful tools\n",
                    "to extract geological insights from large, multi-method provenance datasets.\n",
                    "Procrustes analysis is the process by which a combination of\n",
                    "shape-preserving transformations is used to match the shape of\n",
                    "one object with that of another. Generalised Procrustes Analysis\n",
                    "(GPA) is a generalisation of this procedure to multiple objects.\n",
                    "In a provenance context, GPA extracts a single 'consensus' view from\n",
                    "a collection of MDS configurations, by rotating, reflecting and scaling\n",
                    "them to minimise a least squares criterion.\n",
                    "As the name suggests, 3-way MDS is a generalisation of the MDS method from\n",
                    "two-to three-dimensional dissimilarity matrices. provenance() implements\n",
                    "the INdividual Differences SCALing (INDCAL) algorithm, which is the oldest\n",
                    "and still most widely used 3-way MDS algorithm (de Leeuw and Mair, 2009).\n",
                    "In contrast with 2-way MDS and GPA, INDSCAL produces not one but two\n",
                    "pieces of graphical output: the 'group configuration' and the 'source weights'.\n",
                    "\nCiteable references:\n",
                    "de Leeuw, J., Mair, P., 2009. Multidimensional scaling using majorization:\n",
                    "The R package smacof. Journal of Statistical Software 31, 1-30.\n",
                    "URL: http://www.jstatsoft.org/v31/i03\n",
                    "Vermeesch, P., Garzanti, E., 2015. Making geological sense of 'Big Data'\n",
                    "in sedimentary provenance analysis. Chemical Geology 409, 20-27.\n",
                    "Vermeesch, P., Resentini, A. and Garzanti, E., 2016. An R package for\n",
                    "statistical provenance analysis. Sedimentary Geology,\n",
                    "doi:10.1016/j.sedgeo.2016.01.009.\n"
                    )
        } else {
            break
        }
    }
}

gui.open.compositional <- function(){
    message('Open a compositional dataset:')
    fname <- file.choose()
    dat <- read.compositional(fname)
    while (TRUE){
        message("1 - Apply SRD correction\n",
                "2 - Subset components\n",
                "3 - Amalgamate components\n",
                "4 - Subset samples\n",
                "c - Continue")
        response <- readline()
        if (response == "1"){
            dat <- gui.SRD(dat)
        } else if (response == "2"){
            dat <- gui.subset.components(dat)
        } else if (response == "3"){
            dat <- gui.amalgamate.components(dat)
        } else if (response == "4"){
            dat <- gui.subset.samples(dat)
        } else {
            break
        }
    }
    dat
}

gui.open.distributional <- function(){
    message('Open a distributional dataset:')
    fname <- file.choose()
    errorfile <- NA
    dat <- read.distributional(fname,errorfile)
    while (TRUE){
        message("Options:\n",
                "1 - Subset samples\n",
                "2 - Combine samples\n",
                "3 - Load analytical uncertainties\n",
                "c - Continue")
        response <- readline()
        if (response == "1"){
            dat <- gui.subset.samples(dat,TRUE)
        } else if (response == "2"){
            dat <- gui.combine.samples(dat)
        } else if (response == "3"){
            errorfile <- file.choose()
            dat <- read.distributional(fname,errorfile)
        } else {
            return(dat)
        }
    }
}

gui.SRD <- function(dat){
    mydens <- gui.load.densities()
    response <- readline("Enter target density in g/cm3: ")
    restore(dat,mydens,target=as.numeric(response))
}

gui.subset.components <- function(dat){
    complist <- paste(colnames(dat$x),collapse=',')
    message("Select a subset group of components from ",
            "the following list:\n", complist ,"\n",
            "Enter as a comma separated list of labels:")
    response <- readline()
    subcomp <- gsub(" ","",response) # get rid of spaces
    subcomp <- paste0("c('",gsub(",","','",subcomp),"')")  # add apostrophes
    eval(parse(text=paste0("subset(dat,components=",subcomp,")")))
}

gui.combine.samples <- function(dat){
    samples <- names(dat$x)
    combinations <- ""
    while (TRUE){
        samplist <- paste(samples,collapse=',')
        message("Select a group of samples from ",
                "the following list:\n", samplist ,"\n",
                "Enter as a comma separated list of labels ",
                "or click [Return] to exit:")
        response <- readline()
        if (response == ""){
            break
        } else {
            subsamples <- strsplit(response,',')[[1]]
            subsamplist <- gsub(" ","",response) # get rid of spaces
            subsamplist <- paste0("c('",gsub(",","','",subsamplist),"')")  # add apostrophes
            sampname <- readline(paste0("Name of the combined samples? "))
            combinations <- paste0(combinations,sampname,"=",subsamplist,",")
            samples <- c(sampname,samples[which(!samples %in% subsamples)])
        }
    }
    leftovers <- samples[which(samples %in% names(dat$x))]
    for (leftover in leftovers){
        combinations <- paste0(combinations,leftover,"=c('",leftover,"'),")
    }
    # replace last comma with bracket
    substr(combinations,start=nchar(combinations),stop=nchar(combinations)) <- ")"
    eval(parse(text=paste0("combine(dat,",combinations)))
}

gui.amalgamate.components <- function(dat){
    components <- colnames(dat$x)
    command <- "amalgamate(dat,"
    while (TRUE){
        complist <- paste(components,collapse=',')
        message("Select a group of components from ",
                "the following list:\n", complist ,"\n",
                "Enter as a comma separated list of labels ",
                "or click [Return] to exit:")
        response <- readline()
        if (response == ""){
            break
        } else {
            subcomponents <- strsplit(response,',')[[1]]
            subcomplist <- gsub(" ","",response) # get rid of spaces
            subcomplist <- paste0("c('",gsub(",","','",subcomplist),"')")  # add apostrophes
            compname <- readline(paste0("Name of the amalgamated component? "))
            command <- paste0(command,compname,"=",subcomplist,",")
            components <- c(compname,components[which(!components %in% subcomponents)])
        }
    }
    leftovers <- components[which(components %in% colnames(dat$x))]
    for (leftover in leftovers){
        command <- paste0(command,leftover,"=c('",leftover,"'),")
    }
    substr(command,start=nchar(command),stop=nchar(command)) <- ")" # replace last comma with bracket
    eval(parse(text=command))
}

gui.subset.samples <- function(dat,returnsubsamp=FALSE){
    samplist <- paste(names(dat),collapse=',')
    message("Select a subset group of samples from ",
            "the following list:\n", samplist ,"\n",
            "Enter as a comma separated list of labels:")
    response <- readline()
    subsamp <- gsub(" ","",response) # get rid of spaces
    subsamp <- paste0("c('",gsub(",","','",subsamp),"')")  # add apostrophes
    if (returnsubsamp) return(subsamp)
    eval(parse(text=paste0("subset(dat,select=",subsamp,")")))
}

gui.get.kdes <- function(dat){
    from <- NA
    to <- NA
    bw <- NA
    adaptive <- TRUE
    logscale <- FALSE
    normalise <- FALSE
    samebandwidth <- FALSE
    message1 <- paste0("Options:\n",
                       "1 - Set minimum age\n",
                       "2 - Set maximum age\n",
                       "3 - Turn off adaptive density estimation\n",
                       "4 - Plot on a log scale\n",
                       "5 - Set bandwidth\n")
    message2 <- paste0("6 - Use the same bandwidth for all samples\n",
                       "7 - Normalise area under the KDEs\n")
    message3 <- "c - Continue"
    while (TRUE){
        if (length(names(dat))>1){
            message(message1,message2,message3)
        } else {
            message(message1,message3)
        }
        response <- readline()
        if (response == "1"){
            from <- as.numeric(readline("Enter a minimum age limit for the time axis: "))
        } else if (response == "2"){
            to <- as.numeric(readline("Enter a maximum age limit for the time axis: "))
        } else if (response == "3"){
            adaptive <- FALSE
        } else if (response == "4"){
            logscale <- TRUE
        } else if (response == "5"){
            bw <- as.numeric(readline("Enter a custom value for the bandwidth: "))
        } else if (response == "6"){
            samebandwidth <- TRUE
        } else if (response == "7"){
            normalise <- TRUE
        } else {
            break
        }
    }
    if (length(names(dat))>1){ # more than one sample?
        out <- KDEs(dat,from=from,to=to,bw=bw,adaptive=adaptive,
                    normalise=normalise,log=logscale,samebandwidth=samebandwidth)
    } else {
        out <- KDE(dat$x[[1]],from=from,to=to,bw=bw,adaptive=adaptive,log=logscale)
    }
    out
}

gui.load.densities <- function(){
    response <- readline("Open a density file [y] or use default values [N]? ")
    if (response %in% c('Y','y')){
        fname <- file.choose()
        out <- read.densities(fname)
    } else {
        out <- provenance::densities
    }
    out
}

gui.get.datasets <- function(multiple=TRUE,kdes=FALSE){
    datasets <- list()
    while (TRUE){
        if (multiple){
            message("1 - Add a compositional dataset\n",
                    "2 - Add a distributional dataset\n",
                    "c - Continue")
        } else {
            message("1 - Load a compositional dataset\n",
                    "2 - Load a distributional dataset")
        }
        response <- readline()
        if (response == '1'){
            dat <- gui.open.compositional()
        } else if (response == '2'){
            dat <- gui.open.distributional()
            if (kdes) dat <- gui.get.kdes(dat)
        } else {
            break
        }
        if (multiple) datasets[[dat$name]] <- dat
        else return(dat)
    }
    datasets
}

gui.get.dataset <- function(){
    gui.get.datasets(multiple=FALSE)
}
