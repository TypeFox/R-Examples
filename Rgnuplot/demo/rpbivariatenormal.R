library(rpanel)
#plots a univariate normal distribution

#Initialize the gnuplot handles
h1 <- Gpinit()
h2 <- Gpinit()

#set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

#change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
#read gnuplot script
gpfile <- system.file('extdata/bivariatenormal.gnu', package='Rgnuplot')
s2 <- Gpfile2string(gpfile)

Gpcmd(h1,'flag=1\nmx=0.0\nsigx=(0.77)\nmy=0.61\nsigy=(1.725)\nsc=0.5\n' %s% s2)
Gpcmd(h2,'flag=0\nmx=0.0\nsigx=(0.77)\nmy=0.61\nsigy=(1.725)\nsc=0.5\n' %s% s2)

mx <- 0.0
sigx <- 0.77
my <- 0.61
sigy <- 1.725
sc <- 0.5

plot.binomial <- function(panel) {
   with(panel, {

s1 <- 'mx = ' %s% mx %s% '
sigx = ' %s% sigx %s% '
my = ' %s% my %s% '
sigy = ' %s% sigy %s% '
sc = ' %s% sc

Gpresetplot(h1)
Gpcmd(h1,'flag = 1\n' %s% s1 %s% '\n' %s% s2)
Gpresetplot(h2)
Gpcmd(h2,'flag = 0\n' %s% s1 %s% '\n' %s% s2)
   
      })
   panel
   }

rp.binomial <- function() {
   #x11(width = 3, height = 3)size = c(400, 400), 
   pname <- rp.control("Bivariate normal distribution", mx = 0.0,
   sigx = 0.77,
my = 0.61,
sigy = 1.725,
sc = 0.5, dev.number = dev.cur())
   rp.slider(pname, mx, -5, 5, initval = 0.0, title = "ux:", 
                   action = plot.binomial)
   rp.slider(pname, sigx, 0, (5), initval = (0.77), title = "var(x):", 
                   action = plot.binomial)               
   rp.slider(pname, my, -5, 5, initval = 0.61, title = "uy:", 
                   action = plot.binomial)
   rp.slider(pname, sigy, 0, (5), initval = (1.725), title = "var(y):", 
                   action = plot.binomial)               
   rp.slider(pname, sc, -1, 1, initval = .5, title = "p:", 
                   action = plot.binomial)    
rp.button(pname, { h1<-Gpclose(h1);h2<-Gpclose(h2);break } , title = 'Close', quitbutton=TRUE)
   }

rp.binomial()

#close gnuplot handles
#h1<-Gpclose(h1);h2<-Gpclose(h2)   
