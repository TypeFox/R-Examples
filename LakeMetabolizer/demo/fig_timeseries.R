#figure generation code for time series of Sparkling Lake (fig. 1)
library(LakeMetabolizer)

data.path = system.file('extdata/', package="LakeMetabolizer")
sp.data = load.all.data('sparkling', data.path)
ts.data = sp.data$data #pull out just the timeseries data


#calculate U10 and add it back onto the original
u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

o2.sat = o2.at.sat(ts.data[,c('datetime','wtr_0')])

ts.data = merge(o2.sat, ts.data)
z.mix = ts.meta.depths(get.vars(ts.data, 'wtr'), seasonal=TRUE)
names(z.mix) = c('datetime','z.mix', 'bottom')

#set z.mix to bottom of lake when undefined
z.mix[z.mix$z.mix <=0 | is.na(z.mix$z.mix), 'z.mix'] = 20 
ts.data = merge(ts.data, z.mix[,c('datetime','z.mix')])

wtr.data<-get.vars(ts.data,var.names = c('datetime','wtr')) # pulling out just water data 
data.path = system.file('extdata', package="rLakeAnalyzer") 
bathy<-load.bathy(fPath = file.path(data.path,'Sparkling.bth')) 

ts.data$schmidt<-ts.schmidt.stability(wtr = wtr.data,bathy = bathy)$schmidt.stability # calculating schmidt stability 
ts.data$meta.bot<-ts.meta.depths(wtr=wtr.data)$bot
ts.data$thermo.depth<-ts.thermo.depth(wtr.data)$thermo.depth



add_axes <- function(xlim, ylim, ylabel = pretty(ylim,10), panel.txt, no.x=TRUE){
  prc_x = 0.14 # position for text relative to axes
  prc_y = 0.07
  tick_len <- 0.15
  ext_x <- c(xlim[1]-86400, pretty(xlim,3), xlim[2]+86400)
  ext_y <- c(ylim[1]-10, pretty(ylim,10), ylim[2]+10)
  ylab <- c("",ylabel,"")
  if (is.na(ylabel[1])) ylab = NA
  #if(no.x) ext_x=NA
  if(no.x){
    axis(side = 1, at = ext_x , labels = FALSE, tcl = tick_len)
  }else{
    axis(side = 1, at = ext_x , labels = strftime(ext_x,'%m-%d'), tcl = tick_len)  
  }
  axis(side = 2, at = ext_y, labels = ylab, tcl = tick_len)
  axis(side = 3, at = ext_x, labels = NA, tcl = tick_len)
  axis(side = 4, at = ext_y, labels = NA, tcl = tick_len)
  x_txt <- (xlim[2] - xlim[1])*prc_x+xlim[1]
  y_txt <- ylim[2]-(ylim[2] - ylim[1])*prc_y
  text(x = x_txt, y_txt,labels = panel.txt)
}


# Figure 1 
cols <- c("#1b9e77", "#d95f02", "black", "#e7298a", "DodgerBlue", "#e6ab02", "grey50")
width = 3.37 # single column width for journal
night_col = 'grey90'
height = 8
l_mar = 0.35
b_mar = 0.1
t_mar = 0.05
r_mar= 0.08
gapper = 0.15 # space between panels

# figure saved at home folder 
png('~/fig_1.png', res=300, width=width, height=height, units = 'in')

#layout(matrix(c(rep(1,10),rep(2,9)),ncol=1)) # 55% on the left panel
par(mai=c(b_mar,l_mar,t_mar,0), omi = c(0.1,0,0,r_mar),xpd=FALSE,
    mgp = c(1.15,.05,0), mfrow=c(5,1))

#Plot the  time series data 
ylim = c(max(ts.data$meta.bot),0)
xlim = as.POSIXct(c('2009-07-01 16:00', '2009-07-11'))
plot(ts.data$thermo.depth~ts.data$datetime, type='l',
     col=cols[3], ylim=ylim,xaxt = 'n', ylab=expression(Mixed~Layer~'&'~Thermocline~(m)),
     xlab='', axes=FALSE)
lines(ts.data$z.mix~ts.data$datetime,col=cols[7],lty=2)
lines(ts.data$meta.bot~ts.data$datetime,col=cols[7],lty=2)
add_axes(xlim, ylim, panel.txt='a)')
add_axes(xlim,rev(ylim),ylabel=NA,panel.txt='')


# schmidt 
ylim = c(min(ts.data$schmidt),max(ts.data$schmidt))
plot(ts.data$schmidt~ts.data$datetime, type='l',
     col=cols[1], ylim=ylim,xaxt = 'n', ylab=expression(Schmidt~Stability~(J~m^-2)),
     xlab='', axes=FALSE)
add_axes(xlim, ylim, panel.txt='b)')

# PAR 
ylim = c(min(ts.data$par),max(ts.data$par))
plot(ts.data$par~ts.data$datetime, type='l',
     col=cols[6], ylim=ylim,xaxt = 'n', ylab=expression(PAR~(paste(mu,mol,sep='')~m^-2~s^-1)),
     xlab='', axes=FALSE)
add_axes(xlim, ylim, panel.txt='c)')

# U10 
ylim = c(min(ts.data$wnd_10),max(ts.data$wnd_10))
plot(ts.data$wnd_10~ts.data$datetime, type='l',
     col=cols[7], ylim=ylim,xaxt = 'n', ylab=expression(U10~(m~s^-1)),
     xlab='', axes=FALSE)
add_axes(xlim, ylim, panel.txt='d)')

# DO deviation from saturation 
ts.data$do.dev<-ts.data$doobs_0.5-ts.data$do.sat # deviation of DO obs from saturation 
ylim = c(min(ts.data$do.dev),max(ts.data$do.dev))
plot(ts.data$do.dev~ts.data$datetime, type='l',
     col=cols[5], ylim=ylim,xaxt = 'n', ylab=expression(DO~-~O[s]~(mg~O[2]~L^-1)),
     xlab='', axes=FALSE)
abline(0,0,col=rgb(0,0,0,0.5),lty=2)
add_axes(xlim, ylim, panel.txt='e)',no.x=F)

dev.off()

