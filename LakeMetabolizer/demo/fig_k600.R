
# Generates a figure showing the results of the different k600 models
library(LakeMetabolizer)

data.path = system.file('extdata', package="LakeMetabolizer")

sp.data = load.all.data('sparkling', data.path)

ts.data = sp.data$data #pull out just the timeseries data

#calculate U10 and add it back onto the original

u10 = wind.scale(ts.data)
ts.data = rmv.vars(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset

#Calculate k600 using the k.cole wind-based model
k600_cole = k.cole(ts.data)

k600_crusius = k.crusius(ts.data)

ha2m2 <- 10000
kd        = sp.data$metadata$averagekd
wnd.z      = 10   #because we converted to u10
atm.press  = 1018
lat       = sp.data$metadata$latitude
lake.area = sp.data$metadata$lakearea*ha2m2

#for k.read and k.macIntyre, we need LW_net.
#Calculate from the observations we have available.
lwnet = calc.lw.net(ts.data, lat, atm.press)
ts.data = merge(ts.data, lwnet)

k600_read = k.read(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press,
									 lat=lat, lake.area=lake.area)

k600_soloviev = k.read.soloviev(ts.data, wnd.z=wnd.z, Kd=kd,
																atm.press=atm.press, lat=lat, lake.area=lake.area)

k600_macIntyre = k.macIntyre(ts.data, wnd.z=wnd.z, Kd=kd, atm.press=atm.press)

k600_heiskanen = k.heiskanen(ts.data, wnd.z, kd, atm.press)


k600_vachon = k.vachon(ts.data, lake.area)


# ---- figure ----
cols <- c("#1b9e77", "#d95f02", "black", "#e7298a", "DodgerBlue", "#e6ab02", "grey50")

models <- list(
  list('name'="MacIntyre", data = k600_macIntyre, col = cols[1], lty = 6, lwd = 1.7),
  list('name'="Cole", data = k600_cole, col = cols[2], lty = 1, lwd = 1.2),
  list('name'="Vachon", data = k600_vachon, col = cols[3], lty = 1, lwd = 1.1),
  list('name'="Read", data = k600_read, col = cols[4], lty = 1, lwd = 1.2),
  list('name'="Soloviev", data = k600_soloviev, col = cols[5], lty = 6, lwd = 1.7),
  list('name'="Heiskanen", data = k600_heiskanen, col = cols[6], lty = 1, lwd = 1.2),
  list('name'="Crusius", data = k600_crusius, col = cols[7], lty = 1, lwd = 1.2))
  
add_axes <- function(xlim, ylim, ylabel = pretty(ylim,10), panel.txt){
  prc_x = 0.1 # position for text relative to axes
  prc_y = 0.07
  tick_len <- 0.15
  ext_x <- c(xlim[1]-86400, pretty(xlim,3), xlim[2]+86400)
  ext_y <- c(ylim[1]-10, pretty(ylim,10), ylim[2]+10)
  ylab <- c("",ylabel,"")
  if (is.na(ylabel[1])) ylab = NA
  axis(side = 1, at = ext_x , labels = strftime(ext_x,'%H:%M'), tcl = tick_len)
  axis(side = 2, at = ext_y, labels = ylab, tcl = tick_len)
  axis(side = 3, at = ext_x, labels = NA, tcl = tick_len)
  axis(side = 4, at = ext_y, labels = NA, tcl = tick_len)
  x_txt <- (xlim[2] - xlim[1])*prc_x+xlim[1]
  y_txt <- ylim[2]-(ylim[2] - ylim[1])*prc_y
  text(x = x_txt, y_txt,labels = panel.txt)
}
add_night <- function(xlim, ylim){
  rise.set = sun.rise.set(xlim[1]+43200, lat) # add mid-day
  polygon(x = c(xlim[1], rise.set[1], rise.set[1], xlim[1]), y = c(ylim[1],ylim[1],ylim[2],ylim[2]), col = night_col,border = NA)
  polygon(x = c(xlim[2], rise.set[2], rise.set[2], xlim[2]), y = c(ylim[1],ylim[1],ylim[2],ylim[2]), col = night_col,border = NA)
}
moving_ave <- function(df, window = 18){
  out <- df[,2]*NA
  l_win <- floor(window/2)
  r_win <- window-l_win
  end_win <- nrow(df)-r_win
  for (i in l_win:end_win){
    strt <- i-l_win+1
    out[i] <- mean(df[strt:(strt+r_win), 2])
  }
  df[,2] <- out
  return(df)
}

add_models <- function(models){
  .empty = sapply(X = models, FUN = function(x) {
    lines(moving_ave(x$data), col=x$col, lty = x$lty, lwd = x$lwd)
    })
}

add_legend <- function(models, xlim, ylim, prc_x = 0.26, prc_y = 0.06){
  
  y_strt <- ylim[2]-(ylim[2] - ylim[1])*prc_y
  y_spc <- (ylim[2] - ylim[1])*0.05
  x_len <- (xlim[2] - xlim[1])*0.16
  x <- c((xlim[2] - xlim[1])*prc_x+xlim[1], (xlim[2] - xlim[1])*prc_x+xlim[1] + x_len)
  
  for (i in 1:length(models)){
    y = y_strt-(i-1)*y_spc
    lines(x, c(y,y), 
          col =models[[i]]$col, 
          lty = models[[i]]$lty, 
          lwd = models[[i]]$lwd)
    text(x[2],y, models[[i]]$name, pos = 4, cex = 0.65)
  }
}

width = 3.37 # single column width for journal
night_col = 'grey90'
height = 1.9
l_mar = 0.35
b_mar = 0.3
t_mar = 0.05
r_mar= 0.05
gapper = 0.15 # space between panels
ylim = c(0,7.55)
xlim = as.POSIXct(c('2009-07-02 22:00', '2009-07-04 02:00', '2009-07-08 22:00', '2009-07-10 02:00'))


#Create plot and save in user home directory (on Windows, Documents folder, on Mac, Home folder)
png('~/k600_figure.png', res=300, width=width, height=height, units = 'in')

#


layout(matrix(c(rep(1,10),rep(2,9)),nrow=1)) # 55% on the left panel
par(mai=c(b_mar,l_mar,t_mar,0), omi = c(0,0,0,r_mar),xpd=FALSE,
    mgp = c(1.15,.05,0))

plot(c(0,NA),c(0,NA), type='l', 
     axes = FALSE,
     xaxs = 'i', yaxs = 'i',
     ylim=ylim, 
     ylab=expression(k[600]~(m~day^-1)),
     xlab=strftime(mean(xlim[1:2]), '%Y-%m-%d'),
     xlim=xlim[1:2])
add_night(xlim[1:2], ylim)
add_models(models)
add_legend(models, xlim, ylim)
add_axes(xlim[1:2], ylim, panel.txt = 'a)')

par(mai=c(b_mar,gapper,t_mar,0))

plot(c(0,NA),c(0,NA), type='l', 
     axes = FALSE,
     xaxs = 'i', yaxs = 'i',
     ylim=ylim, 
     ylab=NA,
     xlab=strftime(mean(xlim[3:4]), '%Y-%m-%d'),
     xlim=xlim[3:4])
add_night(xlim[3:4], ylim)
add_models(models)
add_axes(xlim[3:4], ylim, ylabel = NA, panel.txt = 'b)')
#par(xpd=TRUE)

dev.off()

#quick summary of k600 values for results section

cat('Cole:', mean(k600_cole[,2]), '\n')
cat('Crusius:', mean(k600_crusius[,2]), '\n')
cat('MacIntyre:', mean(k600_macIntyre[,2]), '\n')
cat('Read:', mean(k600_read[,2]), '\n')
cat('Read_soloviev:', mean(k600_soloviev[,2]), '\n')
cat('heiskanen:', mean(k600_heiskanen[,2]), '\n')
cat('vachon:', mean(k600_vachon[,2]), '\n')


