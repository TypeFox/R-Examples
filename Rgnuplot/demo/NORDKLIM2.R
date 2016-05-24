
if (!(file.exists('NORDKLIM-Helsinki-prec-columns.dat'))) stop('Please run demo NORDKLIM1 first, to create the data file')

lbls<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

#Initialize the gnuplot handle
h1<-Gpinit()
#Gpcmd(h1,'unset key; set grid xtics; set grid ytics; set ztics scale 5; set xmtics; 
#set pm3d ; set xrange[1:12]\nsplot \"NORDKLIM304-601FIN3.dat\" using 1:2:3  w l  linewidth 2 linecolor 3' )

Gpcmd(h1,'unset key; set grid xtics; set grid ytics; set ztics scale 10;set pm3d ;set xrange[1:12];
set zrange [0:2000]
set ticslevel 0 #place z-axis at the xy level
set size ratio 1,.1
set xtics ("Jan" 1,"Feb" 2,"Mar" 3,"Apr" 4,"May" 5,"Jun" 6,"Jul" 7,"Aug" 8,"Sep" 9,"Oct" 10,"Nov" 11,"Dec" 12)
splot \"NORDKLIM-Helsinki-prec-columns.dat\" using 1:2:3 w l  linewidth 2 linecolor 3' )

if (interactive()) {
   plot.hist <- function(panel) {
         with(panel, {
scounter=0
sxtics=''
sSelect=''
         for (n in 1:12)
         if (panel$cbox[n])
         {
         scounter<-scounter+1
sxtics<-sxtics %s% '\"' %s% lbls[n] %s% '\" ' %s% scounter
if (n<12) sxtics<-sxtics %s% ','
sSelect<-sSelect %s% '($1==' %s% n %s% ')?' %s% scounter %s% ':'
         if (n==12) sSelect<-sSelect %s% '0'
         }
         else
         {
       sSelect<-sSelect %s% '($1==' %s% n %s% ')?0:'
         if (n==12) sSelect<-sSelect %s% '0'
         }

notsure<<-'unset key; set grid xtics; set grid ytics; set ztics scale 5;set pm3d ;set xrange[1:' %s% scounter %s% '];
set xtics (' %s% sxtics %s% ')
splot \"NORDKLIM-Helsinki-prec-columns.dat\" using (' %s% sSelect %s% '):2:3 w l  linewidth 2 linecolor 3'
         
Gpresetplot(h1)
Gpcmd(h1,'unset key; set grid xtics; set grid ytics; set ztics scale 10;set pm3d ;set xrange[1:' %s% scounter %s% '];
set zrange [0:2000]
set ticslevel 0 #place z-axis at the xy level
set size ratio 1,.1
set xtics (' %s% sxtics %s% ')
splot \"NORDKLIM-Helsinki-prec-columns.dat\" using (' %s% sSelect %s% '):2:3 w l  linewidth 2 linecolor 3' )

#\"Jan\" 1,\"Feb\" 2,\"Mar\" 3,\"Apr\" 4,\"May\" 5,\"Jun\" 6,\"Aug\" 8,\"Sep\" 9,\"Oct\" 10,\"Nov\" 11,\"Dec\" 12
#( ($1==1)?1:($1==2)?2:($1==3)?3:($1==4)?0:($1==5)?4:($1==6)?5:($1==7)?6:($1==8)?7:($1==9)?0:($1==10)?8:($1==11)?9:($1==12)?10:0 )
#which(panel$labels[n]==lbls)
#set xtics (1,2,3,4,5,6,8,9,10,11,12);
#set xtics ("April" 1, "May" 2, "June" 3, "July" 4),'Jul' 7

           })
         panel
         }
   panel <- rp.control()
   rp.checkbox(panel, cbox, plot.hist, initval =c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
      labels = lbls, title = "Months")
   rp.do(panel, plot.hist)#
   }

#Gppause()
#close gnuplot handles
#h1<-Gpclose(h1)    
