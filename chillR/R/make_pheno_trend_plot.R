#' Make image showing phenology response to temperatures during two phases
#' 
#' The timing of many developmental stages of temperate plants is understood to
#' depend on temperatures during two phases. This function seeks to illustrate
#' this dependency by plotting phenological dates as colored surface, as a
#' function of mean temperatures during both phases, which are indicated on the
#' x and y axes.
#' 
#' The generation of the color surface is based on the Kriging technique, which
#' is typically used for interpolation of spatial data. The use for this
#' particular purpose is a bit experimental. The function uses the Krig
#' function from the fields package.
#' 
#' @param weather_data_frame a dataframe containing daily minimum and maximum
#' temperature data (in columns called Tmin and Tmax, respectively), and/or
#' mean daily temperature (in a column called Tmean). There also has to be a
#' column for Year and one for JDay (the Julian date, or day of the year).
#' Alternatively, the date can also be given in three columns (Year, Month and
#' Day).
#' @param split_month the procedure analyzes data by phenological year, which
#' can start and end in any month during the calendar year (currently only at
#' the beginning of a month). This variable indicates the last month (e.g. 5
#' for May) that should be included in the record for a given phenological
#' year. All subsequent months are assigned to the following phenological year.
#' @param pheno a data frame that contains information on the timing of
#' phenological events by year. It should consist of two columns called Year
#' and pheno. Data in the pheno column should be in Julian date (day of the
#' year).
#' @param use_Tmean boolean variable indicating whether or not the column Tmean
#' from the weather_data_frame should be used as input for the PLS analysis. If
#' this is set to FALSE, Tmean is calculated as the arithmetic mean of Tmin and
#' Tmax.
#' @param Start_JDay_chill the Julian date, on which the first relevant period
#' (e.g. the chilling phase) starts
#' @param End_JDay_chill the Julian date, on which the first relevant period
#' (e.g. the chilling phase) ends
#' @param Start_JDay_heat the Julian date, on which the second relevant period
#' (e.g. the forcing phase) starts
#' @param End_JDay_heat the Julian date, on which the second relevant period
#' (e.g. the forcing phase) ends
#' @param outpath the output path
#' @param file_name the output file name
#' @param plot_title the title of the plot
#' @param ylabel the label for the y-axis. There is a default, but in many
#' cases, it may be desirable to customize this
#' @param xlabel the label for the x-axis. There is a default, but in many
#' cases, it may be desirable to customize this
#' @param legend_label the label for the legend (color scheme). There is a
#' default, but in many cases, it may be desirable to customize this
#' @param image_type the type of image to produce. This currently has only two
#' options: "tiff" or anything else (the default). If this is not "tiff", a png
#' image is produced. The "tiff" option was added to produce publishable
#' figures that adhere to the requirements of most scientific journals.
#' @param colorscheme the color scheme for the figure. This currently has only
#' two options: "bw" or anythings else (the default). "bw" produces a grayscale
#' image, otherwise the figure will be in color
#' @param fonttype font style to be used for the figure. Can be 'serif'
#' (default) or 'sans'.
#' @return \item{pheno}{data frame containing all data needed for reproducing
#' the plot: Year (during which the phenological event occurred - the year in
#' which the phenological season indicated by split_month ended), pheno (the
#' date on which the phenological event occurred), Chill_Tmean (mean
#' temperature during the first relevant phase), Heat_Tmean (mean temperature
#' during the second relevant phase) and Year_Tmean (mean annual temperature -
#' not actually used in the plot)} \item{ylabel}{character string used for
#' labeling the y axis} \item{xlabel}{character string used for labeling the x
#' axis}
#' @author Eike Luedeling
#' @references Guo L, Dai J, Wang M, Xu J, Luedeling E, 2015. Responses of
#' spring phenology in temperate zone trees to climate warming: a case study of
#' apricot flowering in China. Agricultural and Forest Meteorology 201, 1-7.
#' 
#' Guo L, Dai J, Ranjitkar S, Xu J, Luedeling E, 2013. Response of chestnut
#' phenology in China to climate variation and change. Agricultural and Forest
#' Meteorology 180, 164-172.
#' 
#' Luedeling E, Guo L, Dai J, Leslie C, Blanke M, 2013. Differential responses
#' of trees to temperature variation during the chilling and forcing phases.
#' Agricultural and Forest Meteorology 181, 33-42.
#' 
#' the interpolation was done according to:
#' 
#' Furrer, R., Nychka, D. and Sain, S., 2012. Fields: Tools for spatial data. R
#' package version 6.7.
#' 
#' Reference to the chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords phenology analysis
#' @examples
#' 
#' weather<-fix_weather(KA_weather)
#' 
#' #the output of the PLS function (PLS_pheno, plotted with plot_PLS) can be used to select
#' #phases that are likely relevant for plant phase timing. See respective examples for running
#' #these functions.
#' 
#' file_path<-paste(getwd(),"/",sep="")
#' 
#' mpt<-make_pheno_trend_plot(weather_data_frame = weather$weather, split_month = 6,
#'      pheno = KA_bloom, use_Tmean = FALSE, Start_JDay_chill = 260, 
#'      End_JDay_chill = 64, Start_JDay_heat = 44, End_JDay_heat = 103,
#'      outpath = file_path, file_name = "pheno_trend_plot",
#'      plot_title = "Impacts of chilling and forcing temperatures on cherry phenology",
#'      ylabel = NA, xlabel = NA, legend_label = NA, image_type = "png", 
#'      colorscheme = "normal")
#' 
#' @export make_pheno_trend_plot
make_pheno_trend_plot <-
function(weather_data_frame,
                                split_month=6,   #last month in same year
                                pheno,
                                use_Tmean=FALSE,
                                Start_JDay_chill,
                                End_JDay_chill,
                                Start_JDay_heat,
                                End_JDay_heat,
                                outpath,
                                file_name,
                                plot_title,
                                ylabel=NA,
                                xlabel=NA,
                                legend_label=NA,
                                image_type="png",
                                colorscheme="normal",
                                fonttype="serif")

{

  
weather_file<-weather_data_frame
weather_file$Tmax<-suppressWarnings(as.numeric(as.character(weather_file$Tmax)))
weather_file$Tmin<-suppressWarnings(as.numeric(as.character(weather_file$Tmin)))
suppressWarnings(weather_file[which(weather_file$Tmin>weather_file$Tmax),c("Tmin","Tmax")]<-NA)

Tmin_gaps<-interpolate_gaps(weather_file$Tmin)
weather_file$Tmin<-Tmin_gaps[[1]]
weather_file[,"Tmin_interpolated"]<-0
suppressWarnings(weather_file[Tmin_gaps[[2]],"Tmin_interpolated"]<-1)
Tmax_gaps<-interpolate_gaps(weather_file$Tmax)
weather_file$Tmax<-Tmax_gaps[[1]]
weather_file[,"Tmax_interpolated"]<-0
suppressWarnings(weather_file[Tmax_gaps[[2]],"Tmax_interpolated"]<-1)

if(use_Tmean)
{Tmean_gaps<-interpolate_gaps(weather_file$Tmean)
 weather_file$Tmean<-Tmean_gaps[[1]]
 weather_file[,"Tmean_interpolated"]<-0
 weather_file[Tmean_gaps[[2]],"Tmean_interpolated"]<-1} else
   
 {weather_file[,"Tmean"]<-(weather_file$Tmax+weather_file$Tmin)/2}

weather_file[which(weather_file$Month<=split_month),"Season"]<-weather_file[which(weather_file$Month<=split_month),"Year"]
weather_file[which(weather_file$Month>split_month),"Season"]<-weather_file[which(weather_file$Month>split_month),"Year"]+1
weather_file[,"Date"]<-weather_file$Month*100+weather_file$Day
weather_file[,"JDay"]<-strptime(paste(weather_file$Month,"/",weather_file$Day,"/",weather_file$Year,sep=""),"%m/%d/%Y")$yday+1

#running mean

ww<-weather_file[,"Tmean"]
rr<-weather_file[,"Tmean"]        

sea<-unique(weather_file$Season)
res<-data.frame(Season=sea,Chill_Tmean=NA,Heat_Tmean=NA,Year_Tmean=NA)

if(Start_JDay_chill>End_JDay_chill) chill_days<-c(Start_JDay_chill:366,1:End_JDay_chill) else
  chill_days<-Start_JDay_chill:End_JDay_chill
  

if(Start_JDay_heat>End_JDay_heat) heat_days<-c(Start_JDay_heat:366,1:End_JDay_heat) else
  heat_days<-Start_JDay_heat:End_JDay_heat


for (s in sea)
{res[which(res$Season==s),"Chill_Tmean"]<-mean(weather_file[which(weather_file$Season==s&weather_file$JDay %in% chill_days),"Tmean"])
 res[which(res$Season==s),"Heat_Tmean"]<-mean(weather_file[which(weather_file$Season==s&weather_file$JDay %in% heat_days),"Tmean"])
 res[which(res$Season==s),"Year_Tmean"]<-mean(weather_file[which(weather_file$Season==s),"Tmean"])
}


for (i in 1:nrow(pheno))
  {if(pheno[i,1] %in% res$Season)
   {pheno[i,"Chill_Tmean"]<-res[which(res$Season==pheno[i,1]),"Chill_Tmean"]
   pheno[i,"Heat_Tmean"]<-res[which(res$Season==pheno[i,1]),"Heat_Tmean"]
  pheno[i,"Year_Tmean"]<-res[which(res$Season==pheno[i,1]),"Year_Tmean"]
  }}


Make_date<-function(Jday)
{
  leg<-strptime(strptime(paste(Jday,"/2001",sep=""),format="%j/%Y"),"%Y-%m-%d")
 comps<-strsplit(as.character(leg),"-")
  JM<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[as.numeric(comps[[1]][2])]
  return(paste(comps[[1]][3],JM))
}


if (is.na(xlabel)) xlabel<-paste("Mean temperature during the chilling phase (",Make_date(Start_JDay_chill)," to ",Make_date(End_JDay_chill),"; deg. C)",sep="")
if (is.na(ylabel)) ylabel<-paste("Mean temperature during the forcing phase (",Make_date(Start_JDay_heat)," to ",Make_date(End_JDay_heat),"; deg. C)",sep="")


pheno[,2]<-as.numeric(as.character(pheno[,2]))
pheno<-pheno[which(!is.na(pheno[,"pheno"])),]
pheno<-pheno[which(!is.na(pheno[,"Chill_Tmean"])),]
pheno<-pheno[which(!is.na(pheno[,"Heat_Tmean"])),]
pheno<-pheno[which(!is.na(pheno[,"Year_Tmean"])),]

k<-Krig(x=as.matrix(pheno[,c("Chill_Tmean","Heat_Tmean")]),Y=pheno$pheno)
#plot(k)

if (is.na(legend_label)) legend_label<-"Flowering date (Julian Day)"

if(image_type=="tiff") tiff(paste(outpath,file_name,".tif",sep=""),width=1000,height=1000) else
png(paste(outpath,file_name,".png",sep=""),width=1000,height=1000)
par(family=fonttype)
par(list(oma=c(0,2,0,2.2),mar=c(5.1,5.1,4.1,2.1)))
surface( k, type="C",xlab=xlabel,ylab=ylabel,cex.lab=2,cex.axis=1.5,labcex=1.5,asp=1,axis.args=list(cex.axis=2),legend.args=list(text=legend_label,side=4,cex=2,line=4.5))

if (colorscheme=="bw") colors<-colorRampPalette(c("#F9F9F9","#343434"))(255) else
  colors<-tim.colors(255)

surface( k, col=colors,type="C",xlab=xlabel,ylab=ylabel,cex.lab=2,cex.axis=1.5,labcex=1.5,asp=1,axis.args=list(cex.axis=2),legend.args=list(text=legend_label,side=4,cex=2,line=4.5))
         
mtext(text=plot_title,side=3,cex=2.5,line=2.3)
points(pheno[,c("Chill_Tmean","Heat_Tmean")],pch=16)
#if (is.na(legend_label)) mtext("Flowering date (Julian Day)",side=4,cex=2,line=6) else
#  mtext(legend_label,side=4,cex=2,line=6)
dev.off()

pheno<-pheno[,which(!is.na(pheno[1,]))]
write.csv(pheno,paste(outpath,file_name,".csv",sep=""),row.names=FALSE)
return(list(pheno=pheno,ylabel=ylabel,xlabel=xlabel))
}
