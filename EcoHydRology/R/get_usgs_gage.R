get_usgs_gage<-function(flowgage_id,begin_date="1979-01-01",end_date="2013-01-01"){
#
# Grabs USGS stream flow data for 1979 to present... to align with CFSR datasets.
#

url = paste("http://waterdata.usgs.gov/nwis/inventory?search_site_no=", flowgage_id, "&search_site_no_match_type=exact&sort_key=site_no&group_key=NONE&format=sitefile_output&sitefile_output_format=rdb&column_name=station_nm&column_name=site_tp_cd&column_name=dec_lat_va&column_name=dec_long_va&column_name=alt_va&column_name=drain_area_va&column_name=contrib_drain_area_va&column_name=rt_bol&list_of_search_criteria=search_site_no",sep="")

gage_tsv=readLines(url)
gage_tsv=gage_tsv[grep("^#",gage_tsv,invert=T)][c(1,3)]
tempdf=read.delim(text=gage_tsv,sep="\t",header=T,colClasses = c("character", "character", "numeric", "numeric", "character", "character", "numeric", "numeric", "character", "numeric"))
area = tempdf$drain_area_va* 1.6^2
if(is.na(area)) {area=0;print("warning, no area associated with gage, setting to 0\n")}
declat = tempdf$dec_lat_va
declon = tempdf$dec_long_va
elev = tempdf$alt_va* 12/25.4
if(is.na(elev)) {elev=0;print("warning, no elevation associated with gage, setting to 0\n")}

gagename = tempdf$station_nm
begin_date = as.character(begin_date)
end_date = as.character(end_date)

url = paste("http://nwis.waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&begin_date=",begin_date,"&end_date=",end_date,"&site_no=", flowgage_id, sep = "")
flowdata_tsv=gage_tsv=readLines(url)
flowdata_tsv=flowdata_tsv[grep("^#",flowdata_tsv,invert=T)][c(3:length(flowdata_tsv))]
flowdata = read.delim(text=flowdata_tsv,header=F,sep="\t",col.names = c("agency", "site_no", "date", "flow", "quality"), colClasses = c("character", "numeric", "character", "character", "character"), fill = T)
flowdata$mdate = as.Date(flowdata$date, format = "%Y-%m-%d")
flowdata$flow = as.numeric(as.character(flowdata$flow)) * 12^3 * 2.54^3/100^3 * 24 * 3600
flowdata = na.omit(flowdata)
returnlist = list(declat = declat, declon = declon, flowdata = flowdata, area = area, elev = elev, gagename = gagename)
return(returnlist)

}


