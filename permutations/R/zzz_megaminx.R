jj <- array(c(
    010,012,014,016,018,  011,013,015,017,019,  021,033,045,057,069,  022,034,046,058,060,  023,035,047,059,061,
    020,022,024,026,028,  021,023,025,027,029,  015,067,091,081,035,  016,068,092,082,036,  017,069,093,083,037,
    030,032,034,036,038,  031,033,035,037,039,  017,029,089,079,047,  018,020,080,070,048,  019,021,081,071,049,
    040,042,044,046,048,  041,043,045,047,049,  010,032,078,118,050,  011,033,079,119,051,  019,031,077,117,059,
    050,052,054,056,058,  051,053,055,057,059,  011,043,115,105,061,  012,044,116,106,062,  013,045,117,107,063,
    060,062,064,066,068,  061,063,065,067,069,  013,055,103,093,023,  014,056,104,094,024,  015,057,105,095,025,
    070,072,074,076,078,  071,073,075,077,079,  030,088,120,110,040,  031,089,121,111,041,  039,087,129,119,049,
    080,082,084,086,088,  081,083,085,087,089,  027,099,121,071,037,  028,090,122,072,038,  029,091,123,073,039,
    090,092,094,096,098,  091,093,095,097,099,  025,065,101,123,083,  026,066,102,124,084,  027,067,103,125,085,
    100,102,104,106,108,  101,103,105,107,109,  053,113,125,095,063,  054,114,126,096,064,  055,115,127,097,065,
    110,112,114,116,118,  111,113,115,117,119,  041,075,127,107,051,  042,076,128,108,052,  043,077,129,109,053,
    120,122,124,126,128,  121,123,125,127,129,  073,085,097,109,111,  074,086,098,100,112,  075,087,099,101,113
    ),c(5,5,12))

megaminx <- rep(id,12)

for(i in seq_len(12)){
    for(j in seq_len(5)){
        megaminx[i] <- megaminx[i] + as.cycle(jj[,j,i])
    }
}
rm(i,j,jj)


names(megaminx) <- c("White", "Purple", "DarkYellow", "DarkBlue", "Red",
                   "DarkGreen", "LightGreen", "Orange", "LightBlue",
                   "LightYellow", "Pink", "Grey")



"W"  <- megaminx[01]  # "White" 
"Pu" <- megaminx[02]  # "Purple" 
"DY" <- megaminx[03]  # "DarkYellow" 
"DB" <- megaminx[04]  # "DarkBlue" 
"R"  <- megaminx[05]  # "Red" 
"DG" <- megaminx[06]  # "DarkGreen" 
"LG" <- megaminx[07]  # "LightGreen" 
"O"  <- megaminx[08]  # "Orange" 
"LB" <- megaminx[09]  # "LightBlue" 
"LY" <- megaminx[10]  # "LightYellow" 
"Pi" <- megaminx[11]  # "Pink" 
"Gy" <- megaminx[12]  # "Grey"; this cannot be 'Gr' as confused with Green. 

megaminx_colours <- noquote(rep("Black",126))
megaminx_colours[010:019] <- "White"
megaminx_colours[020:029] <- "Purple"
megaminx_colours[030:039] <- "DarkYellow"
megaminx_colours[040:049] <- "DarkBlue"
megaminx_colours[050:059] <- "Red"
megaminx_colours[060:069] <- "DarkGreen"
megaminx_colours[070:079] <- "LightGreen"
megaminx_colours[080:089] <- "Orange"
megaminx_colours[090:099] <- "LightBlue"
megaminx_colours[100:109] <- "LightYellow"
megaminx_colours[110:119] <- "Pink"
megaminx_colours[110:129] <- "Grey"
