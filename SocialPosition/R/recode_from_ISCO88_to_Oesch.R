recode_from_ISCO88_to_Oesch <-
function (ISCO88, EMP_STA, SE_zero_emp, SE_one_to_nine_emp, SE_ten_plus_emp, not_SE, data) 
{
#creation variable number of employee for self-employed
data$emplnum[(EMP_STA %in% c(SE_zero_emp))] <-0
data$emplnum[(EMP_STA %in% c(SE_one_to_nine_emp))] <-1
data$emplnum[(EMP_STA %in% c(SE_ten_plus_emp))] <-2
data$emplnum[(EMP_STA %in% c(not_SE))] <-NA
##self-employed / independent work logic
data$Oesch17[(data$emplnum %in% c(2))] <- 1 #large employers (>9) //// empstat = 1
data$Oesch17[((data$emplnum %in% c(0:1)) & (ISCO88 %in% c(2,20:24,200:247,2000:2479)))] <- 2 #self-employed professionals
data$Oesch17[(data$emplnum %in% c(1)) & (ISCO88 %in% c(1, 10:19, 100:199, 1000:1999, 3:9, 30:93, 300:939, 3000:9399))] <- 3 #small proprietors, artisans, with employees (<9) //// empstat = 2 & is3maj99 not in 200:247 /// 2,20:24,200:247,2000:2479
data$Oesch17[(data$emplnum %in% c(0)) & (ISCO88 %in% c(1, 10:19, 100:199, 1000:1999, 3:9, 30:93, 300:939, 3000:9399))] <- 4 #small proprietirs, artisans, without employees //// empstat = 3 and not in 200:247
##employees / technical work logic
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(2,21,200:214,221,2100:2213))] <- 5 #technical experts 
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(3,31,300:315,321,3100:3213,3434))] <- 6 #technicians
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(60,61,600:615,70:74,700:744,83,831,834,6100:7442,8311,8323:8324,8332:8340))] <- 7 #skilled crafts
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(80:82,800:829,832,93,916,930:933,8000:8290,8312:8322,9153:9162,9300:9330))] <- 8 #routine operatives
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(92,921,833,8330:8331,9211:9213))] <- 9 #routine agriculture
##employees / organizational work logic
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(10:11,12,100:111,114:123,241:242,247,1000:1110,1141:1251,2400:2429,2441,2470))] <- 10 #higher-grade managers
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(13,131,34,340:345,1310:1319,3400:3433,3441:3450))] <- 11 #associate managers
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(40:42,400:422,516,4000:4112,4114:4190,4200:4210,4212:4222,5161:5169))] <- 12 #skilled office
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(4113,4211,4223))] <- 13 #routine office
##employees / interpersonal work logic
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(22:23,222:232,235,244,246,2221:2320,2351:2359,2442,2443,2445:2451,2460))] <- 14 #socio-cultural professionals
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(233:234,243,245,32:33,322:334,346:348,2331:2340,2431,2432,2444,2452:2455,3221:3224,3226,3229:3340,3460:3472,3480,3000))] <- 15 #socio-cultural semi-professionals
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(51:52,500,511,521:522,3225,3227,3228,3473:3475,5111:5113,5122:5141,5143,5149,5210,5220))] <- 16 #skilled service
data$Oesch17[(data$emplnum %in% NA) & (ISCO88 %in% c(512:514,91,900:915,5120:5121,5142,9100:9152))] <- 17 #routine service
##grouping
data$Oesch8[data$Oesch17 %in% c(1,2)] <- 1
data$Oesch8[data$Oesch17 %in% c(3,4)] <- 2
data$Oesch8[data$Oesch17 %in% c(5,6)] <- 3
data$Oesch8[data$Oesch17 %in% c(7,8,9)] <- 4
data$Oesch8[data$Oesch17 %in% c(10,11)] <- 5
data$Oesch8[data$Oesch17 %in% c(12,13)] <- 6
data$Oesch8[data$Oesch17 %in% c(14,15)] <- 7
data$Oesch8[data$Oesch17 %in% c(16,17)] <- 8
return(data)
}
