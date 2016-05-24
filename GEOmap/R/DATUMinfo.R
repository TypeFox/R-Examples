DATUMinfo<-function()
   {
  dnames = c('Datum', 'Equat Radius, m(a)', 'Polar Radius, m(b)','Flat a:(a-b)/a', 'Use')   
  HEY =  c('NAD83/WGS84',6378137, 6356752.3142, 298.257223563 , 'Global',
'GRS 80 ', 6378137, 6356752.3141,298.257222101,'US',
'WGS72' , 6378135,6356750.5 , 298.26 ,'NASA DOD',
'Australian 1965',6378160 , 6356774.7,298.25 , 'Australia',
'Krasovsky 1940',6378245 , 6356863.0 , 298.3 , 'Soviet Union',
'International (1924) -Hayford (1909)',6378388,6356911.9 , 297, 'Global except as listed',
'Clake 1880' ,	6378249.1 , 6356514.9 ,	293.46 , 'France Africa',
'Clarke 1866',6378206.4 , 6356583.8 , 294.98 , 'North America',
'Airy 1830' ,	6377563.4 , 6356256.9 , 299.32 , 'Great Britain',
'Bessel 1841',6377397.2 , 6356079.0 ,	299.15 , 'Central Europe Chile Indonesia',
 'Everest 1830', 6377276.3 , 6356075.4 , 300.80 ,  'South Asia')

h = vector(mode='list')

   for(i in 1:length(dnames))
   {
h[[i]] = HEY[seq(from=i, by=length(dnames), to=length(HEY))  ]

}
h[[2]] = as.numeric(h[[2]])
h[[3]] = as.numeric(h[[3]])
h[[4]] = as.numeric(h[[4]])
   
  names(h) = dnames
return(h)
}

