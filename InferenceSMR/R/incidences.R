incidences <-
function(age_min, age_max, year_min, year_max, follow_up, start_age, start_year, case)
{
incid<-array(0,c(length(age_min:age_max),length(year_min:year_max),2));
for(i in 1:length(follow_up))
{
X = follow_up[i];
dep_t = start_year[i];
age_t = start_age[i];
while(dep_t - start_year[i] < X)
{
Y = min(floor(age_t) + 1 - age_t, floor(dep_t) + 1 - dep_t);
if(dep_t - start_year[i] + Y > X)
{
incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,1] = incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,1] + X - dep_t + start_year[i];
incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,2] = incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,2] + case[i];
}
else
{
incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,1] = incid[floor(age_t) - age_min + 1, floor(dep_t) - year_min + 1,1] + Y;
}
dep_t = dep_t + Y;
age_t = age_t + Y;
}
}
incidences2<-incid[,,2]/incid[,,1];
incidences2[incidences2=='NaN'] = 0;

pre_incid<-matrix(0,nrow=length(seq(age_min,age_max))*length(seq(year_min,year_max)),ncol=4)
k=1;
for(i in seq(age_min,age_max))
{
for(j in seq(year_min,year_max))
{
pre_incid[k,]<-c(incidences2[i - age_min + 1, j - year_min + 1], j, i, incid[i - age_min + 1, j - year_min + 1, 1]);
k = k + 1;
}
}
return(pre_incid);
}
