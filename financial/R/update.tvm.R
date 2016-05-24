"update.tvm" <-
function (object,i = NULL,n = NULL,pv = NULL,fv = NULL,pmt = NULL,days = NULL,adv = NULL,pyr = NULL,cyr = NULL,...) 
{

options(warn=-1);

	x=object;

if (is.null(i)) i = x[,1]
if (is.null(n)) n = x[,2]
if (is.null(pv)) pv = x[,3]
if (is.null(fv)) fv = x[,4]
if (is.null(pmt)) pmt = x[,5]
if (is.null(days)) days = x[,6]
if (is.null(adv)) adv = x[,7]
if (is.null(pyr)) pyr = x[,8]
if (is.null(cyr)) cyr = x[,9]

ii = i;
nn = n;

res =c()

l = max(c(length(i),length(n),length(pv),length(fv),length(pmt),length(days),
length(adv),length(pyr),length(cyr)));

ii = ii+rep(0,l);
nn = nn+rep(0,l);
pv = pv+rep(0,l);
fv = fv+rep(0,l);
pmt = pmt+rep(0,l);
days = days+rep(0,l);
adv = adv+rep(0,l);
pyr = pyr+rep(0,l);
cyr = cyr+rep(0,l);

for (k in 1:l) {

b=pv[k]; p=pmt[k]; f=fv[k]; d=days[k]; 
a=adv[k]; i=ii[k]; n=nn[k]; py=pyr[k]; cy=cyr[k];

if (py!=cy) i = irnom(ireff(i,cy),py);

i=i/(100*py);

s=1;

q=d/(360/py);


if (sum(is.na(c(b,p,f,i,n,a)))!=1)
{  stop("Incorrect number of Not Available values"); }

if (is.na(b)) b = (i + 1)^(-n)*(p*(i + 1)^a*(i*s + 1) - p*(i + 1)^n*(a*i + 1)*(i*s + 1) - f*i)/(i*(i*q + 1));

if (is.na(p)) p = i*(b*(i + 1)^n*(i*q + 1) + f)/(((i + 1)^a - (i + 1)^n*(a*i + 1))*(i*s + 1));

if (is.na(f)) f = (p*(i + 1)^(2*a)*(i*s + 1) - (i + 1)^(a + n)*(2*a*i*p*(i*s + 1) + b*i*(i*q + 1) + 2*p*(i*s + 1))
  + (i + 1)^(2*n)*(a*i + 1)*(a*i*p*(i*s + 1) + b*i*(i*q + 1) + p*(i*s + 1)))/(i*((i + 1)^a - (i + 1)^n*(a*i + 1)));

if (is.na(n)) n = log((p*(i + 1)^a*(i*s + 1) - f*i)/(a*i*p*(i*s + 1) + b*i*(i*q + 1) + p*(i*s + 1)))/log(i + 1)

if (is.na(i)) {
fi = function (i) {
 (- b*(1 + i*q) - f*(1 + i)^(-n))/((1 - (1 + i)^(- (n - a)))/i + a) - p*(1 + i*s)
}
r = try(uniroot(fi,c(1e-10,1/py),tol=1e-10),silent=T);
if (inherits(r,"try-error")) i = NA else i = r$root;

}
if (is.na(a)) {
fa = function (a) {
 (- b*(1 + i*q) - f*(1 + i)^(-n))/((1 - (1 + i)^(- (n - a)))/i + a) - p*(1 + i*s)
}
r = try(uniroot(fa,c(0,n),tol=1e-10),silent=T);
if (inherits(r,"try-error")) a = NA else a = r$root;
}
i=i*py*100;
res=rbind(res,c(i,n,b,f,p,d,a,py,cy));

}

class(res)="tvm";
rownames(res)=1:k;
colnames(res)=c("I%","#N","PV","FV","PMT","Days","#Adv","P/YR","C/YR");
return(res);

}

