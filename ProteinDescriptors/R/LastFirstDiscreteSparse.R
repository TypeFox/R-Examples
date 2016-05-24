#' Last first discrete sparse descriptor.
#'
#' \code{LastFirstDiscreteSparse} returns the concatenation of the sum of sparse
#' descriptors of amino acids in a protein sequence and the sum of combination
#' of last part of an amino acid descriptor with the first part of its neighbour
#' amino acid descriptor.
#'
#' @param x A string of amino acid letters
#' @return A 40 dimensional numeric vector
#'
#' @export LastFirstDiscreteSparse
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' LastFirstDiscreteSparse(x)


LastFirstDiscreteSparse<-function (x)
{
 AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x1=strsplit(x, split = "")
result1=summary(factor(x1[[1]],levels=AAs),maxsum=21)/n

    sparse = t(data.frame(
	A = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	R = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	N = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	D = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	C = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	E = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	Q = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
      G = c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
	H = c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
	I = c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
	L = c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
	K = c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
	M = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
	F = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
      P = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
	S = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
	T = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
	W = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
	Y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
	V = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)))


DCDict = c("AA", "RA", "NA", "DA", "CA", "EA", "QA", "GA",
        "HA", "IA", "LA", "KA", "MA", "FA", "PA", "SA", "TA",
        "WA", "YA", "VA", "AR", "RR", "NR", "DR", "CR", "ER",
        "QR", "GR", "HR", "IR", "LR", "KR", "MR", "FR", "PR",
        "SR", "TR", "WR", "YR", "VR", "AN", "RN", "NN", "DN",
        "CN", "EN", "QN", "GN", "HN", "IN", "LN", "KN", "MN",
        "FN", "PN", "SN", "TN", "WN", "YN", "VN", "AD", "RD",
        "ND", "DD", "CD", "ED", "QD", "GD", "HD", "ID", "LD",
        "KD", "MD", "FD", "PD", "SD", "TD", "WD", "YD", "VD",
        "AC", "RC", "NC", "DC", "CC", "EC", "QC", "GC", "HC",
        "IC", "LC", "KC", "MC", "FC", "PC", "SC", "TC", "WC",
        "YC", "VC", "AE", "RE", "NE", "DE", "CE", "EE", "QE",
        "GE", "HE", "IE", "LE", "KE", "ME", "FE", "PE", "SE",
        "TE", "WE", "YE", "VE", "AQ", "RQ", "NQ", "DQ", "CQ",
        "EQ", "QQ", "GQ", "HQ", "IQ", "LQ", "KQ", "MQ", "FQ",
        "PQ", "SQ", "TQ", "WQ", "YQ", "VQ", "AG", "RG", "NG",
        "DG", "CG", "EG", "QG", "GG", "HG", "IG", "LG", "KG",
        "MG", "FG", "PG", "SG", "TG", "WG", "YG", "VG", "AH",
        "RH", "NH", "DH", "CH", "EH", "QH", "GH", "HH", "IH",
        "LH", "KH", "MH", "FH", "PH", "SH", "TH", "WH", "YH",
        "VH", "AI", "RI", "NI", "DI", "CI", "EI", "QI", "GI",
        "HI", "II", "LI", "KI", "MI", "FI", "PI", "SI", "TI",
        "WI", "YI", "VI", "AL", "RL", "NL", "DL", "CL", "EL",
        "QL", "GL", "HL", "IL", "LL", "KL", "ML", "FL", "PL",
        "SL", "TL", "WL", "YL", "VL", "AK", "RK", "NK", "DK",
        "CK", "EK", "QK", "GK", "HK", "IK", "LK", "KK", "MK",
        "FK", "PK", "SK", "TK", "WK", "YK", "VK", "AM", "RM",
        "NM", "DM", "CM", "EM", "QM", "GM", "HM", "IM", "LM",
        "KM", "MM", "FM", "PM", "SM", "TM", "WM", "YM", "VM",
        "AF", "RF", "NF", "DF", "CF", "EF", "QF", "GF", "HF",
        "IF", "LF", "KF", "MF", "FF", "PF", "SF", "TF", "WF",
        "YF", "VF", "AP", "RP", "NP", "DP", "CP", "EP", "QP",
        "GP", "HP", "IP", "LP", "KP", "MP", "FP", "PP", "SP",
        "TP", "WP", "YP", "VP", "AS", "RS", "NS", "DS", "CS",
        "ES", "QS", "GS", "HS", "IS", "LS", "KS", "MS", "FS",
        "PS", "SS", "TS", "WS", "YS", "VS", "AT", "RT", "NT",
        "DT", "CT", "ET", "QT", "GT", "HT", "IT", "LT", "KT",
        "MT", "FT", "PT", "ST", "TT", "WT", "YT", "VT", "AW",
        "RW", "NW", "DW", "CW", "EW", "QW", "GW", "HW", "IW",
        "LW", "KW", "MW", "FW", "PW", "SW", "TW", "WW", "YW",
        "VW", "AY", "RY", "NY", "DY", "CY", "EY", "QY", "GY",
        "HY", "IY", "LY", "KY", "MY", "FY", "PY", "SY", "TY",
        "WY", "YY", "VY", "AV", "RV", "NV", "DV", "CV", "EV",
        "QV", "GV", "HV", "IV", "LV", "KV", "MV", "FV", "PV",
        "SV", "TV", "WV", "YV", "VV")

   x2 = strsplit(x, split = "")[[1]]
    DC = summary(factor(paste(x2[-n], x2[-1], sep = ""),
        levels = DCDict), maxsum = 401)

    disparse = t(data.frame(
AA = c(sparse[1,11:20],sparse[1,1:10]),
RA = c(sparse[2,11:20],sparse[1,1:10]),
NA1 = c(sparse[3,11:20],sparse[1,1:10]),
DA = c(sparse[4,11:20],sparse[1,1:10]),
CA=c(sparse[5,11:20],sparse[1,1:10]),
EA=c(sparse[6,11:20],sparse[1,1:10]),
QA=c(sparse[7,11:20],sparse[1,1:10]),
GA=c(sparse[8,11:20],sparse[1,1:10]  ),
HA=c(sparse[9,11:20],sparse[1,1:10]  ),
IA=c(sparse[10,11:20],sparse[1,1:10]  ),
LA=c(sparse[11,11:20],sparse[1,1:10]  ),
KA=c(sparse[12,11:20],sparse[1,1:10]  ),
MA=c(sparse[13,11:20],sparse[1,1:10]  ),
FA=c(sparse[14,11:20],sparse[1,1:10]  ),
PA=c(sparse[15,11:20],sparse[1,1:10]  ),
SA=c(sparse[16,11:20],sparse[1,1:10]  ),
TA=c(sparse[17,11:20],sparse[1,1:10]  ),
WA=c(sparse[18,11:20],sparse[1,1:10]  ),
YA=c(sparse[19,11:20],sparse[1,1:10]  ),
VA=c(sparse[20,11:20],sparse[1,1:10]  ),

AR=c(sparse[1,11:20],sparse[2,1:10]  ),
RR=c(sparse[2,11:20],sparse[2,1:10]  ),
NR=c(sparse[3,11:20],sparse[2,1:10]  ),
DR=c(sparse[4,11:20],sparse[2,1:10]  ),
CR=c(sparse[5,11:20],sparse[2,1:10]  ),
ER=c(sparse[6,11:20],sparse[2,1:10]  ),
QR=c(sparse[7,11:20],sparse[2,1:10]  ),
GR=c(sparse[8,11:20],sparse[2,1:10]  ),
HR=c(sparse[9,11:20],sparse[2,1:10]  ),
IR=c(sparse[10,11:20],sparse[2,1:10]  ),
LR=c(sparse[11,11:20],sparse[2,1:10]  ),
KR=c(sparse[12,11:20],sparse[2,1:10]  ),
MR=c(sparse[13,11:20],sparse[2,1:10]  ),
FR=c(sparse[14,11:20],sparse[2,1:10]  ),
PR=c(sparse[15,11:20],sparse[2,1:10]  ),
SR=c(sparse[16,11:20],sparse[2,1:10]  ),
TR=c(sparse[17,11:20],sparse[2,1:10]  ),
WR=c(sparse[18,11:20],sparse[2,1:10]  ),
YR=c(sparse[19,11:20],sparse[2,1:10]  ),
VR=c(sparse[20,11:20],sparse[2,1:10]  ),

AN=c(sparse[1,11:20],sparse[3,1:10]  ),
RN=c(sparse[2,11:20],sparse[3,1:10]  ),
NN=c(sparse[3,11:20],sparse[3,1:10]  ),
DN=c(sparse[4,11:20],sparse[3,1:10]  ),
CN=c(sparse[5,11:20],sparse[3,1:10]  ),
EN=c(sparse[6,11:20],sparse[3,1:10]  ),
QN=c(sparse[7,11:20],sparse[3,1:10]  ),
GN=c(sparse[8,11:20],sparse[3,1:10]  ),
HN=c(sparse[9,11:20],sparse[3,1:10]  ),
IN=c(sparse[10,11:20],sparse[3,1:10]  ),
LN=c(sparse[11,11:20],sparse[3,1:10]  ),
KN=c(sparse[12,11:20],sparse[3,1:10]  ),
MN=c(sparse[13,11:20],sparse[3,1:10]  ),
FN=c(sparse[14,11:20],sparse[3,1:10]  ),
PN=c(sparse[15,11:20],sparse[3,1:10]  ),
SN=c(sparse[16,11:20],sparse[3,1:10]  ),
TN=c(sparse[17,11:20],sparse[3,1:10]  ),
WN=c(sparse[18,11:20],sparse[3,1:10]  ),
YN=c(sparse[19,11:20],sparse[3,1:10]  ),
VN=c(sparse[20,11:20],sparse[3,1:10]  ),

AD=c(sparse[1,11:20],sparse[4,1:10]  ),
RD=c(sparse[2,11:20],sparse[4,1:10]  ),
ND=c(sparse[3,11:20],sparse[4,1:10]  ),
DD=c(sparse[4,11:20],sparse[4,1:10]  ),
CD=c(sparse[5,11:20],sparse[4,1:10]  ),
ED=c(sparse[6,11:20],sparse[4,1:10]  ),
QD=c(sparse[7,11:20],sparse[4,1:10]  ),
GD=c(sparse[8,11:20],sparse[4,1:10]  ),
HD=c(sparse[9,11:20],sparse[4,1:10]  ),
ID=c(sparse[10,11:20],sparse[4,1:10]  ),
LD=c(sparse[11,11:20],sparse[4,1:10]  ),
KD=c(sparse[12,11:20],sparse[4,1:10]  ),
MD=c(sparse[13,11:20],sparse[4,1:10]  ),
FD=c(sparse[14,11:20],sparse[4,1:10]  ),
PD=c(sparse[15,11:20],sparse[4,1:10]  ),
SD=c(sparse[16,11:20],sparse[4,1:10]  ),
TD=c(sparse[17,11:20],sparse[4,1:10]  ),
WD=c(sparse[18,11:20],sparse[4,1:10]  ),
YD=c(sparse[19,11:20],sparse[4,1:10]  ),
VD=c(sparse[20,11:20],sparse[4,1:10]  ),


AC=c(sparse[1,11:20],sparse[5,1:10]  ),
RC=c(sparse[2,11:20],sparse[5,1:10]  ),
NC=c(sparse[3,11:20],sparse[5,1:10]  ),
DC=c(sparse[4,11:20],sparse[5,1:10]  ),
CC=c(sparse[5,11:20],sparse[5,1:10]  ),
EC=c(sparse[6,11:20],sparse[5,1:10]  ),
QC=c(sparse[7,11:20],sparse[5,1:10]  ),
GC=c(sparse[8,11:20],sparse[5,1:10]  ),
HC=c(sparse[9,11:20],sparse[5,1:10]  ),
IC=c(sparse[10,11:20],sparse[5,1:10]  ),
LC=c(sparse[11,11:20],sparse[5,1:10]  ),
KC=c(sparse[12,11:20],sparse[5,1:10]  ),
MC=c(sparse[13,11:20],sparse[5,1:10]  ),
FC=c(sparse[14,11:20],sparse[5,1:10]  ),
PC=c(sparse[15,11:20],sparse[5,1:10]  ),
SC=c(sparse[16,11:20],sparse[5,1:10]  ),
TC=c(sparse[17,11:20],sparse[5,1:10]  ),
WC=c(sparse[18,11:20],sparse[5,1:10]  ),
YC=c(sparse[19,11:20],sparse[5,1:10]  ),
VC=c(sparse[20,11:20],sparse[5,1:10]  ),


AE=c(sparse[1,11:20],sparse[6,1:10]  ),
RE=c(sparse[2,11:20],sparse[6,1:10]  ),
NE=c(sparse[3,11:20],sparse[6,1:10]  ),
DE=c(sparse[4,11:20],sparse[6,1:10]  ),
CE=c(sparse[5,11:20],sparse[6,1:10]  ),
EE=c(sparse[6,11:20],sparse[6,1:10]  ),
QE=c(sparse[7,11:20],sparse[6,1:10]  ),
GE=c(sparse[8,11:20],sparse[6,1:10]  ),
HE=c(sparse[9,11:20],sparse[6,1:10]  ),
IE=c(sparse[10,11:20],sparse[6,1:10]  ),
LE=c(sparse[11,11:20],sparse[6,1:10]  ),
KE=c(sparse[12,11:20],sparse[6,1:10]  ),
ME=c(sparse[13,11:20],sparse[6,1:10]  ),
FE=c(sparse[14,11:20],sparse[6,1:10]  ),
PE=c(sparse[15,11:20],sparse[6,1:10]  ),
SE=c(sparse[16,11:20],sparse[6,1:10]  ),
TE=c(sparse[17,11:20],sparse[6,1:10]  ),
WE=c(sparse[18,11:20],sparse[6,1:10]  ),
YE=c(sparse[19,11:20],sparse[6,1:10]  ),
VE=c(sparse[20,11:20],sparse[6,1:10]  ),



AQ=c(sparse[1,11:20],sparse[7,1:10]  ),
RQ=c(sparse[2,11:20],sparse[7,1:10]  ),
NQ=c(sparse[3,11:20],sparse[7,1:10]  ),
DQ=c(sparse[4,11:20],sparse[7,1:10]  ),
CQ=c(sparse[5,11:20],sparse[7,1:10]  ),
EQ=c(sparse[6,11:20],sparse[7,1:10]  ),
QQ=c(sparse[7,11:20],sparse[7,1:10]  ),
GQ=c(sparse[8,11:20],sparse[7,1:10]  ),
HQ=c(sparse[9,11:20],sparse[7,1:10]  ),
IQ=c(sparse[10,11:20],sparse[7,1:10]  ),
LQ=c(sparse[11,11:20],sparse[7,1:10]  ),
KQ=c(sparse[12,11:20],sparse[7,1:10]  ),
MQ=c(sparse[13,11:20],sparse[7,1:10]  ),
FQ=c(sparse[14,11:20],sparse[7,1:10]  ),
PQ=c(sparse[15,11:20],sparse[7,1:10]  ),
SQ=c(sparse[16,11:20],sparse[7,1:10]  ),
TQ=c(sparse[17,11:20],sparse[7,1:10]  ),
WQ=c(sparse[18,11:20],sparse[7,1:10]  ),
YQ=c(sparse[19,11:20],sparse[7,1:10]  ),
VQ=c(sparse[20,11:20],sparse[7,1:10]  ),

AG=c(sparse[1,11:20],sparse[8,1:10]  ),
RG=c(sparse[2,11:20],sparse[8,1:10]  ),
NG=c(sparse[3,11:20],sparse[8,1:10]  ),
DG=c(sparse[4,11:20],sparse[8,1:10]  ),
CG=c(sparse[5,11:20],sparse[8,1:10]  ),
EG=c(sparse[6,11:20],sparse[8,1:10]  ),
QG=c(sparse[7,11:20],sparse[8,1:10]  ),
GG=c(sparse[8,11:20],sparse[8,1:10]  ),
HG=c(sparse[9,11:20],sparse[8,1:10]  ),
IG=c(sparse[10,11:20],sparse[8,1:10]  ),
LG=c(sparse[11,11:20],sparse[8,1:10]  ),
KG=c(sparse[12,11:20],sparse[8,1:10]  ),
MG=c(sparse[13,11:20],sparse[8,1:10]  ),
FG=c(sparse[14,11:20],sparse[8,1:10]  ),
PG=c(sparse[15,11:20],sparse[8,1:10]  ),
SG=c(sparse[16,11:20],sparse[8,1:10]  ),
TG=c(sparse[17,11:20],sparse[8,1:10]  ),
WG=c(sparse[18,11:20],sparse[8,1:10]  ),
YG=c(sparse[19,11:20],sparse[8,1:10]  ),
VG=c(sparse[20,11:20],sparse[8,1:10]  ),

AH=c(sparse[1,11:20],sparse[9,1:10]  ),
RH=c(sparse[2,11:20],sparse[9,1:10]  ),
NH=c(sparse[3,11:20],sparse[9,1:10]  ),
DH=c(sparse[4,11:20],sparse[9,1:10]  ),
CH=c(sparse[5,11:20],sparse[9,1:10]  ),
EH=c(sparse[6,11:20],sparse[9,1:10]  ),
QH=c(sparse[7,11:20],sparse[9,1:10]  ),
GH=c(sparse[8,11:20],sparse[9,1:10]  ),
HH=c(sparse[9,11:20],sparse[9,1:10]  ),
IH=c(sparse[10,11:20],sparse[9,1:10]  ),
LH=c(sparse[11,11:20],sparse[9,1:10]  ),
KH=c(sparse[12,11:20],sparse[9,1:10]  ),
MH=c(sparse[13,11:20],sparse[9,1:10]  ),
FH=c(sparse[14,11:20],sparse[9,1:10]  ),
PH=c(sparse[15,11:20],sparse[9,1:10]  ),
SH=c(sparse[16,11:20],sparse[9,1:10]  ),
TH=c(sparse[17,11:20],sparse[9,1:10]  ),
WH=c(sparse[18,11:20],sparse[9,1:10]  ),
YH=c(sparse[19,11:20],sparse[9,1:10]  ),
VH=c(sparse[20,11:20],sparse[9,1:10]  ),

AI=c(sparse[1,11:20],sparse[10,1:10]  ),
RI=c(sparse[2,11:20],sparse[10,1:10]  ),
NI=c(sparse[3,11:20],sparse[10,1:10]  ),
DI=c(sparse[4,11:20],sparse[10,1:10]  ),
CI=c(sparse[5,11:20],sparse[10,1:10]  ),
EI=c(sparse[6,11:20],sparse[10,1:10]  ),
QI=c(sparse[7,11:20],sparse[10,1:10]  ),
GI=c(sparse[8,11:20],sparse[10,1:10]  ),
HI=c(sparse[9,11:20],sparse[10,1:10]  ),
II=c(sparse[10,11:20],sparse[10,1:10]  ),
LI=c(sparse[11,11:20],sparse[10,1:10]  ),
KI=c(sparse[12,11:20],sparse[10,1:10]  ),
MI=c(sparse[13,11:20],sparse[10,1:10]  ),
FI=c(sparse[14,11:20],sparse[10,1:10]  ),
PI=c(sparse[15,11:20],sparse[10,1:10]  ),
SI=c(sparse[16,11:20],sparse[10,1:10]  ),
TI=c(sparse[17,11:20],sparse[10,1:10]  ),
WI=c(sparse[18,11:20],sparse[10,1:10]  ),
YI=c(sparse[19,11:20],sparse[10,1:10]  ),
VI=c(sparse[20,11:20],sparse[10,1:10]  ),

AL=c(sparse[1,11:20],sparse[11,1:10]  ),
RL=c(sparse[2,11:20],sparse[11,1:10]  ),
NL=c(sparse[3,11:20],sparse[11,1:10]  ),
DL=c(sparse[4,11:20],sparse[11,1:10]  ),
CL=c(sparse[5,11:20],sparse[11,1:10]  ),
EL=c(sparse[6,11:20],sparse[11,1:10]  ),
QL=c(sparse[7,11:20],sparse[11,1:10]  ),
GL=c(sparse[8,11:20],sparse[11,1:10]  ),
HL=c(sparse[9,11:20],sparse[11,1:10]  ),
IL=c(sparse[10,11:20],sparse[11,1:10]  ),
LL=c(sparse[11,11:20],sparse[11,1:10]  ),
KL=c(sparse[12,11:20],sparse[11,1:10]  ),
ML=c(sparse[13,11:20],sparse[11,1:10]  ),
FL=c(sparse[14,11:20],sparse[11,1:10]  ),
PL=c(sparse[15,11:20],sparse[11,1:10]  ),
SL=c(sparse[16,11:20],sparse[11,1:10]  ),
TL=c(sparse[17,11:20],sparse[11,1:10]  ),
WL=c(sparse[18,11:20],sparse[11,1:10]  ),
YL=c(sparse[19,11:20],sparse[11,1:10]  ),
VL=c(sparse[20,11:20],sparse[11,1:10]  ),

AK=c(sparse[1,11:20],sparse[12,1:10]  ),
RK=c(sparse[2,11:20],sparse[12,1:10]  ),
NK=c(sparse[3,11:20],sparse[12,1:10]  ),
DK=c(sparse[4,11:20],sparse[12,1:10]  ),
CK=c(sparse[5,11:20],sparse[12,1:10]  ),
EK=c(sparse[6,11:20],sparse[12,1:10]  ),
QK=c(sparse[7,11:20],sparse[12,1:10]  ),
GK=c(sparse[8,11:20],sparse[12,1:10]  ),
HK=c(sparse[9,11:20],sparse[12,1:10]  ),
IK=c(sparse[10,11:20],sparse[12,1:10]  ),
LK=c(sparse[11,11:20],sparse[12,1:10]  ),
KK=c(sparse[12,11:20],sparse[12,1:10]  ),
MK=c(sparse[13,11:20],sparse[12,1:10]  ),
FK=c(sparse[14,11:20],sparse[12,1:10]  ),
PK=c(sparse[15,11:20],sparse[12,1:10]  ),
SK=c(sparse[16,11:20],sparse[12,1:10]  ),
TK=c(sparse[17,11:20],sparse[12,1:10]  ),
WK=c(sparse[18,11:20],sparse[12,1:10]  ),
YK=c(sparse[19,11:20],sparse[12,1:10]  ),
VK=c(sparse[20,11:20],sparse[12,1:10]  ),

AM=c(sparse[1,11:20],sparse[13,1:10]  ),
RM=c(sparse[2,11:20],sparse[13,1:10]  ),
NM=c(sparse[3,11:20],sparse[13,1:10]  ),
DM=c(sparse[4,11:20],sparse[13,1:10]  ),
CM=c(sparse[5,11:20],sparse[13,1:10]  ),
EM=c(sparse[6,11:20],sparse[13,1:10]  ),
QM=c(sparse[7,11:20],sparse[13,1:10]  ),
GM=c(sparse[8,11:20],sparse[13,1:10]  ),
HM=c(sparse[9,11:20],sparse[13,1:10]  ),
IM=c(sparse[10,11:20],sparse[13,1:10]  ),
LM=c(sparse[11,11:20],sparse[13,1:10]  ),
KM=c(sparse[12,11:20],sparse[13,1:10]  ),
MM=c(sparse[13,11:20],sparse[13,1:10]  ),
FM=c(sparse[14,11:20],sparse[13,1:10]  ),
PM=c(sparse[15,11:20],sparse[13,1:10]  ),
SM=c(sparse[16,11:20],sparse[13,1:10]  ),
TM=c(sparse[17,11:20],sparse[13,1:10]  ),
WM=c(sparse[18,11:20],sparse[13,1:10]  ),
YM=c(sparse[19,11:20],sparse[13,1:10]  ),
VM=c(sparse[20,11:20],sparse[13,1:10]  ),


AF=c(sparse[1,11:20],sparse[14,1:10]  ),
RF=c(sparse[2,11:20],sparse[14,1:10]  ),
NF=c(sparse[3,11:20],sparse[14,1:10]  ),
DF=c(sparse[4,11:20],sparse[14,1:10]  ),
CF=c(sparse[5,11:20],sparse[14,1:10]  ),
EF=c(sparse[6,11:20],sparse[14,1:10]  ),
QF=c(sparse[7,11:20],sparse[14,1:10]  ),
GF=c(sparse[8,11:20],sparse[14,1:10]  ),
HF=c(sparse[9,11:20],sparse[14,1:10]  ),
IF=c(sparse[10,11:20],sparse[14,1:10]  ),
LF=c(sparse[11,11:20],sparse[14,1:10]  ),
KF=c(sparse[12,11:20],sparse[14,1:10]  ),
MF=c(sparse[13,11:20],sparse[14,1:10]  ),
FF=c(sparse[14,11:20],sparse[14,1:10]  ),
PF=c(sparse[15,11:20],sparse[14,1:10]  ),
SF=c(sparse[16,11:20],sparse[14,1:10]  ),
TF=c(sparse[17,11:20],sparse[14,1:10]  ),
WF=c(sparse[18,11:20],sparse[14,1:10]  ),
YF=c(sparse[19,11:20],sparse[14,1:10]  ),
VF=c(sparse[20,11:20],sparse[14,1:10]  ),


AP=c(sparse[1,11:20],sparse[15,1:10]  ),
RP=c(sparse[2,11:20],sparse[15,1:10]  ),
NP=c(sparse[3,11:20],sparse[15,1:10]  ),
DP=c(sparse[4,11:20],sparse[15,1:10]  ),
CP=c(sparse[5,11:20],sparse[15,1:10]  ),
EP=c(sparse[6,11:20],sparse[15,1:10]  ),
QP=c(sparse[7,11:20],sparse[15,1:10]  ),
GP=c(sparse[8,11:20],sparse[15,1:10]  ),
HP=c(sparse[9,11:20],sparse[15,1:10]  ),
IP=c(sparse[10,11:20],sparse[15,1:10]  ),
LP=c(sparse[11,11:20],sparse[15,1:10]  ),
KP=c(sparse[12,11:20],sparse[15,1:10]  ),
MP=c(sparse[13,11:20],sparse[15,1:10]  ),
FP=c(sparse[14,11:20],sparse[15,1:10]  ),
PP=c(sparse[15,11:20],sparse[15,1:10]  ),
SP=c(sparse[16,11:20],sparse[15,1:10]  ),
TP=c(sparse[17,11:20],sparse[15,1:10]  ),
WP=c(sparse[18,11:20],sparse[15,1:10]  ),
YP=c(sparse[19,11:20],sparse[15,1:10]  ),
VP=c(sparse[20,11:20],sparse[15,1:10]  ),


AS=c(sparse[1,11:20],sparse[16,1:10]  ),
RS=c(sparse[2,11:20],sparse[16,1:10]  ),
NS=c(sparse[3,11:20],sparse[16,1:10]  ),
DS=c(sparse[4,11:20],sparse[16,1:10]  ),
CS=c(sparse[5,11:20],sparse[16,1:10]  ),
ES=c(sparse[6,11:20],sparse[16,1:10]  ),
QS=c(sparse[7,11:20],sparse[16,1:10]  ),
GS=c(sparse[8,11:20],sparse[16,1:10]  ),
HS=c(sparse[9,11:20],sparse[16,1:10]  ),
IS=c(sparse[10,11:20],sparse[16,1:10]  ),
LS=c(sparse[11,11:20],sparse[16,1:10]  ),
KS=c(sparse[12,11:20],sparse[16,1:10]  ),
MS=c(sparse[13,11:20],sparse[16,1:10]  ),
FS=c(sparse[14,11:20],sparse[16,1:10]  ),
PS=c(sparse[15,11:20],sparse[16,1:10]  ),
SS=c(sparse[16,11:20],sparse[16,1:10]  ),
TS=c(sparse[17,11:20],sparse[16,1:10]  ),
WS=c(sparse[18,11:20],sparse[16,1:10]  ),
YS=c(sparse[19,11:20],sparse[16,1:10]  ),
VS=c(sparse[20,11:20],sparse[16,1:10]  ),


AT=c(sparse[1,11:20],sparse[17,1:10]  ),
RT=c(sparse[2,11:20],sparse[17,1:10]  ),
NT=c(sparse[3,11:20],sparse[17,1:10]  ),
DT=c(sparse[4,11:20],sparse[17,1:10]  ),
CT=c(sparse[5,11:20],sparse[17,1:10]  ),
ET=c(sparse[6,11:20],sparse[17,1:10]  ),
QT=c(sparse[7,11:20],sparse[17,1:10]  ),
GT=c(sparse[8,11:20],sparse[17,1:10]  ),
HT=c(sparse[9,11:20],sparse[17,1:10]  ),
IT=c(sparse[10,11:20],sparse[17,1:10]  ),
LT=c(sparse[11,11:20],sparse[17,1:10]  ),
KT=c(sparse[12,11:20],sparse[17,1:10]  ),
MT=c(sparse[13,11:20],sparse[17,1:10]  ),
FT=c(sparse[14,11:20],sparse[17,1:10]  ),
PT=c(sparse[15,11:20],sparse[17,1:10]  ),
ST=c(sparse[16,11:20],sparse[17,1:10]  ),
TT=c(sparse[17,11:20],sparse[17,1:10]  ),
WT=c(sparse[18,11:20],sparse[17,1:10]  ),
YT=c(sparse[19,11:20],sparse[17,1:10]  ),
VT=c(sparse[20,11:20],sparse[17,1:10]  ),

AW=c(sparse[1,11:20],sparse[18,1:10]  ),
RW=c(sparse[2,11:20],sparse[18,1:10]  ),
NW=c(sparse[3,11:20],sparse[18,1:10]  ),
DW=c(sparse[4,11:20],sparse[18,1:10]  ),
CW=c(sparse[5,11:20],sparse[18,1:10]  ),
EW=c(sparse[6,11:20],sparse[18,1:10]  ),
QW=c(sparse[7,11:20],sparse[18,1:10]  ),
GW=c(sparse[8,11:20],sparse[18,1:10]  ),
HW=c(sparse[9,11:20],sparse[18,1:10]  ),
IW=c(sparse[10,11:20],sparse[18,1:10]  ),
LW=c(sparse[11,11:20],sparse[18,1:10]  ),
KW=c(sparse[12,11:20],sparse[18,1:10]  ),
MW=c(sparse[13,11:20],sparse[18,1:10]  ),
FW=c(sparse[14,11:20],sparse[18,1:10]  ),
PW=c(sparse[15,11:20],sparse[18,1:10]  ),
SW=c(sparse[16,11:20],sparse[18,1:10]  ),
TW=c(sparse[17,11:20],sparse[18,1:10]  ),
WW=c(sparse[18,11:20],sparse[18,1:10]  ),
YW=c(sparse[19,11:20],sparse[18,1:10]  ),
VW=c(sparse[20,11:20],sparse[18,1:10]  ),

AY=c(sparse[1,11:20],sparse[19,1:10]  ),
RY=c(sparse[2,11:20],sparse[19,1:10]  ),
NY=c(sparse[3,11:20],sparse[19,1:10]  ),
DY=c(sparse[4,11:20],sparse[19,1:10]  ),
CY=c(sparse[5,11:20],sparse[19,1:10]  ),
EY=c(sparse[6,11:20],sparse[19,1:10]  ),
QY=c(sparse[7,11:20],sparse[19,1:10]  ),
GY=c(sparse[8,11:20],sparse[19,1:10]  ),
HY=c(sparse[9,11:20],sparse[19,1:10]  ),
IY=c(sparse[10,11:20],sparse[19,1:10]  ),
LY=c(sparse[11,11:20],sparse[19,1:10]  ),
KY=c(sparse[12,11:20],sparse[19,1:10]  ),
MY=c(sparse[13,11:20],sparse[19,1:10]  ),
FY=c(sparse[14,11:20],sparse[19,1:10]  ),
PY=c(sparse[15,11:20],sparse[19,1:10]  ),
SY=c(sparse[16,11:20],sparse[19,1:10]  ),
TY=c(sparse[17,11:20],sparse[19,1:10]  ),
WY=c(sparse[18,11:20],sparse[19,1:10]  ),
YY=c(sparse[19,11:20],sparse[19,1:10]  ),
VY=c(sparse[20,11:20],sparse[19,1:10]  ),


AV=c(sparse[1,11:20],sparse[20,1:10]  ),
RV=c(sparse[2,11:20],sparse[20,1:10]  ),
NV=c(sparse[3,11:20],sparse[20,1:10]  ),
DV=c(sparse[4,11:20],sparse[20,1:10]  ),
CV=c(sparse[5,11:20],sparse[20,1:10]  ),
EV=c(sparse[6,11:20],sparse[20,1:10]  ),
QV=c(sparse[7,11:20],sparse[20,1:10]  ),
GV=c(sparse[8,11:20],sparse[20,1:10]  ),
HV=c(sparse[9,11:20],sparse[20,1:10]  ),
IV=c(sparse[10,11:20],sparse[20,1:10]  ),
LV=c(sparse[11,11:20],sparse[20,1:10]  ),
KV=c(sparse[12,11:20],sparse[20,1:10]  ),
MV=c(sparse[13,11:20],sparse[20,1:10]  ),
FV=c(sparse[14,11:20],sparse[20,1:10]  ),
PV=c(sparse[15,11:20],sparse[20,1:10]  ),
SV=c(sparse[16,11:20],sparse[20,1:10]  ),
TV=c(sparse[17,11:20],sparse[20,1:10]  ),
WV=c(sparse[18,11:20],sparse[20,1:10]  ),
YV=c(sparse[19,11:20],sparse[20,1:10]  ),
VV=c(sparse[20,11:20],sparse[20,1:10]  )))



DC1=DC[1]*disparse[1,]
DC2=DC[2]*disparse[2,]
DC3=DC[3]*disparse[3,]
DC4=DC[4]*disparse[4,]
DC5=DC[5]*disparse[5,]
DC6=DC[6]*disparse[6,]
DC7=DC[7]*disparse[7,]
DC8=DC[8]*disparse[8,]
DC9=DC[9]*disparse[9,]
DC10=DC[10]*disparse[10,]
DC11=DC[11]*disparse[11,]
DC12=DC[12]*disparse[12,]
DC13=DC[13]*disparse[13,]
DC14=DC[14]*disparse[14,]
DC15=DC[15]*disparse[15,]
DC16=DC[16]*disparse[16,]
DC17=DC[17]*disparse[17,]
DC18=DC[18]*disparse[18,]
DC19=DC[19]*disparse[19,]
DC20=DC[20]*disparse[20,]
DC21=DC[21]*disparse[21,]
DC22=DC[22]*disparse[22,]
DC23=DC[23]*disparse[23,]
DC24=DC[24]*disparse[24,]
DC25=DC[25]*disparse[25,]
DC26=DC[26]*disparse[26,]
DC27=DC[27]*disparse[27,]
DC28=DC[28]*disparse[28,]
DC29=DC[29]*disparse[29,]
DC30=DC[30]*disparse[30,]
DC31=DC[31]*disparse[31,]
DC32=DC[32]*disparse[32,]
DC33=DC[33]*disparse[33,]
DC34=DC[34]*disparse[34,]
DC35=DC[35]*disparse[35,]
DC36=DC[36]*disparse[36,]
DC37=DC[37]*disparse[37,]
DC38=DC[38]*disparse[38,]
DC39=DC[39]*disparse[39,]
DC40=DC[40]*disparse[40,]
DC41=DC[41]*disparse[41,]
DC42=DC[42]*disparse[42,]
DC43=DC[43]*disparse[43,]
DC44=DC[44]*disparse[44,]
DC45=DC[45]*disparse[45,]
DC46=DC[46]*disparse[46,]
DC47=DC[47]*disparse[47,]
DC48=DC[48]*disparse[48,]
DC49=DC[49]*disparse[49,]
DC50=DC[50]*disparse[50,]
DC51=DC[51]*disparse[51,]
DC52=DC[52]*disparse[52,]
DC53=DC[53]*disparse[53,]
DC54=DC[54]*disparse[54,]
DC55=DC[55]*disparse[55,]
DC56=DC[56]*disparse[56,]
DC57=DC[57]*disparse[57,]
DC58=DC[58]*disparse[58,]
DC59=DC[59]*disparse[59,]
DC60=DC[60]*disparse[60,]
DC61=DC[61]*disparse[61,]
DC62=DC[62]*disparse[62,]
DC63=DC[63]*disparse[63,]
DC64=DC[64]*disparse[64,]
DC65=DC[65]*disparse[65,]
DC66=DC[66]*disparse[66,]
DC67=DC[67]*disparse[67,]
DC68=DC[68]*disparse[68,]
DC69=DC[69]*disparse[69,]
DC70=DC[70]*disparse[70,]
DC71=DC[71]*disparse[71,]
DC72=DC[72]*disparse[72,]
DC73=DC[73]*disparse[73,]
DC74=DC[74]*disparse[74,]
DC75=DC[75]*disparse[75,]
DC76=DC[76]*disparse[76,]
DC77=DC[77]*disparse[77,]
DC78=DC[78]*disparse[78,]
DC79=DC[79]*disparse[79,]
DC80=DC[80]*disparse[80,]
DC81=DC[81]*disparse[81,]
DC82=DC[82]*disparse[82,]
DC83=DC[83]*disparse[83,]
DC84=DC[84]*disparse[84,]
DC85=DC[85]*disparse[85,]
DC86=DC[86]*disparse[86,]
DC87=DC[87]*disparse[87,]
DC88=DC[88]*disparse[88,]
DC89=DC[89]*disparse[89,]
DC90=DC[90]*disparse[90,]
DC91=DC[91]*disparse[91,]
DC92=DC[92]*disparse[92,]
DC93=DC[93]*disparse[93,]
DC94=DC[94]*disparse[94,]
DC95=DC[95]*disparse[95,]
DC96=DC[96]*disparse[96,]
DC97=DC[97]*disparse[97,]
DC98=DC[98]*disparse[98,]
DC99=DC[99]*disparse[99,]
DC100=DC[100]*disparse[100,]
DC101=DC[101]*disparse[101,]
DC102=DC[102]*disparse[102,]
DC103=DC[103]*disparse[103,]
DC104=DC[104]*disparse[104,]
DC105=DC[105]*disparse[105,]
DC106=DC[106]*disparse[106,]
DC107=DC[107]*disparse[107,]
DC108=DC[108]*disparse[108,]
DC109=DC[109]*disparse[109,]
DC110=DC[110]*disparse[110,]
DC111=DC[111]*disparse[111,]
DC112=DC[112]*disparse[112,]
DC113=DC[113]*disparse[113,]
DC114=DC[114]*disparse[114,]
DC115=DC[115]*disparse[115,]
DC116=DC[116]*disparse[116,]
DC117=DC[117]*disparse[117,]
DC118=DC[118]*disparse[118,]
DC119=DC[119]*disparse[119,]
DC120=DC[120]*disparse[120,]
DC121=DC[121]*disparse[121,]
DC122=DC[122]*disparse[122,]
DC123=DC[123]*disparse[123,]
DC124=DC[124]*disparse[124,]
DC125=DC[125]*disparse[125,]
DC126=DC[126]*disparse[126,]
DC127=DC[127]*disparse[127,]
DC128=DC[128]*disparse[128,]
DC129=DC[129]*disparse[129,]
DC130=DC[130]*disparse[130,]
DC131=DC[131]*disparse[131,]
DC132=DC[132]*disparse[132,]
DC133=DC[133]*disparse[133,]
DC134=DC[134]*disparse[134,]
DC135=DC[135]*disparse[135,]
DC136=DC[136]*disparse[136,]
DC137=DC[137]*disparse[137,]
DC138=DC[138]*disparse[138,]
DC139=DC[139]*disparse[139,]
DC140=DC[140]*disparse[140,]
DC141=DC[141]*disparse[141,]
DC142=DC[142]*disparse[142,]
DC143=DC[143]*disparse[143,]
DC144=DC[144]*disparse[144,]
DC145=DC[145]*disparse[145,]
DC146=DC[146]*disparse[146,]
DC147=DC[147]*disparse[147,]
DC148=DC[148]*disparse[148,]
DC149=DC[149]*disparse[149,]
DC150=DC[150]*disparse[150,]
DC151=DC[151]*disparse[151,]
DC152=DC[152]*disparse[152,]
DC153=DC[153]*disparse[153,]
DC154=DC[154]*disparse[154,]
DC155=DC[155]*disparse[155,]
DC156=DC[156]*disparse[156,]
DC157=DC[157]*disparse[157,]
DC158=DC[158]*disparse[158,]
DC159=DC[159]*disparse[159,]
DC160=DC[160]*disparse[160,]
DC161=DC[161]*disparse[161,]
DC162=DC[162]*disparse[162,]
DC163=DC[163]*disparse[163,]
DC164=DC[164]*disparse[164,]
DC165=DC[165]*disparse[165,]
DC166=DC[166]*disparse[166,]
DC167=DC[167]*disparse[167,]
DC168=DC[168]*disparse[168,]
DC169=DC[169]*disparse[169,]
DC170=DC[170]*disparse[170,]
DC171=DC[171]*disparse[171,]
DC172=DC[172]*disparse[172,]
DC173=DC[173]*disparse[173,]
DC174=DC[174]*disparse[174,]
DC175=DC[175]*disparse[175,]
DC176=DC[176]*disparse[176,]
DC177=DC[177]*disparse[177,]
DC178=DC[178]*disparse[178,]
DC179=DC[179]*disparse[179,]
DC180=DC[180]*disparse[180,]
DC181=DC[181]*disparse[181,]
DC182=DC[182]*disparse[182,]
DC183=DC[183]*disparse[183,]
DC184=DC[184]*disparse[184,]
DC185=DC[185]*disparse[185,]
DC186=DC[186]*disparse[186,]
DC187=DC[187]*disparse[187,]
DC188=DC[188]*disparse[188,]
DC189=DC[189]*disparse[189,]
DC190=DC[190]*disparse[190,]
DC191=DC[191]*disparse[191,]
DC192=DC[192]*disparse[192,]
DC193=DC[193]*disparse[193,]
DC194=DC[194]*disparse[194,]
DC195=DC[195]*disparse[195,]
DC196=DC[196]*disparse[196,]
DC197=DC[197]*disparse[197,]
DC198=DC[198]*disparse[198,]
DC199=DC[199]*disparse[199,]
DC200=DC[200]*disparse[200,]
DC201=DC[201]*disparse[201,]
DC202=DC[202]*disparse[202,]
DC203=DC[203]*disparse[203,]
DC204=DC[204]*disparse[204,]
DC205=DC[205]*disparse[205,]
DC206=DC[206]*disparse[206,]
DC207=DC[207]*disparse[207,]
DC208=DC[208]*disparse[208,]
DC209=DC[209]*disparse[209,]
DC210=DC[210]*disparse[210,]
DC211=DC[211]*disparse[211,]
DC212=DC[212]*disparse[212,]
DC213=DC[213]*disparse[213,]
DC214=DC[214]*disparse[214,]
DC215=DC[215]*disparse[215,]
DC216=DC[216]*disparse[216,]
DC217=DC[217]*disparse[217,]
DC218=DC[218]*disparse[218,]
DC219=DC[219]*disparse[219,]
DC220=DC[220]*disparse[220,]
DC221=DC[221]*disparse[221,]
DC222=DC[222]*disparse[222,]
DC223=DC[223]*disparse[223,]
DC224=DC[224]*disparse[224,]
DC225=DC[225]*disparse[225,]
DC226=DC[226]*disparse[226,]
DC227=DC[227]*disparse[227,]
DC228=DC[228]*disparse[228,]
DC229=DC[229]*disparse[229,]
DC230=DC[230]*disparse[230,]
DC231=DC[231]*disparse[231,]
DC232=DC[232]*disparse[232,]
DC233=DC[233]*disparse[233,]
DC234=DC[234]*disparse[234,]
DC235=DC[235]*disparse[235,]
DC236=DC[236]*disparse[236,]
DC237=DC[237]*disparse[237,]
DC238=DC[238]*disparse[238,]
DC239=DC[239]*disparse[239,]
DC240=DC[240]*disparse[240,]
DC241=DC[241]*disparse[241,]
DC242=DC[242]*disparse[242,]
DC243=DC[243]*disparse[243,]
DC244=DC[244]*disparse[244,]
DC245=DC[245]*disparse[245,]
DC246=DC[246]*disparse[246,]
DC247=DC[247]*disparse[247,]
DC248=DC[248]*disparse[248,]
DC249=DC[249]*disparse[249,]
DC250=DC[250]*disparse[250,]
DC251=DC[251]*disparse[251,]
DC252=DC[252]*disparse[252,]
DC253=DC[253]*disparse[253,]
DC254=DC[254]*disparse[254,]
DC255=DC[255]*disparse[255,]
DC256=DC[256]*disparse[256,]
DC257=DC[257]*disparse[257,]
DC258=DC[258]*disparse[258,]
DC259=DC[259]*disparse[259,]
DC260=DC[260]*disparse[260,]
DC261=DC[261]*disparse[261,]
DC262=DC[262]*disparse[262,]
DC263=DC[263]*disparse[263,]
DC264=DC[264]*disparse[264,]
DC265=DC[265]*disparse[265,]
DC266=DC[266]*disparse[266,]
DC267=DC[267]*disparse[267,]
DC268=DC[268]*disparse[268,]
DC269=DC[269]*disparse[269,]
DC270=DC[270]*disparse[270,]
DC271=DC[271]*disparse[271,]
DC272=DC[272]*disparse[272,]
DC273=DC[273]*disparse[273,]
DC274=DC[274]*disparse[274,]
DC275=DC[275]*disparse[275,]
DC276=DC[276]*disparse[276,]
DC277=DC[277]*disparse[277,]
DC278=DC[278]*disparse[278,]
DC279=DC[279]*disparse[279,]
DC280=DC[280]*disparse[280,]
DC281=DC[281]*disparse[281,]
DC282=DC[282]*disparse[282,]
DC283=DC[283]*disparse[283,]
DC284=DC[284]*disparse[284,]
DC285=DC[285]*disparse[285,]
DC286=DC[286]*disparse[286,]
DC287=DC[287]*disparse[287,]
DC288=DC[288]*disparse[288,]
DC289=DC[289]*disparse[289,]
DC290=DC[290]*disparse[290,]
DC291=DC[291]*disparse[291,]
DC292=DC[292]*disparse[292,]
DC293=DC[293]*disparse[293,]
DC294=DC[294]*disparse[294,]
DC295=DC[295]*disparse[295,]
DC296=DC[296]*disparse[296,]
DC297=DC[297]*disparse[297,]
DC298=DC[298]*disparse[298,]
DC299=DC[299]*disparse[299,]
DC300=DC[300]*disparse[300,]
DC301=DC[301]*disparse[301,]
DC302=DC[302]*disparse[302,]
DC303=DC[303]*disparse[303,]
DC304=DC[304]*disparse[304,]
DC305=DC[305]*disparse[305,]
DC306=DC[306]*disparse[306,]
DC307=DC[307]*disparse[307,]
DC308=DC[308]*disparse[308,]
DC309=DC[309]*disparse[309,]
DC310=DC[310]*disparse[310,]
DC311=DC[311]*disparse[311,]
DC312=DC[312]*disparse[312,]
DC313=DC[313]*disparse[313,]
DC314=DC[314]*disparse[314,]
DC315=DC[315]*disparse[315,]
DC316=DC[316]*disparse[316,]
DC317=DC[317]*disparse[317,]
DC318=DC[318]*disparse[318,]
DC319=DC[319]*disparse[319,]
DC320=DC[320]*disparse[320,]
DC321=DC[321]*disparse[321,]
DC322=DC[322]*disparse[322,]
DC323=DC[323]*disparse[323,]
DC324=DC[324]*disparse[324,]
DC325=DC[325]*disparse[325,]
DC326=DC[326]*disparse[326,]
DC327=DC[327]*disparse[327,]
DC328=DC[328]*disparse[328,]
DC329=DC[329]*disparse[329,]
DC330=DC[330]*disparse[330,]
DC331=DC[331]*disparse[331,]
DC332=DC[332]*disparse[332,]
DC333=DC[333]*disparse[333,]
DC334=DC[334]*disparse[334,]
DC335=DC[335]*disparse[335,]
DC336=DC[336]*disparse[336,]
DC337=DC[337]*disparse[337,]
DC338=DC[338]*disparse[338,]
DC339=DC[339]*disparse[339,]
DC340=DC[340]*disparse[340,]
DC341=DC[341]*disparse[341,]
DC342=DC[342]*disparse[342,]
DC343=DC[343]*disparse[343,]
DC344=DC[344]*disparse[344,]
DC345=DC[345]*disparse[345,]
DC346=DC[346]*disparse[346,]
DC347=DC[347]*disparse[347,]
DC348=DC[348]*disparse[348,]
DC349=DC[349]*disparse[349,]
DC350=DC[350]*disparse[350,]
DC351=DC[351]*disparse[351,]
DC352=DC[352]*disparse[352,]
DC353=DC[353]*disparse[353,]
DC354=DC[354]*disparse[354,]
DC355=DC[355]*disparse[355,]
DC356=DC[356]*disparse[356,]
DC357=DC[357]*disparse[357,]
DC358=DC[358]*disparse[358,]
DC359=DC[359]*disparse[359,]
DC360=DC[360]*disparse[360,]
DC361=DC[361]*disparse[361,]
DC362=DC[362]*disparse[362,]
DC363=DC[363]*disparse[363,]
DC364=DC[364]*disparse[364,]
DC365=DC[365]*disparse[365,]
DC366=DC[366]*disparse[366,]
DC367=DC[367]*disparse[367,]
DC368=DC[368]*disparse[368,]
DC369=DC[369]*disparse[369,]
DC370=DC[370]*disparse[370,]
DC371=DC[371]*disparse[371,]
DC372=DC[372]*disparse[372,]
DC373=DC[373]*disparse[373,]
DC374=DC[374]*disparse[374,]
DC375=DC[375]*disparse[375,]
DC376=DC[376]*disparse[376,]
DC377=DC[377]*disparse[377,]
DC378=DC[378]*disparse[378,]
DC379=DC[379]*disparse[379,]
DC380=DC[380]*disparse[380,]
DC381=DC[381]*disparse[381,]
DC382=DC[382]*disparse[382,]
DC383=DC[383]*disparse[383,]
DC384=DC[384]*disparse[384,]
DC385=DC[385]*disparse[385,]
DC386=DC[386]*disparse[386,]
DC387=DC[387]*disparse[387,]
DC388=DC[388]*disparse[388,]
DC389=DC[389]*disparse[389,]
DC390=DC[390]*disparse[390,]
DC391=DC[391]*disparse[391,]
DC392=DC[392]*disparse[392,]
DC393=DC[393]*disparse[393,]
DC394=DC[394]*disparse[394,]
DC395=DC[395]*disparse[395,]
DC396=DC[396]*disparse[396,]
DC397=DC[397]*disparse[397,]
DC398=DC[398]*disparse[398,]
DC399=DC[399]*disparse[399,]
DC400=DC[400]*disparse[400,]

result2=(DC1+DC2+DC3+DC4+DC5+DC6+DC7+DC8+DC9+DC10+DC11+DC12+DC13+DC14+DC15+DC16+DC17+DC18+DC19+DC20+DC21+DC22+DC23+DC24+DC25+DC26+DC27+DC28+DC29+DC30+DC31+DC32+DC33+DC34+DC35+DC36+DC37+DC38+DC39+DC40+DC41+DC42+DC43+DC44+DC45+DC46+DC47+DC48+DC49+DC50+DC51+DC52+DC53+DC54+DC55+DC56+DC57+DC58+DC59+DC60+DC61+DC62+DC63+DC64+DC65+DC66+DC67+DC68+DC69+DC70+DC71+DC72+DC73+DC74+DC75+DC76+DC77+DC78+DC79+DC80+DC81+DC82+DC83+DC84+DC85+DC86+DC87+DC88+DC89+DC90+DC91+DC92+DC93+DC94+DC95+DC96+DC97+DC98+DC99+DC100+
DC101+DC102+DC103+DC104+DC105+DC106+DC107+DC108+DC109+DC110+DC111+DC112+DC113+DC114+DC115+DC116+DC117+DC118+DC119+DC120+DC121+DC122+DC123+DC124+DC125+DC126+DC127+DC128+DC129+DC130+DC131+DC132+DC133+DC134+DC135+DC136+DC137+DC138+DC139+DC140+DC141+DC142+DC143+DC144+DC145+DC146+DC147+DC148+DC149+DC150+DC151+DC152+DC153+DC154+DC155+DC156+DC157+DC158+DC159+DC160+DC161+DC162+DC163+DC164+DC165+DC166+DC167+DC168+DC169+DC170+DC171+DC172+DC173+DC174+DC175+DC176+DC177+DC178+DC179+DC180+DC181+DC182+DC183+DC184+DC185+DC186+DC187+DC188+DC189+DC190+DC191+DC192+DC193+DC194+DC195+DC196+DC197+DC198+DC199+DC200+
DC201+DC202+DC203+DC204+DC205+DC206+DC207+DC208+DC209+DC210+DC211+DC212+DC213+DC214+DC215+DC216+DC217+DC218+DC219+DC220+DC221+DC222+DC223+DC224+DC225+DC226+DC227+DC228+DC229+DC230+DC231+DC232+DC233+DC234+DC235+DC236+DC237+DC238+DC239+DC240+DC241+DC242+DC243+DC244+DC245+DC246+DC247+DC248+DC249+DC250+DC251+DC252+DC253+DC254+DC255+DC256+DC257+DC258+DC259+DC260+DC261+DC262+DC263+DC264+DC265+DC266+DC267+DC268+DC269+DC270+DC271+DC272+DC273+DC274+DC275+DC276+DC277+DC278+DC279+DC280+DC281+DC282+DC283+DC284+DC285+DC286+DC287+DC288+DC289+DC290+DC291+DC292+DC293+DC294+DC295+DC296+DC297+DC298+DC299+DC300+
DC301+DC302+DC303+DC304+DC305+DC306+DC307+DC308+DC309+DC310+DC311+DC312+DC313+DC314+DC315+DC316+DC317+DC318+DC319+DC320+DC321+DC322+DC323+DC324+DC325+DC326+DC327+DC328+DC329+DC330+DC331+DC332+DC333+DC334+DC335+DC336+DC337+DC338+DC339+DC340+DC341+DC342+DC343+DC344+DC345+DC346+DC347+DC348+DC349+DC350+DC351+DC352+DC353+DC354+DC355+DC356+DC357+DC358+DC359+DC360+DC361+DC362+DC363+DC364+DC365+DC366+DC367+DC368+DC369+DC370+DC371+DC372+DC373+DC374+DC375+DC376+DC377+DC378+DC379+DC380+DC381+DC382+DC383+DC384+DC385+DC386+DC387+DC388+DC389+DC390+DC391+DC392+DC393+DC394+DC395+DC396+DC397+DC398+DC399+DC400)
result<-c(result1,result2)/n
return(result)
}


