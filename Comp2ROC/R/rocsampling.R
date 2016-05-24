rocsampling <-
function(curve1.fpr,curve1.tpr,curve2.fpr,curve2.tpr, K=100) {
  # n. de linhas de amostragem
  #K=100
  # ponto de refer?ncia
  ref.point=c(1,0)

  # calcula declives das segmentos de recta das curvas ROC     OK!!!
  curve1.segslope=curvesegslope(curve1.fpr,curve1.tpr)
  curve2.segslope=curvesegslope(curve2.fpr,curve2.tpr)
  
  # calcula declives das segmentos de recta que unem os pontos das curvas ROC ao ponto de referencia OK!!!
  curve1.slope=curvesegsloperef(curve1.fpr,curve1.tpr,ref.point)
  curve2.slope=curvesegsloperef(curve2.fpr,curve2.tpr,ref.point)

  # calcula os declives das rectas de amostragem  OK!!!
  line.slope=lineslope(K)

  # calcula os pontos de intersec??o das linhas de amostragem com os segmento de recta das curvas ROC  OK!!!
  # as distancias ao ponto de refer?ncia  OK!!!
  line.distance1=linedistance(curve1.fpr,curve1.tpr,curve1.segslope,curve1.slope,line.slope,ref.point)
  line.distance2=linedistance(curve2.fpr,curve2.tpr,curve2.segslope,curve2.slope,line.slope,ref.point)
  line.dist1=line.distance1$dist
  line.dist2=line.distance2$dist
  line.xint1=line.distance1$x
  line.xint2=line.distance2$x
  line.yint1=line.distance1$y
  line.yint2=line.distance2$y

  # calcula as ?reas dos K+1 triangulos e a area total com base nos triangulos  OK!!!
  area1=areatriangles(line.slope,line.dist1)
  area2=areatriangles(line.slope,line.dist2)
  area.tri1=area1$areatri
  area.tri2=area2$areatri

  # calcula as diferen?as das ?reas  OK!!!
  diff=diffareatriangles(area.tri1,area.tri2)
  diff.areas=diff$diffareas

# calcula as proporcoes do espa?o (com base nas areas)  OK!!!
prop.curve1=length(diff.areas[diff.areas>0])/length(diff.areas)
prop.curve2=length(diff.areas[diff.areas<0])/length(diff.areas)
prop.ties=length(diff.areas[diff.areas==0])/length(diff.areas)

# calcula as localizacoes no espaco   NOT OK???
# graficos NOT OK???
limits.curve1=c(1,line.xint1)
loc.curve1=rbind(limits.curve1[which(diff.areas>0)+1],limits.curve1[which(diff.areas>0)])
limits.curve2=c(1,line.xint2)
loc.curve2=rbind(limits.curve2[which(diff.areas<0)+1],limits.curve2[which(diff.areas<0)])
loc.ties1=rbind(limits.curve1[which(diff.areas==0)+1],limits.curve1[which(diff.areas==0)])
loc.ties2=rbind(limits.curve2[which(diff.areas==0)+1],limits.curve2[which(diff.areas==0)])

ang=sin(-atan(line.slope[1])+pi/2)

answer=list(AUC1=sum(area.tri1),AUC2=sum(area.tri2),propc1=prop.curve1,propc2=prop.curve2,propties=prop.ties,
  locc1=loc.curve1,locc2=loc.curve2,locties=loc.ties1,K=K,lineslope=seq(ang,ang*101,by=ang),diffareas=diff.areas,dist1=line.dist1,dist2=line.dist2)
return(answer)
}
