# ###############   Trim EEG Head & Creat Cap Mesh   ###############
# 
# ######**###### Trim Head Vertices to Cap ######**######
# 
# ### load head and electrode coordinate (for reference)
# data(eegcoord)
# data(eeghead)
# 
# ### convert vertices to Euclidean
# library(rgl)
# newdense=t(apply(eeghead$vb,2,asEuclidean))
# 
# 
# ### find eye points
# widx=c(which(rownames(eegcoord)=="F5"),which(rownames(eegcoord)=="F6"),which(rownames(eegcoord)=="FPZ"))
# eyeidx=which((newdense[,1]-.1)>eegcoord[widx[1],1] & (newdense[,1]+.1)<eegcoord[widx[2],1] & (newdense[,3]+1.5)<eegcoord[widx[3],3] & newdense[,2]>0)
# 
# 
# ### find lower face and neck
# tidx=c(which(rownames(eegcoord)=="NZ"),which(rownames(eegcoord)=="F9"),
#        which(rownames(eegcoord)=="FT9"),which(rownames(eegcoord)=="T9"),
#        which(rownames(eegcoord)=="TP9"),which(rownames(eegcoord)=="P9"),
#        which(rownames(eegcoord)=="PO9"),which(rownames(eegcoord)=="I1"),
#        which(rownames(eegcoord)=="IZ"),which(rownames(eegcoord)=="I2"),
#        which(rownames(eegcoord)=="F10"),which(rownames(eegcoord)=="FT10"),
#        which(rownames(eegcoord)=="T10"),which(rownames(eegcoord)=="TP10"),
#        which(rownames(eegcoord)=="P10"),which(rownames(eegcoord)=="PO10"))
# zmin=min(eegcoord[tidx,3])
# zidx=which((newdense[,3]+1)<zmin)
# 
# 
# ### find first part of ear
# widx=which(rownames(eegcoord)=="T9")
# widx=c(widx,which(rownames(eegcoord)=="TP9"))
# ear1idx=which((newdense[,2]-.1)>eegcoord[widx[2],2] & (newdense[,2]+.1)<eegcoord[widx[1],2] & newdense[,3]<(-2.95))
# 
# 
# ### find second part of ear
# ee=rbind(c(-0.7833908, 0.08973989, -0.372643),
#       c(-0.7905217, 0.08175375, -0.2937119),
#       c(-0.7965502, 0.0651892, -0.2231785),
#       c(-0.8011097, 0.03303456, -0.1711342),
#       c(-0.8033209, -0.01559477, -0.1354716),
#       c(-0.8057978, -0.07360887, -0.1147273),
#       c(-0.809656, -0.1303143, -0.1056324),
#       c(-0.8135467, -0.181978, -0.106038),
#       c(-0.8151776, -0.225369, -0.1192436),
#       c(-0.8134203, -0.2630854, -0.1444437),
#       c(-0.8101834, -0.291853, -0.1810548),
#       c(-0.8053645, -0.3132131, -0.2226144),
#       c(-0.7992824, -0.3247121, -0.266523),
#       c(-0.7911068, -0.3304008, -0.3090448),
#       c(-0.7815523, -0.3292334, -0.3521771))*10
# mex=max(abs(ee[,1]))
# ear2idx=which(abs(newdense[,1])>mex)
# 
# 
# ### maunually remove remaining points in ears
# 
# # lidx=NULL
# # f=select3d()
# # idx=f(newdense[,1],newdense[,2],newdense[,3])
# # #*# pick points here
# # tidx=which(idx)
# # tidx=tidx[is.na(match(tidx,c(eyeidx,zidx,ear1idx,ear2idx,lidx)))]
# 
# ############################### DO NOT DELETE THIS ###############################
# lidx=c(418,  782,  783, 1672, 2457, 2458, 2500, 2501, 2502, 3742, 3373, 3631,
#        3762, 116,  686, 1635, 233, 1391, 1478, 1939, 1940, 1941, 1442, 1443,
#        4258, 4300, 4301, 5255, 5557, 5559, 1267, 1900, 4960, 4613, 2795,
#        1027, 1064, 4265, 4266, 4342, 4545, 421,  529, 2463, 2549, 2557, 2748)
# ############################### DO NOT DELETE THIS ###############################
# 
# 
# 
# ### plot final results
# outidx=unique(c(eyeidx,zidx,ear1idx,ear2idx,lidx))
# open3d(); points3d(newdense[-outidx,1],newdense[-outidx,2],newdense[-outidx,3])
# text3d(eegcoord[,1],eegcoord[,2],eegcoord[,3],texts=rownames(eegcoord),col="blue")
# 
# ### make 2d projection
# rrr=max(eegcoord[,3])-min(eegcoord[,3])
# ttt=(newdense[-outidx,3]+rrr)/rrr
# eeg2d=newdense[-outidx,1:2]/ttt
# x11()
# plot(eeg2d,xlim=c(-16,16),ylim=c(-16,16))
# text(eegcoord[,4],eegcoord[,5],labels=rownames(eegcoord),col="blue",cex=1.25)
# rad=12.5
# xx=rad*cos(seq(0,2*pi,length.out=360))
# yy=rad*sin(seq(0,2*pi,length.out=360))
# lines(xx,yy)
# 
# ### combine into eegdense data frame
# eegdense=data.frame(newdense[-outidx,],eeg2d)
# colnames(eegdense)=c("x","y","z","xproj","yproj")
# 
# 
# 
# 
# ######**###### Create Cap Mesh from Trimmed Vertices ######**######
# 
# ### find vextex indices
# vfun=function(x,outidx){
#   ival=rep(FALSE,3)
#   for(j in 1:3){ival[j]=any(x[j]==outidx)}
#   ifelse(any(ival),FALSE,TRUE)
# }
# 
# vidx=apply(eeghead$it,2,vfun,outidx=outidx)
# 
# ### remap vertex indices
# inidx=(1:5637)[-outidx]
# eegit=eeghead$it[,vidx]
# newit=matrix(NA,nrow(eegit),ncol(eegit))
# for(j in 1:977){
#   idx=inidx[j]
#   widx=which(eegit==idx)
#   newit[widx]=j
# }
# 
# ### define mesh eeg cap
# eegmesh=list(vb=apply(newdense[-outidx,],1,asHomogeneous),
#              it=newit,normals=eeghead$normals[,-outidx],
#              material=list(color=rep(c("red","orange","yellow"),1829)))
# class(eegmesh)="mesh3d"
# 
# ### plot results
# open3d()
# wire3d(eegmesh)
# shade3d(eegmesh)
# text3d(eegcoord[,1],eegcoord[,2],eegcoord[,3],texts=rownames(eegcoord),col="blue")
# 
# 
# ### manually add missing face vertex indices
# 
# # eegmesh$material$color=rep("black",length(eegmesh$material$color))
# # open3d()
# # wire3d(eegmesh)
# # points3d(eegdense[,1],eegdense[,2],eegdense[,3],col="red",size=10)
# # 
# # tidx=NULL
# # f=select3d()
# # idx=f(eegdense[,1],eegdense[,2],eegdense[,3])
# # #*# pick points here
# # tidx=c(tidx,which(idx))
# 
# newface=cbind(c(410, 442, 564), c(774, 802, 934), c(271, 272, 357),
#               c(348, 349, 608), c(42, 43, 44), c(42, 269, 347),
#               c(130, 609, 610), c(129, 130, 267), c(21, 27, 295),
#               c(657, 658, 659), c(253, 672, 977), c(633, 976, 977), 
#               c(632, 634, 715), c(632, 714, 716))
# 
# eegmesh=tmesh3d(eegmesh$vb,cbind(eegmesh$it,newface),
#                 material=list(color=rep(c("black"),1843*3)))
# 
# ### slightly rescale mesh (to be visible over head)
# eegmesh=scale3d(eegmesh,1.01,1.01,1.01)
# 
# open3d()
# shade3d(eegmesh)
# 
# wire3d(eegmesh)
# shade3d(eeghead)
# 
# 
# ### save results
# save(eegdense,file="/Users/Nate/Documents/R/My_R_Code/eegkit/CRAN/CRAN_ver1.0-0/eegkit/data/eegdense.rda")
# save(eegmesh,file="/Users/Nate/Documents/R/My_R_Code/eegkit/CRAN/CRAN_ver1.0-0/eegkit/data/eegmesh.rda")
# 
# 
