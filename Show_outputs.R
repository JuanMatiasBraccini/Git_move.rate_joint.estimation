### Script for combining inputs and outputs of movement rate estimations of gummy and whiskery sharks  #####

library(expm)
require(plotrix)

#DATA SECTION
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

hndl=handl_OneDrive("Analyses/Movement rate estimation/Joint.estim_ind.base.mod")
setwd(paste(hndl,"Show Gummy and whiskery outputs",sep="/"))

#conventional tagging input data
Conv.Gummy=read.csv("Conv.Gummy.csv")
Conv.Whiskery=read.csv("Conv.Whiskery.csv")

Conv.Gummy$Rec.zn=with(Conv.Gummy,ifelse(Rec.zn==1,"West",ifelse(Rec.zn==2,"Zone1",
                    ifelse(Rec.zn==3,"Zone2",NA))))
Conv.Gummy$Rel.zn=with(Conv.Gummy,ifelse(Rel.zn==1,"West",ifelse(Rel.zn==2,"Zone1",
                     ifelse(Rel.zn==3,"Zone2",NA))))

Conv.Whiskery$Rec.zn=with(Conv.Whiskery,ifelse(Rec.zn==1,"West",ifelse(Rec.zn==2,"Zone1",
                     ifelse(Rec.zn==3,"Zone2",NA))))
Conv.Whiskery$Rel.zn=with(Conv.Whiskery,ifelse(Rel.zn==1,"West",ifelse(Rel.zn==2,"Zone1",
                     ifelse(Rel.zn==3,"Zone2",NA))))


#Acoustic tagging input data
Acous.Gummy=read.csv("Acous.Gummy.csv")
Acous.Whiskery=read.csv("Acous.Whiskery.csv")

Acous.Gummy$Rec.zn=with(Acous.Gummy,ifelse(Rec.zn==1,"West",ifelse(Rec.zn==2,"Zone1",
                   ifelse(Rec.zn==3,"Zone2",NA))))
Acous.Gummy$Rel.zn=with(Acous.Gummy,ifelse(Rel.zn==1,"West",ifelse(Rel.zn==2,"Zone1",
                   ifelse(Rel.zn==3,"Zone2",NA))))

Acous.Whiskery$Rec.zn=with(Acous.Whiskery,ifelse(Rec.zn==1,"West",ifelse(Rec.zn==2,"Zone1",
                   ifelse(Rec.zn==3,"Zone2",NA))))
Acous.Whiskery$Rel.zn=with(Acous.Whiskery,ifelse(Rel.zn==1,"West",ifelse(Rel.zn==2,"Zone1",
                   ifelse(Rel.zn==3,"Zone2",NA))))


#Raw data for mapping
Raw.Conv.Gummy=read.csv("GM_Raw.Conv.Tag.csv")
Raw.Conv.Whiskery=read.csv("WH_Raw.Conv.Tag.csv")
Raw.Acous.Gummy=read.csv("Gummy_Raw.Acous.Tag.csv")
Raw.Acous.Whiskery=read.csv("Whiskery_Raw.Acous.Tag.csv")


#Movement rate modelling outputs
Prob.mov.Gummy=read.table(paste(hndl,"Gummy/model.std",sep="/"),header = T)
Prob.mov.Whiskery=read.table(paste(hndl,"Whiskery/model.std",sep="/"),header = T)
 


#PROCEDURE SECTION

#1.Display data inputs

#1.1 plot mat
SS=c("Latitude.prev","Longitude.prev")
names(Raw.Acous.Gummy)[match(SS,names(Raw.Acous.Gummy))]=c("Lat.rels" ,"Long.rels")
names(Raw.Acous.Whiskery)[match(SS,names(Raw.Acous.Whiskery))]=c("Lat.rels", "Long.rels")

library(rgdal)
SDGDLL_zone1=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone1.shp"), layer="SDGDLL_zone1") 
SDGDLL_zone2=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/SDGDLL_zone2.shp"), layer="SDGDLL_zone2") 
WCDGDLL=readOGR(handl_OneDrive("Data/Mapping/Shark_shape_files/WCDGDLL.shp"), layer="WCDGDLL") 



#Plot management zones and receiver lines for TDGDLF
library(PBSmapping)
#receivers
AATAMS.all=read.csv(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/AATAMS.all.csv"))
SMN.all=read.csv(handl_OneDrive("Data/Tagging/Acoustic_tagging/Acoustic_tagging_data/SMN.all.csv"))
Receivers=rbind(SMN.all,AATAMS.all)
Receivers$Station=paste(Receivers$latitude,Receivers$longitude) 
Receivers=Receivers[order(Receivers$Station),]
STATIONS=Receivers[!duplicated(Receivers$Station),]
STATIONS=subset(STATIONS,longitude<=150)
Perth=c(115.866,-31.95)
Rotnest=c(115.50,-32.02)

Zn.col=c("grey80","grey60","grey40")
names(Zn.col)=c("West","Zone1","Zone2" )

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Plot.Map.R"))  
data(worldLLhigh)

fn.map=function(South.WA.lat,South.WA.long,PLATE,dat)
{
  LATT=South.WA.lat[2]:South.WA.lat[1]
  LONGG=South.WA.long[1]:South.WA.long[2]
  S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
  S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
  
  plotmap(South.WA.long,South.WA.lat,PLATE,"grey95",XLIM=South.WA.long,YLIM=South.WA.lat)
  
  plot(WCDGDLL,add=T,col=Zn.col[1])
  plot(SDGDLL_zone1,add=T,col=Zn.col[2])
  plot(SDGDLL_zone2,add=T,col=Zn.col[3])
  text(114.8,-32.1,("1"),col="white", cex=3)
  text(114.8,-34,("2"),col="white", cex=3)
  text(118,-35.3,("3"),col="white", cex=3)
  
  axis(side = 1, at =seq(112,129,1), labels = F, tcl = .4)
  axis(side = 2, at =seq(-36,-25,1), labels = F, tcl = .4)
  axis(side = 1, at =seq(112,129,2), labels = seq(112,129,2), tcl = .7,las=1,cex.axis=1.5,padj=-.5)
  axis(side = 2, at = seq(-36,-25,2), labels = -seq(-36,-25,2), tcl = .7,las=2,cex.axis=1.5,hadj=.3)
  
  #Add tags
  dat$Rel.zone=as.character(with(dat,
                                 ifelse(Long.rels>=116.5 & Lat.rels<=(-26),"Zone2",
                                        ifelse(Long.rels<116.5 & Lat.rels<=(-33),"Zone1",
                                               ifelse(Lat.rels>(-33) & Lat.rels<=(-26) & Long.rels<116.5,"West",
                                                      ifelse(Lat.rels>(-26) & Long.rels<114,"West",
                                                             ifelse(Lat.rels>(-26) & Long.rels>=114 & Long.rels<123.75,"North",
                                                                    ifelse(Lat.rels>(-26) & Long.rels>=123.75,"Joint",NA))))))))
  dat$Rel.zone=with(dat,ifelse(Long.rels>=129 & Lat.rels<=(-26),"SA",Rel.zone))
  dat$Rec.zone=as.character(with(dat,
                                 ifelse(Long.rec>=116.5 & Lat.rec<=(-26),"Zone2",
                                        ifelse(Long.rec<116.5 & Lat.rec<=(-33),"Zone1",
                                               ifelse(Lat.rec>(-33) & Lat.rec<=(-26) & Long.rec<116.5,"West",
                                                      ifelse(Lat.rec>(-26) & Long.rec<114,"West",
                                                             ifelse(Lat.rec>(-26) & Long.rec>=114 & Long.rec<123.75,"North",
                                                                    ifelse(Lat.rec>(-26) & Long.rec>=123.75,"Joint",NA))))))))
  dat$Rec.zone=with(dat,ifelse(Long.rec>=129 & Lat.rec<=(-26),"SA",Rec.zone))
  
  Tab.rel=table(dat$Rel.zone)
  Tab.rel=log(Tab.rel)
  Tab.rec=table(dat$Rec.zone)
  Tab.rec=log(Tab.rec)
  
  La=c(-32.6,-35,-35.3)
  Ln=c(114.6,114.7,116.8)
  names(La)=names(Ln)=c("West","Zone1", "Zone2")
  Coord=data.frame(Zone=names(La),long=Ln,lat=La)
  plt.rel=data.frame(Zone=names(Tab.rel),pch=c(Tab.rel))
  plt.rel=merge(plt.rel,Coord,by="Zone",all.x=T)
  points(plt.rel$long,plt.rel$lat,cex=plt.rel$pch*scaler,pch=21,col="white",bg="black")
  
  plt.rec=data.frame(Zone=names(Tab.rec),pch=c(Tab.rec))
  plt.rec=merge(plt.rec,Coord,by="Zone",all.x=T)
  plt.rec$long=plt.rec$long+.4
  points(plt.rec$long,plt.rec$lat,cex=plt.rec$pch*scaler,pch=21,bg="white",col="black")
  box()
  
}

scaler=1.15

tiff(file="Fig1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(1,1,1,1),mgp=c(1,1,0))
#Plot conventional
fn.map(South.WA.long=c(114.4,118.8),South.WA.lat=c(-35.6,-31.8),PLATE=c(.01,.9,.075,.9),dat=Raw.Conv.Gummy)
mtext("Gummy shark",3,cex=1.75,line=1)

fn.map(South.WA.long=c(114.4,118.8),South.WA.lat=c(-35.6,-31.8),PLATE=c(.01,.9,.075,.9),dat=Raw.Conv.Whiskery)
mtext("Whiskery shark",3,cex=1.75,line=1)
mtext("Conventional tags",4,cex=1.75,line=1)
legend("topright",paste(c(5,"",25,"",50)),bty='n',pch=19,title="Number of tags",cex=1.5,
       col=c("black","transparent","black","transparent","black"),pt.cex=scaler*log(c(5,1,25,1,50)))

#Plot acoustic
fn.map(South.WA.long=c(114.4,118.8),South.WA.lat=c(-35.6,-31.8),PLATE=c(.01,.9,.075,.9),dat=Raw.Acous.Gummy)
points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19,cex=.8)  #receiver location

fn.map(South.WA.long=c(114.4,118.8),South.WA.lat=c(-35.6,-31.8),PLATE=c(.01,.9,.075,.9),dat=Raw.Acous.Whiskery)
points(STATIONS$longitude,STATIONS$latitude,col=1,pch=19,cex=.8)  #receiver location

mtext("Acoustic tags",4,cex=1.75,line=1)

text(117.25,-33,("Western"),col="black", cex=3)
text(117.25,-33.5,("Australia"),col="black", cex=3)

mtext("Latitude (?S)",side=2,line=-0.5,las=3,cex=1.75,outer=T)
mtext("Longitude (?E)",side=1,line=-0.5,cex=1.75,outer=T)

dev.off()



#1.2 plot data inputs
setwd(paste(hndl,"Show Gummy and whiskery outputs",sep="/"))

#Conventional tagging                   
fn.rec.ind.base=function(DAT,where,Ymax)
{
  dat=subset(DAT,Rel.zn==where)
  dat$N=1
  dat=aggregate(N~Rec.zn+DaysAtLarge,dat,sum)
  dummy=data.frame(CL=CL,Rec.zn=names(CL),jitr=jit.zn)
  dat=merge(dat,dummy,by="Rec.zn")
  dat$CL=as.character(dat$CL)
  with(dat,plot(DaysAtLarge,N*jitr,ylim=c(0,Ymax),xlim=c(0,max(DAT$DaysAtLarge)),
                pch=21,cex=2,ylab="",xlab="",cex.axis=2,bg=CL))
  legend("topright",paste("Released in",where,sep=" "),cex=1.8,bty='n')
}

Zns=unique(Conv.Whiskery$Rec.zn)
CL=c("black","grey60","white")
jit.zn=c(.9,1,1.1)
names(CL)=names(jit.zn)=c("West","Zone1","Zone2")

tiff(file="FigA1.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")  
par(mfcol=c(3,2),mai=c(.35,.4,.3,.1),oma=c(2,2,.01,.1),las=1,mgp=c(1.6,.75,0))

#Gummy
fn.rec.ind.base(Conv.Gummy,"West",Ymax=8)
legend("bottomright",c("West", "Zone 1", "Zone 2"),
       bty='n',pch=c(21,21,21),cex=1.65,pt.bg=CL,title="Recaptured in")
mtext("Gummy shark",3,0.5,cex=1.5)
fn.rec.ind.base(Conv.Gummy,"Zone1",Ymax=8)
fn.rec.ind.base(Conv.Gummy,"Zone2",Ymax=8)  

#Whiskery
fn.rec.ind.base(Conv.Whiskery,"West",Ymax=5)
mtext("Whiskery shark",3,0.5,cex=1.5)
fn.rec.ind.base(Conv.Whiskery,"Zone1",Ymax=5)
fn.rec.ind.base(Conv.Whiskery,"Zone2",Ymax=5)  
mtext("Days at large",1,0.5,outer=T,cex=1.75)
mtext("Numbers",2,-0.5,outer=T,las=3,cex=1.75)  

dev.off()


#Acoustic tagging             
#recaptures
fn.rec.ind.base.a=function(DAT,where,Ymax)
{
  dat=subset(DAT,Rel.zn==where)
  dat$DaysAtLarge=30+dat$DaysAtLarge
  YRS=1:3
  if(nrow(dat)==0)
  {
    plot(YRS[1]:YRS[length(YRS)],YRS[1]:YRS[length(YRS)],ylim=c(0,1),
         ylab="",xlab="",cex.axis=2,col='transparent',xaxt='n',yaxt='n')    
  }
  
  if(nrow(dat)>0)
  {
    dat$N=1
    dat=aggregate(N~Rec.zn+DaysAtLarge,dat,sum)
    dummy=data.frame(CL=CL,Rec.zn=names(CL),jitr=jit.zn)
    dat=merge(dat,dummy,by="Rec.zn")
    dat$CL=as.character(dat$CL)
    with(dat,plot(DaysAtLarge,N*jitr,ylim=c(0,Ymax),xlim=c(0,max(DAT$DaysAtLarge)),
                  pch=21,cex=2,ylab="",xlab="",cex.axis=2,bg=CL))
  }
  legend("topright",paste("Released in",where,sep=" "),cex=1.8,bty='n')
}

tiff(file="FigA2.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")  
par(mfcol=c(3,2),mai=c(.35,.4,.3,.1),oma=c(2,2,.01,.1),las=1,mgp=c(1.6,.75,0))

#Gummy
fn.rec.ind.base.a(Acous.Gummy,"West",Ymax=15)
mtext("Gummy shark",3,0.5,cex=1.5)
fn.rec.ind.base.a(Acous.Gummy,"Zone1",Ymax=15)
fn.rec.ind.base.a(Acous.Gummy,"Zone2",Ymax=15)  

#Whiskery
fn.rec.ind.base.a(Acous.Whiskery,"West",Ymax=5)
mtext("Whiskery shark",3,0.5,cex=1.5)
legend("bottomright",c("West", "Zone 1", "Zone 2"),
       bty='n',pch=c(21,21,21),cex=1.65,pt.bg=CL,title="Recaptured in")
fn.rec.ind.base.a(Acous.Whiskery,"Zone1",Ymax=5)
fn.rec.ind.base.a(Acous.Whiskery,"Zone2",Ymax=5)  
mtext("Days at large",1,0.5,outer=T,cex=1.75)
mtext("Numbers",2,-0.25,outer=T,las=3,cex=1.75)  

dev.off()


#Show acoustic tags detections by TagID
fn.hist.TAGID=function(dat)
{
  dat$TagID
  TLB=table(dat$TagID)
  a=barplot(TLB,col="grey60",names.arg=1:length(unique(dat$TagID)))
  axis(1,a,F)
  box()
}

tiff(file="Fig.A3.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")  
par(mfcol=c(2,1),mai=c(.6,.8,0.1,.1),oma=c(1.5,.1,1.5,.1),las=1,mgp=c(1,.65,0))

#Gummy
fn.hist.TAGID(Acous.Gummy)
mtext("Gummy shark",3,0.5,cex=1.5)

#Whiskery
fn.hist.TAGID(Acous.Whiskery)
mtext("Whiskery shark",3,0.5,cex=1.5)

mtext("Number of individuals",1,0,outer=T,cex=1.75)
mtext("Number of detection events",2,-2,outer=T,las=3,cex=1.75)  

dev.off()




#2.Display model outputs

#2.1 MLE table
a=Prob.mov.Gummy
names(a)=paste("Gummy",names(a),sep="_")
b=Prob.mov.Whiskery
names(b)=paste("Whiskery",names(b),sep="_")
TAbl1=merge(a,b,by.x="Gummy_name",by.y="Whiskery_name")
TAbl1=TAbl1[,-match(c("Gummy_index","Whiskery_index"),names(TAbl1))]
names(TAbl1)[1]="Parameter"
write.csv(TAbl1,"Table1.csv",row.names=F)

#2.2 Annual transition matrix

    #put parameters in matrix in normal space
p11=exp(Prob.mov.Gummy$value[1])
p21=exp(Prob.mov.Gummy$value[3])
p22=exp(Prob.mov.Gummy$value[2])
p33=exp(Prob.mov.Gummy$value[4])
p32=1-p33
p12=1-p11
p13=0
p23=1-(p21+p22)
p31=0
Mov.mat.Gummy=matrix(c(p11,p12,p13,p21,p22,p23,p31,p32,p33),ncol=3,byrow=T)
colnames(Mov.mat.Gummy)=rownames(Mov.mat.Gummy)=c("West","Zone 1","Zone 2")

p11=exp(Prob.mov.Whiskery$value[1])
p21=exp(Prob.mov.Whiskery$value[3])
p22=exp(Prob.mov.Whiskery$value[2])
p33=exp(Prob.mov.Whiskery$value[4])
p32=1-p33
p12=1-p11
p13=0
p23=1-(p21+p22)
p31=0
Mov.mat.Whiskery=matrix(c(p11,p12,p13,p21,p22,p23,p31,p32,p33),ncol=3,byrow=T)
colnames(Mov.mat.Whiskery)=rownames(Mov.mat.Whiskery)=rownames(Mov.mat.Gummy)

#calculate annual matrix
Mov.mat.Gummy=Mov.mat.Gummy%^%365
Mov.mat.Whiskery=Mov.mat.Whiskery%^%365

#plot model outputs
#colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
colfunc <- colorRampPalette(c("grey10", "grey60","white"))
fn.show.MTM=function(MATRIZ,N.int,TITLE,add.legend)
{
  couleurs=rev(colfunc(N.int))
  BREAKS=seq(0,1,length.out=N.int+1)
  xx=1:nrow(MATRIZ)
  yy=1:ncol(MATRIZ)
  image(xx,yy,t(MATRIZ),ylab="",xlab="",xaxt='n',yaxt='n',col =couleurs,breaks=BREAKS)
  axis(1,xx,F,tck=0.025)
  axis(2,yy,F,tck=0.025)
  mtext(TITLE,3,.25,cex=1.5)
  box()
  SQ=rev(seq(BREAKS[1],BREAKS[length(BREAKS)],.25))
  if(add.legend=="YES")color.legend(3,1,3.35,2,SQ,rect.col=rev(couleurs),gradient="y",col=1,cex=1)
  axis(1,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],las=2,cex.axis=1.25,hadj=0.75,tck=0.025)
  axis(2,1:nrow(MATRIZ),colnames(MATRIZ)[1:nrow(MATRIZ)],cex.axis=1.25,tck=0.025,hadj=.75,las=1) 
}
tiff(file="Movement matrix.tiff",width = 1200, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,1),mai=c(.8,0.85,0.3,.1),oma=c(1.5,.1,.1,.1),las=1)
fn.show.MTM(MATRIZ=Mov.mat.Gummy,N.int=50,"Gummy shark",'NO')
fn.show.MTM(MATRIZ=Mov.mat.Whiskery,N.int=50,"Whiskery shark",'YES')
mtext("To",1,0.25,cex=1.85,outer=T)
mtext("From",2,-1.5,cex=1.85,outer=T,las=3) 
dev.off()

