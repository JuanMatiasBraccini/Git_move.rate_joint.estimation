#Movement model conditioned on recaptures only (McGarvey line of thought)

#note: individual-based model that compares the probability of occurrying in a particular zone
#      after an given time at liberty with the observed recapture zone.

rm(list=ls(all=TRUE))
library(expm)
library(ggplot2)
library(reshape)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Population dynamics/fn.fig.R"))


#----DATA SECTION ----

#number of sharks released per zone in year 1 (for conventional and acoustic tags)
N.sharks_conv=1000  
N.sharks_acous=100  

#note: all sharks released in year 1

#years with recaptures/detections 
Yrs=1:3
Yr.rec_daily=365*Yrs   # in days

#number of samples taken (observed recaptures or detections)
n.obs.conv=200
n.obs.acous=250


#----PARAMETERS SECTION ----

#movement probabilities (for time step of 1 day)
From_1_to_1=0.9995
From_1_to_2=1-From_1_to_1
move.from.1=c(From_1_to_1,From_1_to_2,0)  #assumption of no movement to non-adjacent zone in 1 day

From_2_to_2=0.9995
From_2_to_1=0.00025
move.from.2=c(From_2_to_1,From_2_to_2,1-(From_2_to_1+From_2_to_2))


From_3_to_3=0.9995
From_3_to_2=1-From_3_to_3
move.from.3=c(0,From_3_to_2,From_3_to_3)  #assumption of no movement to non-adjacent zone in 1 day

Nzone=3

Pin.pars=c(From_1_to_1,From_2_to_1,From_2_to_2,From_3_to_3)
names(Pin.pars)=c("P11","P21","P22","P33")

#recapture/detection rate
#note: this captures mortality and recapture/detection probabilities
see.rate_conv=.7  # rate of recapturing shark with conventional tags
see.rate_acous=.9 #lower  rate for acoustic due to just mortality

#number of simulations 
n.sims=1000


#----PROCEDURE SECTION ----

#population movement transition matrix, daily time step, i.e. probability of movement in 1 day (rows sum to 1)
Mov.Mat=matrix(c(move.from.1,move.from.2,move.from.3),ncol=Nzone,nrow=Nzone,byrow=T)   

#Functions

fn.logit=function(x) log(x/(1-x))
fn.inv.logit=function(x) exp(x)/(1+exp(x))

posfun=function(x,eps)
{
  if(x>=eps) return (x)
  if(x<eps) 
  {
    fpen<<-fpen+0.01*(x-eps)^2
    return(eps/(2-(x/eps)))
  }
  
}

#simulate sampling (recaptures / detections)
fn.sample=function(D,nObs)
{
  From=1:nrow(D[[1]])
  To=1:ncol(D[[1]])
  From_to=expand.grid(From,To)
  From_to=From_to[order(From_to$Var1),]
  From_to=paste(From_to$Var1,From_to$Var2,sep="_")
  
  Dat=D
  for(i in 1:length(D))
  {
    x=D[[i]]
    Obs=c(t(x))
    names(Obs)=From_to
    Obs=rep(names(Obs),Obs)
    Dat[[i]]=data.frame(Days.liberty=as.numeric(names(D)[i]),From=substr(Obs,1,1),To=substr(Obs,3,3))
  }
  Dat=do.call(rbind,Dat)
  
  index=sample(1:nrow(Dat),nObs)    #simulate sampling of tagged population
  
  return(Dat[index,])
}

#likelihood
fn.like=function(MAT,DATA)
{
  NLL=0
  for(i in 1:nrow(DATA))
  {
    d=DATA[i,]
    Move=MAT %^% d$Days.liberty
    id.rel=which(rownames(Move)==d$From)
    id.rec=which(colnames(Move)==d$To)
    Pred.Prob=Move[id.rel,id.rec]
    NLL=NLL-log(Pred.Prob)
  }
  return(NLL)
}

#individual-based model
fn.move=function(pars)
{
  #Put movement pars in matrix
  P.1=fn.inv.logit(pars[1])   #inverse logit
  P.2=fn.inv.logit(pars[2:3])
  P.3=fn.inv.logit(pars[4])
  
  #get missing parameter              
  P.1=c(P.1,1-P.1,0)    
  P.2=c(P.2,1-sum(P.2))
  P.3=c(0,1-P.3,P.3)     
  
  #put pars in matrix
  Mov.mat=as.matrix(rbind(P.1,P.2,P.3))
  colnames(Mov.mat)=rownames(Mov.mat)=1:nrow(Mov.mat)
  
  #Calculate likelihood of convetional tagging observations
  neg.LL.conv=0
  neg.LL.conv=fn.like(MAT=Mov.mat,DATA=Data_conv)
  
  #Calculate likelihood of acoustic tagging observations
  neg.LL.acous=0
  neg.LL.acous=fn.like(MAT=Mov.mat,DATA=Data_acous)
  
  #total neg log like
  neg.LL=neg.LL.conv+neg.LL.acous
  
  return(neg.LL)
}



#Run simulation testing
Store.sims=vector('list',n.sims)
system.time(for(n in 1:n.sims)   #takes 37 seconds per iteration
{
  #1. Generate recaptures/detections based on movement matrix and number of sharks released
  
  #conventional tagging
  Rec.list=vector('list',length(Yr.rec_daily))
  names(Rec.list)=Yr.rec_daily
  Rec.list[[1]]=floor(N.sharks_conv*(Mov.Mat%^% Yr.rec_daily[1])*see.rate_conv) #release (in day 1) and move sharks in year 1
  for(i in 2:length(Yr.rec_daily))  #subsequent years
  {
    Move=Mov.Mat %^% 365      #move shark one more year
    Rec.list[[i]]=floor((Rec.list[[i-1]]%*%Move)*see.rate_conv)
  }
  
  #acoustic tagging
  Detec.list=vector('list',length(Yr.rec_daily))
  names(Detec.list)=Yr.rec_daily
  Detec.list[[1]]=floor(N.sharks_acous*(Mov.Mat%^% Yr.rec_daily[1])*see.rate_acous)
  for(i in 2:length(Yr.rec_daily))
  {
    Move=Mov.Mat %^% 365     #move sharks one more year
    Detec.list[[i]]=floor((Detec.list[[i-1]]%*%Move)*see.rate_acous)
  }
  
  
  #2. Create observations
  #conventional tagging
  Data_conv=fn.sample(D=Rec.list,nObs=n.obs.conv)   
  
  #acoustic tagging
  Data_acous=fn.sample(D=Detec.list,nObs=n.obs.acous)   
  #Data_acous$From_to=with(Data_acous,paste(From,To,sep="_")); with(Data_acous,table(Days.liberty,From_to))
  
  
  #3. Estimate parameters
  Pars= c(0.99, 0.005, 0.99, 0.99)
  Pars=fn.logit(Pars)
  fit=optim(Pars,fn.move,method="Nelder-Mead", control = list(trace=T))
  
  #4. Store pars
  Store.sims[[n]]=fn.inv.logit(fit$par) 
})

Store.sims=do.call(rbind,Store.sims)


#Plot simulations and original par values
Store.pin.pars=matrix(rep(Pin.pars,n.sims),nrow=n.sims,byrow=T)
ERROR=Store.sims
#ERROR[,match(c("P11","P22","P33"),names(Pin.pars))]=round(ERROR[,match(c("P11","P22","P33"),names(Pin.pars))],6)
#ERROR[,match(c("P21"),names(Pin.pars))]=round(ERROR[,match(c("P21"),names(Pin.pars))],6)
ERROR=Store.pin.pars-ERROR
#ERROR=ERROR/Store.pin.pars
colnames(ERROR)=names(Pin.pars)
setwd(handl_OneDrive("Analyses/Movement rate estimation/Joint.estim_ind.base.mod/Simulation_testing"))
Do.jpeg="YES"
Do.tiff="NO"

fn.fig("boxplot",2400, 1600)
par(las=1)
p=ggplot(data=melt(ERROR), aes(as.factor(X2), value)) + ylim(-4e-4,4e-4)
p=p+ geom_violin(fill="grey80",adjust=2) + geom_jitter(height = 0,width = 0.15,color="grey20")
p=p+labs(y="Error",x="Estimated parameter")+ geom_hline(aes(yintercept=0),lty=2,colour="grey50")
p+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        panel.background = element_rect(fill = "white", colour = "grey50"))
dev.off()

XLI=list(c(.99875,1.000075),c(-3e-4,9e-4),c(.99875,1.000075),c(.99875,1.000075))
fn.fig("density",2400, 2400)
par(mfcol=c(2,2),mai=c(.65,.85,.1,.1),oma=c(.1,.1,.1,.1),las=1,mgp=c(2,.6,0))
for(i in 1:ncol(Store.sims))
{
  plot(density(Store.sims[,i],adjust=10),type='l',lwd=2,ylab="",xlab=names(Pin.pars)[i],main="",
       cex.lab=1.5,cex.axis=1.25,xlim=XLI[[i]])
  abline(v=Pin.pars[i],lwd=2,lty=2,col="grey60")
}
mtext("Density",2,outer=T,las=3,line=-2,cex=2)
dev.off()

