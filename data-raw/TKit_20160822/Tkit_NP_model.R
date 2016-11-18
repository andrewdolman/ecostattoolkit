#Geoff Phillips
#Date:12 August 2016
#File name for script: Tkit_NP_model.R

rm(list= ls())  		# Clear data

library (lmodel2)
library(car)

#Enter file name and path
FName<-"O:\\Work\\Projects\\2016\\JRC_NutrientPhase2\\Data\\Analysis\\2016_08_10\\DataTemplate_L-CB2L_EC1.csv"

FName<-"DataTemplate_L-CB2L_EC1.csv"


data <- read.csv(file = FName, header = TRUE)
dim(data)
names(data)
summary(data)

GM<-0.6
HG<-0.8

class.col<-c("red","orange","yellow","green","blue")	#set up classification colours
ymin<-0
ymax<-1.5

#Get all data with N&P for phytoplankton

cc<-complete.cases(data$N,data$P,data$EQR)
data.cc<-data[cc,]
dim(data.cc)
summary(data.cc)

y<-data.cc$EQR
x1<-log10(data.cc$P)
x2<-log10(data.cc$N)

Class<-as.factor(data.cc$BioClass)
x1.2<-10^(x1)
x2.2<-10^(x2)

# Consider linear portions, although for this model
ymin<-P.minUsed<-1
ymax<-P.maxUsed<-5000
xmin<-N.minUsed<-100
xmax<-N.maxUsed<-50000


y.u<-data.cc$EQR[data.cc$P>=P.minUsed & data.cc$P<=P.maxUsed & data.cc$N<=N.maxUsed 
	& data.cc$Exclude_N==FALSE & data.cc$Exclude_P==FALSE]
x1.u<-log10(data.cc$P[data.cc$P>=P.minUsed & data.cc$P<=P.maxUsed & data.cc$N<=N.maxUsed
	& data.cc$Exclude_N==FALSE & data.cc$Exclude_P==FALSE])
x2.u<-log10(data.cc$N[data.cc$P>=P.minUsed & data.cc$P<=P.maxUsed & data.cc$N<=N.maxUsed
	& data.cc$Exclude_N==FALSE & data.cc$Exclude_P==FALSE])
Class.u<-as.factor(data.cc$BioClass[data.cc$P>=P.minUsed & data.cc$P<=P.maxUsed & data.cc$N<=N.maxUsed
	& data.cc$Exclude_N==FALSE & data.cc$Exclude_P==FALSE])
x1.2.u<-10^(x1.u)
x2.2.u<-10^(x2.u)


length(y.u)
length(x1.u)

summary(NP.mod0<-lm(y.u~x1.u+x2.u))
summary(P.mod0<-lm(y.u~x1.u))
summary(N.mod0<-lm(y.u~x2.u))

AIC(N.mod0,NP.mod0,P.mod0)

#Check model
sqrt(vif(NP.mod0))#Check colinearity >2 a problem ?
plot(NP.mod0)

#------------------------------------------------------------------------------------------
#	Fit RML regression of TP v TN concentration
#------------------------------------------------------------------------------------------

lmod<-lmodel2(x1~x2,range.y="relative",range.x="interval",nperm=99,
	data=data.cc)
lmod						#Summary of all  models including conventional OLS
(RMA.int<-lmod$regression.results[[2]][[4]])
(RMA.slope<-lmod$regression.results[[3]][[4]])	




#==========================================================================================
#	Plot TP v TN and add contour of fitted line for GM & HG boundaries
#==========================================================================================
#------------------------------------------------------------------------------------------
#	Extract model parameters for use in plotting
#------------------------------------------------------------------------------------------
(int.NP<-coef(NP.mod0)[[1]])		#extract coefficients of model  HA lakes
(slope.P<-coef(NP.mod0)[[2]])		#use [[]] to extract numeric part of named number
(slope.N<-coef(NP.mod0)[[3]])

Mod.NP.resid<-resid(NP.mod0)	#Extract residuals of predicted EQRs
Mod.NP.up<-quantile(Mod.NP.resid,0.75)	#Calculate upper 75th quantile of residuals
Mod.NP.lw<-quantile(Mod.NP.resid,0.25)	#Calculate lower 25th quantile of residuals

#------------------------------------------------------------------------------------------
#	Create plots 
#------------------------------------------------------------------------------------------
ylb<-expression(paste("total phosphorus ","(µg ", L^-1,")"))
xlb<-expression(paste("total nitrogen ","(µg ", L^-1,")"))


win.graph(width=12,height=4.5)
par(mfrow=c(1,2))
par(mar=c(4, 4.5, 2.5, 2) + 0.1)

plot(x1.2 ~ x2.2,log="xy",ylab=ylb,xlab=xlb,las=1,
	ylim=c(ymin,ymax),xlim=c(xmin,xmax))
points(x1.2.u ~x2.2.u,	col=class.col[Class.u],pch=19)
#abline(-1,1,lty=3)
abline(-1.176,1,lty=4,col="orange")

#points(P~ N,data=data.cc,pch=1,subset=Group=="L-EC1")	
abline(RMA.int,RMA.slope,lty=2)

#------------------------------------------------------------------------------------------
#	Add contour line of TN and TP at GM boundary +/- residuals bounding 50% data
#------------------------------------------------------------------------------------------
mtext("Good/Moderate")
new.dat<-data.frame(P = seq(ymin, ymax, length.out=20))		#Edit values if needed
new.dat$N.GM<-10^((GM-int.NP-slope.P*log10(new.dat$P))/slope.N)
new.dat$N.GM.up<-10^((GM-int.NP-Mod.NP.up-slope.P*log10(new.dat$P))/slope.N)
new.dat$N.GM.lw<-10^((GM-int.NP-Mod.NP.lw-slope.P*log10(new.dat$P))/slope.N)
lines(new.dat$N.GM,new.dat$P,lty=1,col="green")
lines(new.dat$N.GM.up,new.dat$P,lty=2,col="green")
lines(new.dat$N.GM.lw,new.dat$P,lty=2,col="green")
#------------------------------------------------------------------------------------------
#	Calculate intersection of lines and mark on graph
#------------------------------------------------------------------------------------------ 
LN.GM<-(GM-(slope.P*RMA.int)-int.NP)/((RMA.slope*slope.P)+slope.N)
LP.GM<-LN.GM*RMA.slope + RMA.int
LN.GM.L<-(GM-(slope.P*RMA.int)-int.NP-Mod.NP.lw)/((RMA.slope*slope.P)+slope.N)
LP.GM.L<-LN.GM.L*RMA.slope + RMA.int
LN.GM.U<-(GM-(slope.P*RMA.int)-int.NP-Mod.NP.up)/((RMA.slope*slope.P)+slope.N)
LP.GM.U<-LN.GM.U*RMA.slope + RMA.int

#P Boundary values for GM
(P.GM<-round((10^LP.GM),0))
(P.GM.L<-round((10^LP.GM.L),0))
(P.GM.U<-round((10^LP.GM.U),0))
abline("h"=c(P.GM.L,P.GM,P.GM.U),lty=3)
text(xmin,P.GM,P.GM,pos=4,cex=0.8)
text(xmin,P.GM.L,P.GM.L,pos=4,cex=0.8)
text(xmin,P.GM.U,P.GM.U,pos=4,cex=0.8)

#N Boundary values for GM
(N.GM<-round((10^LN.GM),0))
(N.GM.L<-round((10^LN.GM.L),0))
(N.GM.U<-round((10^LN.GM.U),0))
abline("v"=c(N.GM.L,N.GM,N.GM.U),lty=3)
text(N.GM,ymin,N.GM,pos=3,cex=0.8)
text(N.GM.L,ymin,N.GM.L,pos=2,cex=0.8)
text(N.GM.U,ymin,N.GM.U,pos=4,cex=0.8)

a<-N.GM/P.GM
#------------------------------------------------------------------------------------------
#	Calculate predicted TN for given TP at EQR H/G boundary and add contour line
#------------------------------------------------------------------------------------------
plot(x1.2 ~ x2.2,log="xy",ylab=ylb,xlab=xlb,las=1,
	ylim=c(ymin,ymax),xlim=c(xmin,xmax))
points(x1.2.u ~x2.2.u,	col=class.col[Class.u],pch=19)
#abline(-1,1,lty=3)
abline(-1.176,1,lty=4,col="orange")

#points(P~ N,data=data.cc,pch=1,subset=Group=="L-EC1")	
abline(RMA.int,RMA.slope,lty=2)
mtext("High/Good")
abline(RMA.int,RMA.slope,lty=2)
new.dat<-data.frame(P = seq(ymin, ymax, length.out=20))	
new.dat$N.HG<-10^((HG-int.NP-slope.P*log10(new.dat$P))/slope.N)
new.dat$N.HG.up<-10^((HG-int.NP-Mod.NP.up-slope.P*log10(new.dat$P))/slope.N)
new.dat$N.HG.lw<-10^((HG-int.NP-Mod.NP.lw-slope.P*log10(new.dat$P))/slope.N)
lines(new.dat$N.HG,new.dat$P,lty=1,col="blue")
lines(new.dat$N.HG.up,new.dat$P,lty=2,col="blue")
lines(new.dat$N.HG.lw,new.dat$P,lty=2,col="blue")
#------------------------------------------------------------------------------------------
#	Calculate intersection of HG lines and mark on graph
#------------------------------------------------------------------------------------------ 
LN.HG<-(HG-(slope.P*RMA.int)-int.NP)/((RMA.slope*slope.P)+slope.N)
LP.HG<-LN.HG*RMA.slope + RMA.int
LN.HG.L<-(HG-(slope.P*RMA.int)-int.NP-Mod.NP.lw)/((RMA.slope*slope.P)+slope.N)
LP.HG.L<-LN.HG.L*RMA.slope + RMA.int
LN.HG.U<-(HG-(slope.P*RMA.int)-int.NP-Mod.NP.up)/((RMA.slope*slope.P)+slope.N)
LP.HG.U<-LN.HG.U*RMA.slope + RMA.int

#P Boundary values for HG
(P.HG<-round((10^LP.HG),0))
(P.HG.L<-round((10^LP.HG.L),0))
(P.HG.U<-round((10^LP.HG.U),0))
abline("h"=c(P.HG.L,P.HG,P.HG.U),lty=3)
text(xmin,P.HG,P.HG,pos=4,cex=0.8)
text(xmin,P.HG.L,P.HG.L,pos=4,cex=0.8)
text(xmin,P.HG.U,P.HG.U,pos=4,cex=0.8)

#N Boundary values for HG
(N.HG<-round((10^LN.HG),0))
(N.HG.L<-round((10^LN.HG.L),0))
(N.HG.U<-round((10^LN.HG.U),0))
abline("v"=c(N.HG.L,N.HG,N.HG.U),lty=3)
text(N.HG,ymin,N.HG,pos=3,cex=0.8)
text(N.HG.L,ymin,N.HG.L,pos=2,cex=0.8)
text(N.HG.U,ymin,N.HG.U,pos=4,cex=0.8)



#=====================================================================================
###########################################################################################
#Output boundary values, 
(out1<-data.frame(
	R2=c(round(summary(NP.mod0)[[9]],3),""),
	N=c(length(x1.u),""),
	Nut=c("P","N"),
	slo=c(round(slope.P,3),round(slope.N,3)),
	int=c(round(int.NP,3),""),
	GM =c(P.GM,N.GM),
	GML=c(P.GM.L,N.GM.L),
	GMU=c(P.GM.U,N.GM.U),
	HG =c(P.HG,N.HG),
	HGL=c(P.HG.L,N.HG.L),
	HGU=c(P.HG.U,N.HG.U)))
write.table(out1,"clipboard",sep="\t",row.names=F,col.names=F)
