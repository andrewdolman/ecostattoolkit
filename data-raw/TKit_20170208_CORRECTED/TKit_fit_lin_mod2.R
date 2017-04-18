#Description: Example script to fit linear models to data of sevoral categories
#Exploratory version
#Geoff Phillips
#Date:10 Aug 2016
#File name for script: TKit_fit_lin_mod2.R

rm(list= ls())  		# Clear data

#################################################################################
#		Step 1 get the data		 	 					  #
#################################################################################

#Enter file name and path
FName<-"O:\\OneDrive\\Work\\Projects\\2016\\JRC_NutrientPhase2\\Data\\Analysis\\2016_08_10\\DataTemplate_Example2.csv"

#Get data and check
data <- read.csv(file = FName, header = TRUE)
dim(data)
summary(data)

#Identify the records with both EQR and TP data to produce a data set without missing values
cc<-complete.cases(data$EQR,data$P)
data.cc<-data[cc,]
dim(data.cc)

#------------------------------------------------------------------------------------------
#			Get all data for plotting
#------------------------------------------------------------------------------------------
x<-log10(data.cc$P)
x2<-10^x
y<-data.cc$EQR
z<-as.factor(data.cc$Group)
#------------------------------------------------------------------------------------------
#             Assign data used for models to x.u and y.u
#------------------------------------------------------------------------------------------
nut.min<-7	#min nutrient value used for model
nut.max<-100	#max nutrient value used for model 

# Enter the high/good and good/moderate boundary values
HG<-0.80
GM<-0.60

x.u<-log10(data.cc$P[data.cc$P<=nut.max &
	data.cc$P>=nut.min &
	data.cc$Exclude_P==FALSE])
y.u<-data.cc$EQR[data.cc$P<=nut.max & 
	data.cc$P>=nut.min &
	data.cc$Exclude_P==FALSE]
x2.u<-10^x.u  #back transformed x for plotting on linear scale
z.u<-as.factor(data.cc$Group[data.cc$P<=nut.max & 
	data.cc$P>=nut.min &
	data.cc$Exclude_P==FALSE])
zN.u<-length(levels(z.u))
# Check the data records are the same length
length(x.u)
length(y.u)
length(z.u)
#..........................................................................................
# Plot the data to check and visualise
#..........................................................................................

MetricUsed<-"Phytoplankton"
#	labels for axis
ylb<-"EQR" # Y axis label
xlb<-expression(paste("total phosphorus ","(?g ", L^-1,")"))

plot(y~x2,main=MetricUsed,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
legend("bottomleft",levels(z.u),col=rainbow(zN.u),pch=19,cex=0.8)

#..........................................................................................
# Compare models including grouping variable, with and without interaction
#..........................................................................................
summary(mod1<-lm(y.u~x.u))	 #base model, no group variable
summary(mod1a<-lm(y.u~x.u+z.u))#model with group variable, different intercepts
summary(mod1b<-lm(y.u~x.u*z.u))#model with interaction term, different slopes & intercept
AIC(mod1,mod1a,mod1b)		 #compare models, select model with lowest AIC

#------------------------------------------------------------------------------------------
# Fit models including group variable if AIC mod1a is lowest
#------------------------------------------------------------------------------------------

#..........................................................................................
#	Model 1a OLS Y on X + Z
#..........................................................................................
summary(mod1a<-lm(y.u~x.u+z.u))

# 	predict values of y from model 1 
pred.y1a<-data.frame(predict(mod1a),resid(mod1a),z.u)#Calculate predicted values and residuals

#	Calculate upper and lower quartiles of residuals for 1st factor
upresid.1a<-quantile(subset(pred.y1a$resid.mod1a,z.u==levels(z.u)[1]),0.75,na.rm=T)
lwresid.1a<-quantile(subset(pred.y1a$resid.mod1a,z.u==levels(z.u)[1]),0.25,na.rm=T)

#	Extract the slope for model 1a
slope1a<-summary(mod1a)[[4]][[2]]

#	Extract the intercept for 1st grouping factor of model 1a		
int1a<-summary(mod1a)[[4]][[1]]

#	Calculate GM and HG boundaries with upper and lower estimates for 1st factor	
GM.mod1a<-round(10^((GM-int1a[1])/slope1a),0)
GMU.mod1a<-round(10^((GM-int1a[1]-upresid.1a[1])/slope1a),0)
GML.mod1a<-round(10^((GM-int1a[1]-lwresid.1a[1])/slope1a),0)
HG.mod1a<-round(10^((HG-int1a[1])/slope1a),0)
HGU.mod1a<-round(10^((HG-int1a[1]-upresid.1a[1])/slope1a),0)
HGL.mod1a<-round(10^((HG-int1a[1]-lwresid.1a[1])/slope1a),0)

#	Extract quartiles of residuals, intercepts, calculate boundaries for
#	remaining factors on model 1a
for (i in 2:zN.u){
	upresid.1a<-cbind(upresid.1a,quantile(subset(pred.y1a$resid.mod1a,z.u==levels(z.u)[i]),0.75,na.rm=T))
	lwresid.1a<-cbind(lwresid.1a,quantile(subset(pred.y1a$resid.mod1a,z.u==levels(z.u)[i]),0.25,na.rm=T))	
	int1a<-cbind(int1a,(int1a[1]+summary(mod1a)[[4]][[(i+1)]]))
	GM.mod1a<-cbind(GM.mod1a,round(10^((GM-int1a[i])/slope1a),0))
	GMU.mod1a<-cbind(GMU.mod1a,round(10^((GM-int1a[i]-upresid.1a[i])/slope1a),0))
	GML.mod1a<-cbind(GML.mod1a,round(10^((GM-int1a[i]-lwresid.1a[i])/slope1a),0))
	HG.mod1a<-cbind(HG.mod1a,round(10^((HG-int1a[i])/slope1a),0))
	HGU.mod1a<-cbind(HGU.mod1a,round(10^((HG-int1a[i]-upresid.1a[i])/slope1a),0))
	HGL.mod1a<-cbind(HGL.mod1a,round(10^((HG-int1a[i]-lwresid.1a[i])/slope1a),0))
	}

(Sum.mod1a<-data.frame(
	R2=c(round(summary(mod1a)[[9]],3)),
	N=c(length(x.u)),
	slope=c(slope1a),
	int=c(int1a),
	GM =c(GM.mod1a),
	GML=c(GML.mod1a),
	GMU=c(GMU.mod1a),
	HG =c(HG.mod1a),
	HGL=c(HGL.mod1a),
	HGU=c(HGU.mod1a),
	row.names=levels(z.u)))

#..........................................................................................
#	Model 2a OLS X on Y + Z
#..........................................................................................
summary(mod2a<-lm(x.u~y.u+z.u))

#	Extract the slope & intercept for 1st factor for model 2a
slope2a<-summary(mod2a)[[4]][[2]]
int2a<-summary(mod2a)[[4]][[1]]
y0<-(0-int2a)/slope2a

for (i in 2:zN.u){
	int2a<-cbind(int2a,(int2a[1]+summary(mod2a)[[4]][[(i+1)]]))
	y0<-cbind(y0,(0-int2a[i])/slope2a)
	}
# Set up a data frame with observed x, factor z and add appropriate intercept value y0a for category
# calculate predicted y values (EQR) and calculate residuals

(pred.y2a<-data.frame(x.u,y.u,z.u,y0[1]))
colnames(pred.y2a)[4]<-"y0.i"
for (i in 2:zN.u){
	pred.y2a<-within(pred.y2a,y0.i[z.u==levels(z.u)[i]]<-y0[i])
	}
(pred.y2a<-within(pred.y2a,pred.y2a<-x.u*1/slope2a+y0.i))	#calc predicted values pred.y2
(pred.y2a<-within(pred.y2a,resid.2a<-y.u-pred.y2a))		#calc residuals resid.2

# 	calculate upper and lower 25th 75th quantiles of residuals
upresid.2a<-quantile(subset(pred.y2a$resid.2a,z.u==levels(z.u)[1]),0.75,na.rm=T)
lwresid.2a<-quantile(subset(pred.y2a$resid.2a,z.u==levels(z.u)[1]),0.25,na.rm=T)
for (i in 2:zN.u){
	upresid.2a<-cbind(upresid.2a,quantile(subset(pred.y2a$resid.2a,z.u==levels(z.u)[i]),0.75,na.rm=T))
	lwresid.2a<-cbind(lwresid.2a,quantile(subset(pred.y2a$resid.2a,z.u==levels(z.u)[i]),0.25,na.rm=T))
	}

(GM.mod2a<-round(10^((GM-y0)*slope2a),0))			#Predicted GM boundary value for 1st factor (HA)
(GMU.mod2a<-round(10^((GM-y0-upresid.2a)*slope2a),0))
(GML.mod2a<-round(10^((GM-y0-lwresid.2a)*slope2a),0))
(HG.mod2a<-round(10^((HG-y0)*slope2a),0))			#Predicted GM boundary value for 1st factor (HA)
(HGU.mod2a<-round(10^((HG-y0-upresid.2a)*slope2a),0))
(HGL.mod2a<-round(10^((HG-y0-lwresid.2a)*slope2a),0))

(Sum.mod2a<-data.frame(
	R2=c(round(summary(mod2a)[[9]],3)),
	N=c(length(x.u)),
	slope=c(1/slope2a),
	int=c(y0),
	GM =c(GM.mod2a),
	GML=c(GML.mod2a),
	GMU=c(GMU.mod2a),
	HG =c(HG.mod2a),
	HGL=c(HGL.mod2a),
	HGU=c(HGU.mod2a),
	row.names=levels(z.u)))

#..........................................................................................
#	Model 3a Standard Major Axis (SMA) regression (geometric average of models 1 & 2)
#..........................................................................................
# calculate the average slope and intercept of x on y and y on x

(slope3a<-sign(slope1a)*sqrt(slope1a*1/slope2a))#geometric average or SMA regression
int3a<-mean(y.u[z.u==levels(z.u)[1]])-(slope3a*mean(x.u[z.u==levels(z.u)[1]]))
for (i in 2:zN.u){
	int3a<-cbind(int3a,mean(y.u[z.u==levels(z.u)[i]])-(slope3a*mean(x.u[z.u==levels(z.u)[i]])))
	}
# Set up a data frame with observed x, factor z and add appropriate intercept value y0a for category
# calculate predicted y values (EQR) and calculate residuals

(pred.y3a<-data.frame(x.u,y.u,z.u,int3a[1]))
colnames(pred.y3a)[4]<-"int3a.i"
for (i in 2:zN.u){
	pred.y3a<-within(pred.y3a,int3a.i[z.u==levels(z.u)[i]]<-int3a[i])
	}
(pred.y3a<-within(pred.y3a,pred.y3a<-x.u*slope3a+int3a.i))	#calc predicted values pred.y2
(pred.y3a<-within(pred.y3a,resid.3a<-y.u-pred.y3a))
	
# 	calculate upper and lower 25th 75th quantiles of residuals
upresid.3a<-quantile(subset(pred.y3a$resid.3a,z.u==levels(z.u)[1]),0.75,na.rm=T)
lwresid.3a<-quantile(subset(pred.y3a$resid.3a,z.u==levels(z.u)[1]),0.25,na.rm=T)
for (i in 2:zN.u){
	upresid.3a<-cbind(upresid.3a,quantile(subset(pred.y3a$resid.3a,z.u==levels(z.u)[i]),0.75,na.rm=T))
	lwresid.3a<-cbind(lwresid.3a,quantile(subset(pred.y3a$resid.3a,z.u==levels(z.u)[i]),0.25,na.rm=T))
	}

(GM.mod3a<-round(10^((GM-int3a)/slope3a),0))			#Predicted GM boundary value for 1st factor (HA)
(GMU.mod3a<-round(10^((GM-int3a-upresid.3a)/slope3a),0))
(GML.mod3a<-round(10^((GM-int3a-lwresid.3a)/slope3a),0))
(HG.mod3a<-round(10^((HG-int3a)/slope3a),0))			#Predicted GM boundary value for 1st factor (HA)
(HGU.mod3a<-round(10^((HG-int3a-upresid.3a)/slope3a),0))
(HGL.mod3a<-round(10^((HG-int3a-lwresid.3a)/slope3a),0))

(Sum.mod3a<-data.frame(
	R2="",
	N=c(length(x.u)),
	slope=c(slope3a),
	int=c(int3a),
	GM =c(GM.mod3a),
	GML=c(GML.mod3a),
	GMU=c(GMU.mod3a),
	HG =c(HG.mod3a),
	HGL=c(HGL.mod3a),
	HGU=c(HGU.mod3a),
	row.names=levels(z.u)))

#===========================================================================================
#	re-set scales for plotting
xmin<- 5
xmax<- 200
ymin<- 0
ymax<- 1.5

MetricUsed<-"Phytoplankton"
#	set up scales & labels for axis
ylb<-"EQR" # Y axis label
xlb<-expression(paste("total phosphorus ","(?g ", L^-1,")"))

win.graph(width=14,height=9)					# new graphic window
par(mfrow=c(2,3)) 						# delete this line if separate plots are required
par(mar=c(4,5,3.5,1))
#..........................................................................................
#	Fig A
#..........................................................................................
title<-"Regression of EQR on TP"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*slope1a+int1a[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}
legend("bottomleft",levels(z.u),col=rainbow(zN.u),pch=19,cex=1)
	mtext(paste("R2=", round(summary(mod1a)$r.sq,3),sep=""), side = 3, line=0,adj=0, 
	col="black",cex=0.8)
	pvalue<-pf(summary(mod1a)[[10]][[1]],summary(mod1a)[[10]][[2]],summary(mod1a)[[10]][[3]],lower.tail=F)
	if(pvalue >=0.001) {mtext(paste("p=", round(pvalue,3), sep=""),side = 3, line =0, adj = 0.25,
	col = "black",cex=0.8)}
	if(pvalue <0.001) {mtext("p<0.001", side = 3, line= 0, adj = 0.25, col="black",cex=0.8)}
	mtext(paste("slope= ", round(slope1a,3)," se ",
	round(summary(mod1a)[[4]][[7]],3),sep=""),side = 3,line=0,adj=0.75, col="black",cex=0.8)


abline("h"=GM,lty=2)
abline("v"=GM.mod1a,lty=2,col=rainbow(zN.u))
text(xmin,GM,GM,cex=1,pos=3)
text(GM.mod1a,ymin,GM.mod1a,pos=c(3,4),cex=1)
#..........................................................................................
#	Fig B
#..........................................................................................
title<-"Regression of TP on EQR"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*1/slope2a+y0[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}
	mtext(paste("R2=", round(summary(mod2a)$r.sq,3),sep=""), side = 3, line=0,adj=0, col="black",cex=0.8)
	pvalue<-pf(summary(mod2a)[[10]][[1]],summary(mod2a)[[10]][[2]],summary(mod2a)[[10]][[3]],lower.tail=F)
	if(pvalue >=0.001) {mtext(paste("p=", pvalue, sep=""),side = 3, line =0, adj = 0.25, col = "black",cex=0.8)}
	if(pvalue <0.001) {mtext("p<0.001", side = 3, line= 0, adj = 0.25, col="black",cex=0.8)}
	mtext(paste("slope= ", round(1/slope2a,3)," se ", 
		round(summary(mod2a)[[4]][[7]],3),sep=""), side = 3, line=0,adj=0.75, col="black",cex=0.8)


abline("h"=GM,lty=2)
abline("v"=GM.mod2a,lty=2,col=rainbow(zN.u))
text(xmin,GM,GM,cex=1,pos=3)
text(GM.mod2a,ymin,GM.mod2a,pos=c(3,4),cex=1)

#..........................................................................................
#	Fig C
#..........................................................................................
title<-"Geometric mean (SMA) regression EQR v TP"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb, xlab=xlb,las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*slope3a+int3a[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}
mtext(paste("slope= ", round(slope3a,3),sep=""), side = 3, line=0,adj=0.5, col="black",cex=0.8)
abline("h"=GM,lty=2)
abline("v"=GM.mod3a,lty=2,col=rainbow(zN.u))
text(xmin,GM,GM,cex=1,pos=3)
text(GM.mod3a,ymin,GM.mod3a,pos=c(3,4),cex=1)

#..........................................................................................
#	Fig D
#..........................................................................................
title<-""

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*slope1a+int1a[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}
legend("bottomleft",levels(z.u),col=rainbow(zN.u),pch=19,cex=1)

abline("h"=HG,lty=2)
abline("v"=HG.mod1a,lty=2,col=rainbow(zN.u))
text(xmin,HG,HG,cex=1,pos=3)
text(HG.mod1a,ymin,HG.mod1a,pos=c(3,4),cex=1)
#..........................................................................................
#	Fig E
#..........................................................................................
title<-""

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*1/slope2a+y0[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}
abline("h"=HG,lty=2)
abline("v"=HG.mod2a,lty=2,col=rainbow(zN.u))
text(xmin,HG,HG,cex=1,pos=3)
text(HG.mod2a,ymin,HG.mod2a,pos=c(3,4),cex=1)

#..........................................................................................
#	Fig F
#..........................................................................................
title<-""

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb, xlab=xlb,las=1)
points(y.u~x2.u,pch=19,col=rainbow(zN.u)[z.u])
new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
for (i in 1:zN.u){
	lines(10^(new.x$x),new.x$x*slope3a+int3a[i],col=rainbow(zN.u)[i])
	points(10^(mean(x.u[z.u==levels(z.u)[i]])),
		mean(y.u[z.u==levels(z.u)[i]]),pch=3,cex=3,col=rainbow(zN.u)[i])
	}

abline("h"=HG,lty=2)
abline("v"=HG.mod3a,lty=2,col=rainbow(zN.u))
text(xmin,HG,HG,cex=1,pos=3)
text(HG.mod3a,ymin,HG.mod3a,pos=c(3,4),cex=1)

###########################################################################################
#Output boundary values,into a series of data.frames called
#Out1
#Out2
#to
#OutN, where N is number of groups in data set (zN.u)

for(i in 1:zN.u){
		out<-data.frame(
		Type=c(levels(z.u)[i],"",""),
		R2=c(round(summary(mod1a)[[9]],3),"",round(summary(mod2a)[[9]],3)),
		N=c(length(x.u),"",""),
		slope=c(round(slope1a,3),round(slope3a,3),round(1/slope2a,3)),
		int=c(round(int1a[i],3),round(int3a[i],3),round(y0[i],3)),
		GM =c(GM.mod1a[i],GM.mod3a[i],GM.mod2a[i]),
		GML=c(GML.mod1a[i],GML.mod3a[i],GML.mod2a[i]),
		GMU=c(GMU.mod1a[i],GMU.mod3a[i],GMU.mod2a[i]),
		HG =c(HG.mod1a[i],HG.mod3a[i],HG.mod2a[i]),
		HGL=c(HGL.mod1a[i],HGL.mod3a[i],HGL.mod2a[i]),
		HGU=c(HGU.mod1a[i],HGU.mod3a[i],HGU.mod2a[i]),
		row.names=c("OLS EQR v TP","geometric mean regression","OLS TP v EQR"))
		assign(paste("Out",i,sep=""),out)
	}
# The following lines show output on console and then place in clipboard
print(Out1)
write.table(Out1,"clipboard",sep="\t",row.names=F,col.names=F)
#Now paste to Excel

print(Out2)
write.table(Out2,"clipboard",sep="\t",row.names=F,col.names=F)
#Now paste to Excel

print(Out3)
write.table(Out3,"clipboard",sep="\t",row.names=F,col.names=F)
#Now paste to Excel

print(Out4)
write.table(Out4,"clipboard",sep="\t",row.names=F,col.names=F)
#Now paste to Excel

#Add more rows if needed

