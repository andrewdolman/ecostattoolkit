#Description: Example script to fit linear models to data of a single category
#Geoff Phillips
#Date:09 Aug 2016
#File name for script: TKit_fit_lin_mod1.R

#Load required libraries
library(lmodel2)		# Type 2 regression

rm(list= ls())  		# Clear data

#################################################################################
#		Step 1 get the data		 	 					  #
#################################################################################

#Enter file name and path
#FName<-"DataTemplate_Example1.csv"
FName <- system.file("extdata", "DataTemplate_Example1.csv", package = "ecostattoolkit")

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

#------------------------------------------------------------------------------------------
#             Assign data used for models to x.u and y.u
#------------------------------------------------------------------------------------------
nut.min<-1	#min nutrient value used for model
nut.max<-138	#max nutrient value used for model 

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

length(x.u)
length(y.u)
#..........................................................................................
#	Model 1 OLS Y on X
#..........................................................................................
summary(mod1<-lm(y.u~x.u))		#Fit model
(int1<-summary(mod1)[[4]][[1]])
(slope1<-summary(mod1)[[4]][[2]])

# 	predict values of y from model 1 and calculate upper and lower 25th 75th 
#	quantiles of residuals
pred.y1<-data.frame(predict(mod1),resid(mod1))#Calculate predicted values and residuals
(upresid.1<-quantile(resid(mod1),0.75,na.rm=T))	#All residuals
(lwresid.1<-quantile(resid(mod1),0.25,na.rm=T))
#

(GM.mod1<-round(10^((GM-int1)/slope1),0))		#Predicted GM boundary value 
(GMU.mod1<-round(10^((GM-int1-upresid.1)/slope1),0))
(GML.mod1<-round(10^((GM-int1-lwresid.1)/slope1),0))

(HG.mod1<-round(10^((HG-int1)/slope1),0))			#Predicted HG boundary value 
(HGU.mod1<-round(10^((HG-int1-upresid.1)/slope1),0))
(HGL.mod1<-round(10^((HG-int1-lwresid.1)/slope1),0))

#..........................................................................................
#	Model 2 OLS X on Y
#..........................................................................................
summary(mod2<-lm(x.u~y.u))		#Fit model
(int2<-summary(mod2)[[4]][[1]])
(slope2<-summary(mod2)[[4]][[2]])			

# 	Calculate value of Y when X = 0 to allow average of intercepts of models and predict 
(y0<-(0-int2)/slope2)	
# 	or use the following form
(y0<-mean(y.u)- (mean(x.u)/slope2))


# Set up a data frame with observed x, factor z and add appropriate intercept value y0a for category
# calculate predicted y values (EQR) and calculate residuals
	
(pred.y2<-data.frame(x.u,y.u,y0))			
(pred.y2<-within(pred.y2,pred.y2<-x.u*1/slope2+y0))	#calc predicted values pred.y2
(pred.y2<-within(pred.y2,resid.2<-y.u-pred.y2))		#calc residuals resid.2

# 	predict values of y from model 1 and calculate upper and lower 25th 75th 
#	quantiles of residuals

upresid.2<-quantile(pred.y2$resid.2,0.75,na.rm=T)
lwresid.2<-quantile(pred.y2$resid.2,0.25,na.rm=T)



#
(GM.mod2<-round(10^((GM-y0)*slope2),0))			#Predicted GM boundary value 
(GMU.mod2<-round(10^((GM-y0-upresid.2)*slope2),0))
(GML.mod2<-round(10^((GM-y0-lwresid.2)*slope2),0))

#
(HG.mod2<-round(10^((HG-y0)*slope2),0))			#Predicted HG boundary value 
(HGU.mod2<-round(10^((HG-y0-upresid.2)*slope2),0))
(HGL.mod2<-round(10^((HG-y0-lwresid.2)*slope2),0))
#..........................................................................................
#	Model 4 Ranged Major Axis regression
#..........................................................................................
lmod<-lmodel2(y.u~x.u,range.y="relative",range.x="interval",nperm=99)
lmod						#Summary of all  models including conventional OLS
(RMA.int<-lmod$regression.results[[2]][[4]])
(RMA.slope<-lmod$regression.results[[3]][[4]])

# predict values of y from model and calculate upper and lower 25th 75th quantiles of residuals
pred.yRMA<-x.u*RMA.slope+RMA.int 
(upresid.RMA<-quantile(resid.RMA<-y.u-pred.yRMA,0.75,na.rm=T))
(lwresid.RMA<-quantile(resid.RMA<-y.u-pred.yRMA,0.25,na.rm=T))

(GM.mod4<-round(10^((GM-RMA.int)/RMA.slope),0))			#Predicted GM boundary value 
(GML.mod4<-round(10^((GM-RMA.int-lwresid.RMA)/RMA.slope),0))
(GMU.mod4<-round(10^((GM-RMA.int-upresid.RMA)/RMA.slope),0))

(HG.mod4<-round(10^((HG-RMA.int)/RMA.slope),0))			#Predicted HG boundary value
(HGL.mod4<-round(10^((HG-RMA.int-lwresid.RMA)/RMA.slope),0))
(HGU.mod4<-round(10^((HG-RMA.int-upresid.RMA)/RMA.slope),0))

#===========================================================================================
#	Step 3 
#	Plot models with upper and lower quantiles of residuals
#===========================================================================================
#	set scales for plotting and axis labels
xmin<- 5
xmax<- 500
ymin<- 0
ymax<- 1.5
xlb<-expression(paste("total phosphorus ","(?g ", L^-1,")"))
ylb<-"EQR" # Y axis label

#win.graph(width=14,height=4.5)					# new graphic window
#par(mfrow=c(1,3)) 						# delete this line if separate plots are required
par(mar=c(4,5,3.5,1))
#..........................................................................................
#	Fig A
#..........................................................................................
title<-"Regression of EQR on TP"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19)
points((10^mean(x.u)),mean(y.u),pch=3,cex=3) # plot mean of x and y on linear scale
mtext(paste("R2=", round(summary(mod1)$r.sq,3),sep=""), side = 3, line=0,adj=0, 
	col="black",cex=0.8)
pvalue<-pf(summary(mod1)[[10]][[1]],summary(mod1)[[10]][[2]],summary(mod1)[[10]][[3]],lower.tail=F)
if(pvalue >=0.001) {mtext(paste("p=", round(pvalue,3), sep=""),side = 3, line =0, adj = 0.25,
	col = "black",cex=0.8)}
if(pvalue <0.001) {mtext("p<0.001", side = 3, line= 0, adj = 0.25, col="black",cex=0.8)}
mtext(paste("slope= ", round(summary(mod1)[[4]][[2]],3)," se ",
	round(summary(mod1)[[4]][[4]],3),sep=""),side = 3,line=0,adj=0.75, col="black",cex=0.8)

new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
lines(10^(new.x$x),new.x$x*slope1+int1)
lines(10^(new.x$x),new.x$x*slope1+int1+upresid.1,lty=2)
lines(10^(new.x$x),new.x$x*slope1+int1+lwresid.1,lty=2)
abline("h"=GM,lty=1)
abline("v"=GM.mod1,lty=1)
abline("v"=GMU.mod1,lty=3)
abline("v"=GML.mod1,lty=3)
text(xmin,GM+0.03,GM,pos=4)
text(GM.mod1,ymin,GM.mod1,pos=3)
text(GMU.mod1,ymin,GMU.mod1,pos=4)
text(GML.mod1,ymin,GML.mod1,pos=2)

#..........................................................................................
#	Fig B
#..........................................................................................
title<-"Regression of TP on EQR"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb,xlab=xlb,
	las=1)
points(y.u~x2.u,pch=19)
points((10^mean(x.u)),mean(y.u),pch=3,cex=3) # plot mean of x and y on linear scale
mtext(paste("R2=", round(summary(mod2)$r.sq,3),sep=""), side = 3, line=0,adj=0, col="black",cex=0.8)
pvalue<-pf(summary(mod2)[[10]][[1]],summary(mod2)[[10]][[2]],summary(mod2)[[10]][[3]],lower.tail=F)
if(pvalue >=0.001) {mtext(paste("p=", pvalue, sep=""),side = 3, line =0, adj = 0.25, col = "black",cex=0.8)}
if(pvalue <0.001) {mtext("p<0.001", side = 3, line= 0, adj = 0.25, col="black",cex=0.8)}
mtext(paste("slope= ", round(1/summary(mod2)[[4]][[2]],3)," se ", 
	round(summary(mod2)[[4]][[4]],3),sep=""), side = 3, line=0,adj=0.75, col="black",cex=0.8)

new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = length(x)))#new values of X on log scale
lines(10^(new.x$x),new.x$x*1/slope2+y0)
lines(10^(new.x$x),new.x$x*1/slope2+y0+upresid.2,lty=2)
lines(10^(new.x$x),new.x$x*1/slope2+y0+lwresid.2,lty=2)
abline("h"=GM,lty=1)
abline("v"=GM.mod2,lty=1)
abline("v"=GMU.mod2,lty=3)
abline("v"=GML.mod2,lty=3)
text(xmin,GM+0.03,GM,pos=4)
text(GM.mod2,ymin,GM.mod2,pos=3)
text(GMU.mod2,ymin,GMU.mod2,pos=4)
text(GML.mod2,ymin,GML.mod2,pos=2)


#..........................................................................................
#	Fig C
#..........................................................................................
title<-"Ranged major axis regression EQR v TP"

plot(y~x2,ylim=c(ymin,ymax),xlim=c(xmin,xmax),main=title,log="x",ylab=ylb, xlab=xlb,las=1)
points(y.u~x2.u,pch=19)
points(10^(mean(x.u)),mean(y.u),pch=3,cex=3)
mtext(paste("RMA slope=", round(RMA.slope,3),sep=""), side = 3, line=0,adj=0.6, col="black",cex=0.8)

new.x <- data.frame(x = seq(min(x, na.rm=T), max(x, na.rm=T), length = 500))#new values of X on log scale
lines(10^(new.x$x),new.x$x*RMA.slope+RMA.int)
lines(10^(new.x$x),new.x$x*RMA.slope+RMA.int+upresid.RMA,lty=2)
lines(10^(new.x$x),new.x$x*RMA.slope+RMA.int+lwresid.RMA,lty=2)
abline("h"=GM,lty=1)
abline("v"=GM.mod4,lty=1)
abline("v"=GMU.mod4,lty=3)
abline("v"=GML.mod4,lty=3)
text(xmin,GM+0.03,GM,pos=4)
text(GM.mod4,ymin,GM.mod4,pos=3)
text(GMU.mod4,ymin,GMU.mod4,pos=4)
text(GML.mod4,ymin,GML.mod4,pos=2)


###########################################################################################
#Output boundary values, 

(out1<-data.frame(
	R2=c(round(summary(mod1)[[9]],3),"",""),
	N=c(length(x.u),"",""),
	slope=c(slope1,RMA.slope,1/slope2),
	int=c(int1,RMA.int,y0),
	GM =c(GM.mod1,GM.mod4,GM.mod2),
	GML=c(GML.mod1,GML.mod4,GML.mod2),
	GMU=c(GMU.mod1,GMU.mod4,GMU.mod2),
	HG =c(HG.mod1,HG.mod4,HG.mod2),
	HGL=c(HGL.mod1,HGL.mod4,HGL.mod2),
	HGU=c(HGU.mod1,HGU.mod4,HGU.mod2)))
write.table(out1,"clipboard",sep="\t",row.names=F,col.names=F)


