#Description: Example script to check data for linearity
#Geoff Phillips
#Date:09 August 2016
#File name for script: TKit_check_linearity.R

#Load required libraries
library(mgcv)				# gam modelling
library(segmented)			# segmented regression

rm(list= ls())  		# Clear data

#################################################################################
#		Step 1 get the data		 	 					  #
#################################################################################

#Enter file name and path
FName<-"DataTemplate_Example1.csv"

#Get data and check
data <- read.csv(file = FName, header = TRUE)
dim(data)
summary(data)

#Identify the records with both EQR and TP data to produce a data set without missing values
cc<-complete.cases(data$EQR,data$P)
data.cc<-data[cc,]
dim(data.cc)

#All data including outliers (for information)
x<-log10(data.cc$P)	#log10 transform of independent variable
x2<-10^x			#back transformed value for plotting on linear scale
y<-data.cc$EQR
z<-as.factor(data.cc$Group)# needs to be a factor to act as categorical variable
zN<-length(levels(z))	#number of categories in data

#Data used for modelling, excluding outliers
x.u<-log10(data.cc$P[data.cc$Exclude_P==0])# selects records where Exclude_P = 0 (note R requires ==)
x2.u<-10^x.u
y.u<-data.cc$EQR[data.cc$Exclude_P==0]
z.u<-as.factor(data.cc$Group[data.cc$Exclude_P==0])
z.uN<-length(levels(z.u))	#number of categories in data

#Plot the data

#Set up axis labels and ranges, edit these values for your data
xmin<-1
xmax<-1000
ymin<-0
ymax<-1.5
xlb<-expression(paste("total phosphorus ","(µg ", L^-1,")"))
ylb<-"EQR"

#Plot the data
plot(y~x2,log="x",ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab=xlb,ylab=ylb,cex=1)
points(y.u~x2.u,col=rainbow(z.uN)[z.u],pch=19)
legend("topright",levels(z.u),col=rainbow(z.uN),pch=19,cex=0.8)

#Fit a gam model and add a line to the plot
k<-9 #number of knots used in smoother, can be decreased to produce a smoother curve
fit.gam<- gam(y.u~s(x.u,k=k))
print(summary(fit.gam))
new.x.u <- data.frame(x.u = seq(min(x.u, na.rm=T), max(x.u, na.rm=T), length = length(x.u)))
pred.gam <- predict(fit.gam, newdata=new.x.u, se.fit=TRUE)
ylower <- pred.gam$fit - pred.gam$se.fit
yupper <- pred.gam$fit + pred.gam$se.fit
lines(10^(new.x.u$x.u),(pred.gam$fit), col="black")
mtext("GAM", side = 3, adj=0,line =0, col="black",cex=0.8)
mtext(paste("R2=", round(summary(fit.gam)$r.sq,3),sep=""), side = 3, line=0,adj=.5, col="black",cex=0.8)
pvalue <- summary(fit.gam)$s.pv
if(pvalue >=0.001) {mtext(paste("p=", round(pvalue,3), sep=""),side = 3, line =0, adj = 1, col = "black",cex=0.8)}
if(pvalue <0.001) {mtext("p<0.001", side = 3, line= 0, adj = 1, col="black",cex=0.8)}

#Fit a segmented regression 

#First run the following lines of script which set up 2 user functions

#.................................................................................................................
Myseg<-function(x2,y,EstBk)
{
	Eb<-log10(EstBk)
	x<-log10(x2)
	print(mod<-lm(y ~ x))
	print(o<-segmented(mod,seg.Z=~x,psi=list(x = c(Eb))))#--fit segmented model
	print(summary(o))
	print(bkpt<-10^(o$psi[1,2]))
	segments(10^(o$psi[1,2]-o$psi[1,3]),ymin,10^(o$psi[1,2]+o$psi[1,3]),ymin)
	points(10^(o$psi[1,2]),ymin,pch=3,cex=1.5)
	text(10^(o$psi[1,2]),ymin,round(10^(o$psi[1,2])),pos=3)
	print(10^(o$psi[1,2]-o$psi[1,3])) # lower estimate of break point
	print(10^(o$psi[1,2]+o$psi[1,3]))
	(int<-coef(o)[1])
	(slope<-coef(o)[2])
	(slope2<-coef(o)[2] + coef(o)[3])
	(int2<- int + slope * o$psi[1,2] - slope2 * o$psi[1,2] )
	abline(int,slope, col="blue")
	abline(int2,slope2,col="green")
}
#......................................................................................................................
Myseg2<-function(x2,y,EstBk,EstBk2)
{
	Eb<-log10(EstBk)
	Eb2<-log10(EstBk2)
	x<-log10(x2)
	print(mod<-lm(y ~ x))
	print(o<-segmented(mod,seg.Z=~x,psi=list(x = c(Eb,Eb2))))#--fit segmented model
	print(summary(o))
	#First break point
	print(bkpt<-10^(o$psi[1,2]))
	segments(10^(o$psi[1,2]-o$psi[1,3]),ymin,10^(o$psi[1,2]+o$psi[1,3]),ymin)
	points(10^(o$psi[1,2]),ymin,pch=3,cex=1.5)
	text(10^(o$psi[1,2]),ymin,round(10^(o$psi[1,2])),pos=3)
	print(10^(o$psi[1,2]-o$psi[1,3])) # lower estimate of break point
	print(10^(o$psi[1,2]+o$psi[1,3]))
	#Second break point
	print(bkpt2<-10^(o$psi[2,2]))
	segments(10^(o$psi[2,2]-o$psi[2,3]),ymin,10^(o$psi[2,2]+o$psi[2,3]),ymin)
	points(10^(o$psi[2,2]),ymin,pch=3,cex=1.5)
	text(10^(o$psi[2,2]),ymin,round(10^(o$psi[2,2])),pos=3)
	print(10^(o$psi[2,2]-o$psi[2,3])) # lower estimate of break point
	print(10^(o$psi[2,2]+o$psi[2,3]))
	#Calculate slopes and intercepts for each segment
	#Seg1
	(int<-coef(o)[1])
	(slope<-coef(o)[2])
	abline(int,slope, col="green")
	#Seg2
	(slope2<-coef(o)[2] + coef(o)[3])
	(int2<- int + slope * o$psi[1,2] - slope2 * o$psi[1,2] )
	abline(int2,slope2,col="blue")
	#Seg3
	(slope3<-coef(o)[3] + coef(o)[4])
	(int3<- int2 + slope2 * o$psi[2,2] - slope3 * o$psi[2,2] )
	abline(int3,slope3,col="red")
}
#.................................................................................................................

# For single break point use this section of script

EstBk<-75  # edit this line with estimated break point
Myseg(x2.u,y.u,EstBk)

# For 2 break points use this section of script

EstBk<-25  # edit this line with estimated lower break point
EstBk2<-100 # edit this line with estimated upper break point
Myseg2(x2.u,y.u,EstBk,EstBk2)



