#Description: Example script to check data to be used for establishing
#		nutrient boundary values
#Geoff Phillips
#Date:09 Aug 2016
#File name for script: TKit_check_data.R

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

#################################################################################
#		Step 2 plot data to identify potential outliers				  #
#################################################################################
win.graph(width=8,height=4)	#set up a new graphic window
par(mfrow=c(1,2))			#produce plots next to each other

#Box plots by biological class
boxplot(data$P ~ as.factor(data$BioClass),ylab="P conc")$out	# a plot and printed list outliers
boxplot(data$N ~ as.factor(data$BioClass),ylab="N conc")$out	

#Box plots by biological class and water body type
boxplot(data$P ~ as.factor(data$BioClass)+as.factor(data$Group),
	ylab="P conc",las=3)$out

#Scatter plots
par(mfrow=c(1,2))	
plot(data$EQR ~ data$P,log="x")#note the nutrient axis is set to a log scale
identify(data$P,data$EQR,data$Record)
plot(data$EQR ~ data$N,log="x")#note the nutrient axis is set to a log scale
identify(data$N,data$EQR,data$Record)

#################################################################################
# Now edit the data file marking the outliers that should not be used           #
#################################################################################
