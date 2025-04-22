library(ggplot2)
library(MASS)
library(plyr)
library(data.table)

## Make sure to set working directory so that you can access log files output by BEAST runs
## Names of log files before the final underscore will differ depending on the data set and model
## Please consult XML files for names and change as appropriate 

## Logged samples for siteAssignInd 
## (for each site, there is nonnegative integer value corresponding to evolutionary category)
siteAssignInd = as.matrix(fread("rsva_dp_siteassignind.log"))

## Logged samples for substitution rates by (occupied and unoccupied) evolutionary category
rates = as.matrix(fread("rsva_dp_rates.log"))

## Logged samples for GTR rate matrix parameters by (occupied and unoccupied) evolutionary category
gtrrates = as.matrix(fread("rsva_dp_gtrrates.log"))

## Logged samples for equilibrium frequencies by (occupied and unoccupied) evolutionary category
freq = as.matrix(fread("rsva_dp_freq.log"))

## Total number of alignment sites
numsites = dim(siteAssignInd)[2]-1

## Total number of logged samples (without initial values and before discarding burnin)
initialnumsamples = dim(siteAssignInd)[1]-1

## Number of potential categories (occupied and unoccupied)
numcat = dim(rates)[2]-1

## Store values for different nucleotide base equilibrium frequencies
freq_1 <- freq[,2+0:(numcat-1)*4]
freq_2 <- freq[,3+0:(numcat-1)*4]
freq_3 <- freq[,4+0:(numcat-1)*4]
freq_4 <- freq[,5+0:(numcat-1)*4]

## Store values for different GTR rate matrix parameters
gtrrates_1 <- gtrrates[,2+0:(numcat-1)*6]
gtrrates_2 <- gtrrates[,3+0:(numcat-1)*6]
gtrrates_3 <- gtrrates[,4+0:(numcat-1)*6]
gtrrates_4 <- gtrrates[,5+0:(numcat-1)*6]
gtrrates_5 <- gtrrates[,6+0:(numcat-1)*6]
gtrrates_6 <- gtrrates[,7+0:(numcat-1)*6]

## Number of burnin iterations in millions of states. 
## For example, "burnInMil <- 10" corresponds to a burnin of 10 million states 
## MAKE SURE TO SET APPROPRIATE BURNIN
## DO NOT BLINDLY USE DEFAULT VALUE
burninInMil <- 10

## Converts burnInMil (burnin in terms of millions of states) to number of logged samples
burnin = (burninInMil*1000000)/(10000)

## Number of samples remaining after discarding burnin
numsamples = initialnumsamples - burnin

## Matrix for site-specific substitution rates
## Row corresponds to sample and column to site
siteRateVals <- matrix(0,nrow=numsamples,ncol=numsites)

## Matrices for site-specific equilibrium frequencies
## Rows corresponds to sample and columns to site
siteFreqVals_1 <- matrix(0,nrow=numsamples,ncol=numsites)
siteFreqVals_2 <- matrix(0,nrow=numsamples,ncol=numsites)
siteFreqVals_3 <- matrix(0,nrow=numsamples,ncol=numsites)
siteFreqVals_4 <- matrix(0,nrow=numsamples,ncol=numsites)

## Matrices for site-specific GTR rate matrix parameters
## Rows corresponds to sample and columns to site
siteGTRRateVals_1 <- matrix(0,nrow=numsamples,ncol=numsites)
siteGTRRateVals_2 <- matrix(0,nrow=numsamples,ncol=numsites)
siteGTRRateVals_3 <- matrix(0,nrow=numsamples,ncol=numsites)
siteGTRRateVals_4 <- matrix(0,nrow=numsamples,ncol=numsites)
siteGTRRateVals_5 <- matrix(0,nrow=numsamples,ncol=numsites)
siteGTRRateVals_6 <- matrix(0,nrow=numsamples,ncol=numsites)

## Matrix for number of sites in each potential evolutionary category for each sample
## Row corresponds to sample and column to evolutionary category
catCounts = matrix(0,nrow=numsamples,ncol=numcat)

## Number of active/occupied categories for each sample
numActiveCats <- seq(0,0,length=numsamples)

## Extract site-specific parameter values for each sample
for(i in 1:numsamples){	
	siteRateVals[i,1:numsites] = rates[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+2]
	siteGTRRateVals_1[i,1:numsites] = gtrrates_1[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteGTRRateVals_2[i,1:numsites] = gtrrates_2[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteGTRRateVals_3[i,1:numsites] = gtrrates_3[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteGTRRateVals_4[i,1:numsites] = gtrrates_4[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteGTRRateVals_5[i,1:numsites] = gtrrates_5[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteGTRRateVals_6[i,1:numsites] = gtrrates_6[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteFreqVals_1[i,1:numsites] = freq_1[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteFreqVals_2[i,1:numsites] = freq_2[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteFreqVals_3[i,1:numsites] = freq_3[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	siteFreqVals_4[i,1:numsites] = freq_4[burnin+i,siteAssignInd[burnin+i,2:(numsites+1)]+1]
	catCounts[i,siteAssignInd[burnin+i,2:(numsites+1)]+1] = catCounts[i,siteAssignInd[burnin+i,2:(numsites+1)]+1]+1	
	numActiveCats[i] = numActiveCats[i] + sum((catCounts[i,1:numcat] > 0) == TRUE)
}

## Compute summary statistics
lowerRates <- apply(siteRateVals[,,drop=F],2,quantile,probs=0.025)
upperRates <- apply(siteRateVals[,,drop=F],2,quantile,probs=0.975)
meanRates <- apply(siteRateVals[,,drop=F],2,mean)
medianRates <- apply(siteRateVals[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates1 <- apply(siteGTRRateVals_1[,,drop=F],2,quantile,probs=0.025)
upperGTRRates1 <- apply(siteGTRRateVals_1[,,drop=F],2,quantile,probs=0.975)
meanGTRRates1 <- apply(siteGTRRateVals_1[,,drop=F],2,mean)
medianGTRRates1 <- apply(siteGTRRateVals_1[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates2 <- apply(siteGTRRateVals_2[,,drop=F],2,quantile,probs=0.025)
upperGTRRates2 <- apply(siteGTRRateVals_2[,,drop=F],2,quantile,probs=0.975)
meanGTRRates2 <- apply(siteGTRRateVals_2[,,drop=F],2,mean)
medianGTRRates2 <- apply(siteGTRRateVals_2[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates3 <- apply(siteGTRRateVals_3[,,drop=F],2,quantile,probs=0.025)
upperGTRRates3 <- apply(siteGTRRateVals_3[,,drop=F],2,quantile,probs=0.975)
meanGTRRates3 <- apply(siteGTRRateVals_3[,,drop=F],2,mean)
medianGTRRates3 <- apply(siteGTRRateVals_3[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates4 <- apply(siteGTRRateVals_4[,,drop=F],2,quantile,probs=0.025)
upperGTRRates4 <- apply(siteGTRRateVals_4[,,drop=F],2,quantile,probs=0.975)
meanGTRRates4 <- apply(siteGTRRateVals_4[,,drop=F],2,mean)
medianGTRRates4 <- apply(siteGTRRateVals_4[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates5 <- apply(siteGTRRateVals_5[,,drop=F],2,quantile,probs=0.025)
upperGTRRates5 <- apply(siteGTRRateVals_5[,,drop=F],2,quantile,probs=0.975)
meanGTRRates5 <- apply(siteGTRRateVals_5[,,drop=F],2,mean)
medianGTRRates5 <- apply(siteGTRRateVals_5[,,drop=F],2,quantile,probs=0.5)

lowerGTRRates6 <- apply(siteGTRRateVals_6[,,drop=F],2,quantile,probs=0.025)
upperGTRRates6 <- apply(siteGTRRateVals_6[,,drop=F],2,quantile,probs=0.975)
meanGTRRates6 <- apply(siteGTRRateVals_6[,,drop=F],2,mean)
medianGTRRates6 <- apply(siteGTRRateVals_6[,,drop=F],2,quantile,probs=0.5)

lowerFreq1 <- apply(siteFreqVals_1[,,drop=F],2,quantile,probs=0.025)
upperFreq1 <- apply(siteFreqVals_1[,,drop=F],2,quantile,probs=0.975)
meanFreq1 <- apply(siteFreqVals_1[,,drop=F],2,mean)
medianFreq1 <- apply(siteFreqVals_1[,,drop=F],2,quantile,probs=0.5)
	
lowerFreq2 <- apply(siteFreqVals_2[,,drop=F],2,quantile,probs=0.025)
upperFreq2 <- apply(siteFreqVals_2[,,drop=F],2,quantile,probs=0.975)
meanFreq2 <- apply(siteFreqVals_2[,,drop=F],2,mean)
medianFreq2 <- apply(siteFreqVals_2[,,drop=F],2,quantile,probs=0.5)
		
lowerFreq3 <- apply(siteFreqVals_3[,,drop=F],2,quantile,probs=0.025)
upperFreq3 <- apply(siteFreqVals_3[,,drop=F],2,quantile,probs=0.975)
meanFreq3 <- apply(siteFreqVals_3[,,drop=F],2,mean)
medianFreq3 <- apply(siteFreqVals_3[,,drop=F],2,quantile,probs=0.5)
		
lowerFreq4 <- apply(siteFreqVals_4[,,drop=F],2,quantile,probs=0.025)
upperFreq4 <- apply(siteFreqVals_4[,,drop=F],2,quantile,probs=0.975)
meanFreq4 <- apply(siteFreqVals_4[,,drop=F],2,mean)
medianFreq4 <- apply(siteFreqVals_4[,,drop=F],2,quantile,probs=0.5)


## Create data frame with summary statistics 

site <- seq(1,numsites,length=numsites)

df <- data.frame(site=site,lowerRates=lowerRates,medianRates=medianRates,upperRates=upperRates,lowerFreq1=lowerFreq1,medianFreq1=medianFreq1,upperFreq1=upperFreq1,lowerFreq2=lowerFreq2,medianFreq2=medianFreq2,upperFreq2=upperFreq2,lowerFreq3=lowerFreq3,medianFreq3=medianFreq3,upperFreq3=upperFreq3,lowerFreq4=lowerFreq4,medianFreq4=medianFreq4,upperFreq4=upperFreq4,lowerGTRRates1=lowerGTRRates1,medianGTRRates1=medianGTRRates1,upperGTRRates1=upperGTRRates1,lowerGTRRates2=lowerGTRRates2,medianGTRRates2=medianGTRRates2,upperGTRRates2=upperGTRRates2,lowerGTRRates3=lowerGTRRates3,medianGTRRates3=medianGTRRates3,upperGTRRates3=upperGTRRates3,lowerGTRRates4=lowerGTRRates4,medianGTRRates4=medianGTRRates4,upperGTRRates4=upperGTRRates4,lowerGTRRates5=lowerGTRRates5,medianGTRRates5=medianGTRRates5,upperGTRRates5=upperGTRRates5,lowerGTRRates6=lowerGTRRates6,medianGTRRates6=medianGTRRates6,upperGTRRates6=upperGTRRates6)

df$groups <- 1

# K-Means clustering
km <- kmeans(df[2:(dim(df)[2]-1)],median(numActiveCats))

# Create new categories that are ordered by median substitution rate
newcats <- seq(0,0,length=numsites)
for(i in 1:numsites){
	newcats[i]<- which(order(km$centers[,'medianRates'],decreasing=FALSE) == km$cluster[i])
}

# Create plot
p <- ggplot() + geom_rect(data = df, aes(xmin = as.numeric(site), xmax = as.numeric(site)+1,ymin = (newcats-1)/median(numActiveCats), ymax = (newcats)/median(numActiveCats),fill = as.factor(groups)), show.legend=FALSE)+labs(x = "Site", y = "Category") + scale_color_manual(values=c("black", "black")) +scale_fill_manual(values=c("black", "black")) +scale_y_continuous(expand=c(0.002,0.002),limits=c(0,1), breaks=(seq(1/median(numActiveCats),1,length=median(numActiveCats))- 1/(2*median(numActiveCats))), labels=seq(1,median(numActiveCats),length=median(numActiveCats))) + scale_x_continuous(expand=c(0.006,0.002)) + theme(axis.title=element_text(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + theme_classic() + ggtitle("GTR + DP") + theme(plot.title = element_text(hjust = 0.5)) + theme(text=element_text(size=25)) + theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

p

dev.copy2pdf(file="clusterpattern_rsvagtrdp.pdf",height=6,width=36)