rm(list=ls())
source('/Users/jaywarrick/Desktop/CTC/AnalyzeCTCs_HelperFunctions.R')
library(foreign)

locResults <- reorganizeTable(read.arff('/Users/jaywarrick/Desktop/CTC/JEX Results/Colocalization 3 - x0_y0.arff'))

temp <- calculateSummaryVariables(locResults=locResults, dapiColor=1, ckColor=2, cd45Color=3, arColor=4)
Thresh <- temp$thresholds
results <- temp$summary

CTCs <- subset(results, !CD45.Flag & AR.Flag)
nonCTCs <- subset(results, !CD45.Flag & !AR.Flag)

# Plot the density plots
CTCden <- density(CTCs$AR.NucPer, na.rm=T)
nonCTCden <- density(nonCTCs$AR.NucPer, na.rm=T)
plot(CTCden, col='red', main='')
lines(nonCTCden, col='black')

# Plot Nuclear versus Tot
plot(nonCTCs$AR.Nuc, nonCTCs$AR, pch=21, col=rgb(0,0,0,0.25), bg=rgb(0,0,0,0.25), xlim=c(1,max(nonCTCs$AR.Nuc, CTCs$AR.Nuc)), ylim=c(1,max(nonCTCs$AR, CTCs$AR)), xlab='AR Nuc Mean Intensity [au]', ylab='AR Mean Intensity [au]', log='xy')
points(CTCs$AR.Nuc, CTCs$AR, pch=21, col=rgb(1,0,0,0.7), bg=rgb(1,0,0,0.8))
fracLinesTot(0)
fracLinesTot(0.10)
fracLinesTot(0.20) # comsldlfjlskdjf
fracLinesTot(0.30)
fracLinesTot(0.40)
fracLinesTot(0.50)
fracLinesTot(0.60)
fracLinesTot(0.70)
fracLinesTot(0.80)
fracLinesTot(0.90)
fracLinesTot(1.00)
totLinesTot(10)
totLinesTot(100)
totLinesTot(1000)
totLinesTot(10000)
totLinesTot(100000)
totLinesTot(1000000)

# Plot Nuclear vs Cytoplasmic
plot(nonCTCs$AR.Nuc ,nonCTCs$AR.Cyt, pch=21, col=rgb(0,0,0,0.25), bg=rgb(0,0,0,0.25), xlim=c(1,max(nonCTCs$AR.Nuc, CTCs$AR.Nuc)), ylim=c(1,max(nonCTCs$AR.Cyt, CTCs$AR.Cyt)),xlab='AR Nuc Mean Intensity [au]', ylab='AR Cyt Mean Intensity [au]', log='xy')
points(CTCs$AR.Nuc ,CTCs$AR.Cyt, pch=21, col=rgb(1,0,0,0.7), bg=rgb(1,0,0,0.8))
fracLinesCyt(0)
fracLinesCyt(0.10)
fracLinesCyt(0.20)
fracLinesCyt(0.30)
fracLinesCyt(0.40)
fracLinesCyt(0.50)
fracLinesCyt(0.60)
fracLinesCyt(0.70)
fracLinesCyt(0.80)
fracLinesCyt(0.90)
fracLinesCyt(1.00)
totLinesCyt(10)
totLinesCyt(100)
totLinesCyt(1000)
totLinesCyt(10000)
totLinesCyt(100000)
totLinesCyt(1000000)

plot(results$AR ,results$AR.NucPer, pch=21, col=rgb(0,0,0,0.25), bg=rgb(0,0,0,0.25), xlim=c(1,max(results$AR, CTCs$AR)), ylim=c(0,max(results$AR.NucPer, na.rm=T)),xlab='AR Mean Intensity [au]', ylab='% Nuc AR', log='x')
points(CTCs$AR ,CTCs$AR.NucPer, pch=21, col=rgb(1,0,0,0.7), bg=rgb(1,0,0,0.8))
