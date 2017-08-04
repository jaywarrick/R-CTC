source('~/GitHub/R-General/.Rprofile')
source('~/GitHub/R-Cytoprofiling/PreprocessingHelpers.R')
library(data.table)

# Edit these lines of code for each patient/database
patient <- 'MDV1-GR'
Threshold.CTC <- 1 # This is in log space so the actual ratio threshold is exp(Threshold) = 2.7^Threshold
Dump.max <- 1000
Cyt.min <- 0
Area.min <- 0 # Cells were typically about 2500
Area.max <- 5000 # Cells were typically about 2500. This will at least exclude REALLY big things.

# Then click the 'Source' button to run everything and look
# within the 'J:/Rory/R Output Folder' for all results.
# However, there are some other plots made within the RStudio
# window that might be interesting to look at as well (histograms etc.).

# Distance from y=x line
dist <- function(x, y)
{
  a <- -1
  b <- 1
  c <- 0
  temp <- (a*x + b*y + c)/(sqrt(a^2+b^2))
  return(temp)
}

x <- fread(paste0('J:/Rory/', patient, '/Dataset Name/Cell_x0_y0/File-Converted Coloc Table/x0_y0.csv'))
replaceStringInColNames(x, '$','.')
removeColNamesContaining(x, 'net.imagej.ops.Ops.Stats.SpearmansRankCorrelationCoefficient')
replaceStringInColNames(x, 'net.imagej.ops.Ops.Stats.PearsonsCorrelationCoefficient.', '')

y <- fread(paste0('J:/Rory/', patient, '/Dataset Name/Cell_x0_y0/File-Converted Feature Table/x0_y0.csv'))
replaceStringInColNames(y, '$','.')
removeColNamesContaining(y, '_5')
replaceStringInColNames(y, 'net.imagej.ops.Ops.Stats.', '')
# View(y)

Ids <- c("Id","Label","Experiment","Array.X","Array.Y","Valid")
z <- merge(x, y, by=Ids)

# Calculate the area of the cells
z[, Area_WholeCell:=Sum_WholeCell_1/Mean_WholeCell_1]

# Throw out variables with infinities etc.
z <- lapply.data.table(z, FUN=function(n){if(sum(!is.finite(n)) > 0){return(NULL)}else{return(n)}}, in.place=F)

# Calculate Nuc.Cyt.Ratio
z[, Nuc.Cyt.Ratio:=Similarity_WholeCell_1_2/Similarity_WholeCell_2_3]
# Calculate Nuc.Cyt.Score
z[, Nuc.Cyt.Score:=dist(Similarity_WholeCell_2_3, Similarity_WholeCell_1_2)]

z <- z[Area_WholeCell > Area.min & Area_WholeCell < Area.max]

# This is your table with both colocalization and intensity information.
# View(z)

# Determine which are CTCs
z[, CTCScore:=log(Mean_WholeCell_3/Mean_WholeCell_4)]
clusterResults <- assignToClusters(z$CTCScore, nClusters=2)
z[, Cluster:=clusterResults$data$Cluster.Clean]
z[, CTC:=CTCScore > Threshold.CTC & Mean_WholeCell_3 > Cyt.min & Mean_WholeCell_4 < Dump.max]
plotClusters(z$CTCScore, z$Cluster)

# Plot thresholds for determining CTCs
# This will save the plot to a file
pdf(paste0('J:/Rory/R Output Folder/', patient, '_CTCs.pdf'), width=5, height=4)
plot(z$Mean_WholeCell_4, z$Mean_WholeCell_3, pch=20, xlab='Mean Dump Intensity', ylab='Mean Cytokeratin Intensity', col='black', log='xy')
threshold.x <- seq(min(z$Mean_WholeCell_4), max(z$Mean_WholeCell_4), length.out=100)
threshold.y <- threshold.x*exp(Threshold.CTC)
lines(threshold.x, threshold.y, lty=2, col='red')
abline(v=Dump.max, lty=2, col='red')
abline(h=Cyt.min, lty=2, col='red')
dev.off()

# Store which cells were determined as CTCs
CTCList <- z$Id[z$CTC]
CTCs <- z[CTC==T]

# Summarize percent of cells with preferential colocalization with Nuc rather than Cytokeratin
CTC.summary <- CTCs[, list(Nuc.Cyt.Ratio.Mean=mean(Nuc.Cyt.Ratio), Nuc.Cyt.Ratio.SE=sd(Nuc.Cyt.Ratio)/sqrt(.N), Nuc.Cyt.Score.Mean=mean(Nuc.Cyt.Score), Nuc.Cyt.Score.SE=sd(Nuc.Cyt.Score)/sqrt(.N), Percent.Nuc.AR.Pos=sum(Nuc.Cyt.Ratio > 1)/(.N))]

# See a plot of AR colocalization with cytokeratin vs nuclear
pdf(paste0('J:/Rory/R Output Folder/', patient, '_Localization.pdf'), width=5, height=4)
plot(z$Similarity_WholeCell_2_3, z$Similarity_WholeCell_1_2, ylab='AR Colocalization with Nuc', xlab='AR Colocalization with Cytokeratin', pch=20, col='blue')
points(CTCs$Similarity_WholeCell_2_3, CTCs$Similarity_WholeCell_1_2, pch=20, col='green')
# plot the line where Nuc.Cyt.Score is zero (i.e., Nuc.Cyt.Ratio = 1)
abline(a=0, b=1, lty=2)
#plot the mean Nuc.Cyt.Score for just the CTCs
abline(a=mean(CTCs$Nuc.Cyt.Score)/cosd(45), b=1, lty=2, col='red')
dev.off()

# Make some histograms of the individual measures
# Whole Table
hist(z$Similarity_WholeCell_1_2, breaks=40)
hist(z$Similarity_WholeCell_2_3, breaks=40)
hist(z$Nuc.Cyt.Ratio, breaks=40)
abline(v=1, lty=2, col='red')
hist(z$Nuc.Cyt.Score, breaks=40)
abline(v=1, lty=2, col='red')
# CTCs Only
hist(CTCs$Nuc.Cyt.Ratio, breaks=40)
abline(v=1, lty=2, col='red')
hist(CTCs$Nuc.Cyt.Score, breaks=40)
abline(v=1, lty=2, col='red')

# Add a column in z indicating which are CTCs
z$CTC <- z$Id %in% CTCList

# Save z to a file that can be opened in Excel
fwrite(z, paste0('J:/Rory/R Output Folder/', patient, '.csv'))
# Save the CTC only table to a file that can be opened in Excel
fwrite(CTCs, paste0('J:/Rory/R Output Folder/', patient, '_CTCs.csv'))
# Save the CTC summary table to a file that can be opened in Excel
fwrite(CTC.summary, paste0('J:/Rory/R Output Folder/', patient, '_CTCSummary.csv'))

# This will show the plot in RStudio
plot(z$Similarity_WholeCell_2_3, z$Similarity_WholeCell_1_2, ylab='AR Colocalization with Nuc', xlab='AR Colocalization with Cytokeratin', pch=20, col='blue')
points(CTCs$Similarity_WholeCell_2_3, CTCs$Similarity_WholeCell_1_2, pch=20, col='green')
# plot the line where Nuc.Cyt.Score is zero (i.e., Nuc.Cyt.Ratio = 1)
abline(a=0, b=1, lty=2)
#plot the mean Nuc.Cyt.Score for just the CTCs
abline(a=mean(CTCs$Nuc.Cyt.Score)/cosd(45), b=1, lty=2, col='red')

# This will show the plot in RStudio
plot(z$Mean_WholeCell_4, z$Mean_WholeCell_3, pch=20, xlab='Mean Dump Intensity', ylab='Mean Cytokeratin Intensity', col='black', log='xy')
threshold.x <- seq(min(z$Mean_WholeCell_4), max(z$Mean_WholeCell_4), length.out=100)
threshold.y <- threshold.x*exp(Threshold.CTC)
lines(threshold.x, threshold.y, lty=2, col='red')
abline(v=Dump.max, lty=2, col='red')
abline(h=Cyt.min, lty=2, col='red')

print(CTC.summary)
