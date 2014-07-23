reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE)
{
    library(plyr)
    idCols <- names(data)
    idCols <- idCols[-which(idCols %in% c('Measurement','Value'))]
    newData <- data.frame(stringsAsFactors=FALSE)
    measurements <- unique(data$Measurement)
    #     m <- measurements[1]
    #     browser()
    for(m in measurements)
    {
        if(is.na(baseName))
        {
            newColName <- m
            newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
        }else
        {
            newColName <- paste(baseName,'.',m, sep='')
            newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
        }
        
        temp <- subset(data, Measurement==m)
        temp2 <- temp[,c(idCols,'Value')]
        names(temp2)[names(temp2)=="Value"] <- newColName
        if(nrow(newData) == 0)
        {
            newData <- temp2
        }else
        {
            newData <- merge(newData, temp2, by=idCols)
        }
    }
    
    if(convertToNumeric)
    {
        for(n in idCols)
        {
            newData[,n] <- as.numeric(as.character(newData[,n]))
        }
    }
    
    return(newData)
}

frac <- function(x)
{
    return(length(which(x))/length(x))
}

fracLinesTot <- function(frac)
{
    col <- 'black'
    if(frac == 0.5)
    {
        col <- 'red'
    }
    #     ratio <- frac / (1-frac)
    x <- 0:100000
    x[1] <- 0.00001
    y <- (1/frac)*x
    lines(x,y, col=col)
}

fracLinesCyt <- function(frac)
{
    col <- 'black'
    if(frac == 0.5)
    {
        col <- 'red'
    }
    #     ratio <- frac / (1-frac)
    x <- 0:100000
    x[1] <- 0.00001
    y <- (1/frac - 1)*x
    lines(x,y, col=col)
    
}

totLinesCyt <- function(max)
{
    x <- (max/100)*0:100
    y <- max*sin((pi/2)*((100:0)/100))
    x[1] <- x[1] + 0.0001
    y[length(y)] <- y[length(y)] + 0.0001
    lines(x, y)
}

totLinesTot <- function(max)
{
    abline(h=max)
}

calculateSummaryVariables <- function(locResults, dapiColor, ckColor, cd45Color, arColor)
{
    locResults$Id <- as.numeric(as.character(locResults$Id))
    locResults <- locResults[with(locResults, order(Id)),]
    
    dapiTot <- paste('Sig_', dapiColor, ckColor, '_Tot', sep ='')
    ckTot <- paste('Sig_', ckColor, dapiColor, '_Tot', sep ='')
    cd45Tot <- paste('Sig_', cd45Color, dapiColor, '_Tot', sep ='')
    arTot <- paste('Sig_', arColor, dapiColor, '_Tot', sep ='')
    arWindow <- paste('Sig_', arColor, dapiColor, '_Window', sep ='')
    locResults$DAPI <- locResults[,dapiTot]/locResults$n
    locResults$CK <- locResults[,ckTot]/locResults$n
    locResults$CD45 <- locResults[,cd45Tot]/locResults$n
    locResults$AR <- locResults[,arTot]/locResults$n
    locResults$AR.Nuc <- locResults[,arWindow]/locResults$n
    locResults$AR.Cyt <- (locResults[,arTot]-locResults[,arWindow])/locResults$n
    
    locResults$AR.NucPer <- locResults$AR.Nuc/locResults$AR
    locResults$AR.NucPerFlag <- locResults$AR.NucPer > 0.5
    
    CKThresh <- exp(median(log(locResults$CK)) + 3*mad(log(locResults$CK)))
    DAPIThresh <- exp(median(log(locResults$DAPI)) + 3*mad(log(locResults$DAPI)))
    CD45Thresh <- exp(median(log(locResults$CD45)) + 3*mad(log(locResults$CD45)))
    ARThresh <- exp(median(log(locResults$AR)) + 3*mad(log(locResults$AR)))
    
    locResults$DAPI.Flag <- locResults$DAPI >= DAPIThresh
    locResults$CK.Flag <- locResults$CK >= CKThresh
    locResults$CD45.Flag <- locResults$CD45 >= CD45Thresh
    locResults$AR.Flag <- locResults$AR >= ARThresh
    
    summary <- locResults[,c('Id','DAPI','CK','CD45','AR','AR.Nuc','AR.Cyt','AR.NucPer','AR.NucPerFlag','DAPI.Flag','CK.Flag','CD45.Flag','AR.Flag')]

    return(list(thresholds=list(CK=CKThresh, DAPI=DAPIThresh, CD45=CD45Thresh, AR=ARThresh), summary=summary))
}