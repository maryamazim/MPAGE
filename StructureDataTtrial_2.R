# Script to structure trial 2 raw data coming from the machine into R DFs
# Dec 18, 2021
#
# Read in Illumination Data 
illumData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/illumData_2.csv", header=F)
odData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/odData_2.csv", header=F)
delMatrix <- as.matrix(read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/delMatrix.csv", header=F,
                                check.names=F, col.names=c(2:5, 7:11), row.names=c('B','C','D','E','F','G')))

# Function to return a DF with illum data of all replicates and types at a time point!  
get_1TT <- function(d, delMatrix) {
  illum <- d[-(1:2)]
  type <- c(t(delMatrix)); type <- type[type != ""]
  df <- data.frame(replicate=rep(row.names(delMatrix), each=6), 
                   type=type,
                   time = d[1],
                   temperature = d[2],
                   illumination = illum)
  df
}

# extract and merge illum and od data, one data frame at a time!
illumDF <- NULL
for(i in 2:ncol(illumData)) {
  df <- get_1TT(illumData[,i], delMatrix=delMatrix)
  illumDF <- rbind(illumDF, df)
}
odDF <- NULL
for(i in 2:ncol(odData)) {
  df <- get_1TT(odData[,i], delMatrix=delMatrix)
  odDF <- rbind(odDF, df)
}
names(odDF)[5] <- 'growth'
illumDF$growth <- odDF$growth

# Look at data using poly()
# illum data
illumDF$type <- factor(illumDF$type)[,drop=T]
illumDF.poly <- lm(illumination ~ -1 + type + temperature + poly(time, 7), 
                   data = illumDF[!is.element(illumDF$type, c('LB', 'WT')),])
summary(illumDF.poly)
newD <- illumDF[!is.element(illumDF$type, c('LB', 'WT')),]
newD$replicate <- 'B'
newD$temperature <- mean(newD$temperature)
newD$growth <- mean(newD$growth)
newD$pred <- predict(illumDF.poly, newdata=newD)

# OD data
odDF$type <- factor(odDF$type)[,drop=T]
odDF.poly <- lm(growth ~ -1 + type + temperature + poly(time, 5), 
                data = odDF[!is.element(odDF$type, c('LB', 'WT')),])
summary(odDF.poly)

# adjust illumination for growth - GFP replicates and becomes more intense with greater cell density.
# It is imarative to adjust for the part of illumination that comes from growth so that what remains
# is related only to quorum sensing activity.
illumDF$type <- factor(illumDF$type)[,drop=T]
illumDF.poly <- lm(illumination ~ -1 + type + poly(growth,3) + temperature + poly(time,7), 
                   data = illumDF[!is.element(illumDF$type, c('LB', 'WT')),])
summary(illumDF.poly)

newD$pred_odAdj <- predict(illumDF.poly, newdata=newD)
plot(newD$time, newD$illumination)  # not adjusted for anything
plot(newD$time, newD$pred)  # Not adjusted for OD (type, temp, poly)
plot(newD$time, newD$pred_odAdj)  # adjusted for everything, including OD


plot_type <- function(illumDF, type, average=FALSE, points=FALSE, col='black', ulMargin=0) {
  typeDF <- illumDF[illumDF$type == type,]
  if(average == TRUE) {
    x <- tapply(typeDF[,5], typeDF$time, mean)
    typeDFAvg <- data.frame(time = as.numeric(names(x)), illumAvg = x)
    plotRange = c(0, max(typeDFAvg$illumAvg)+ulMargin)
    if(points == TRUE) {
      points(typeDFAvg$time/60, typeDFAvg$illumAvg, col=col, pch=16, ylim=plotRange)
    } else {
      plot(typeDFAvg$time/60, typeDFAvg$illumAvg, xlab='Time (minutes)', 
           ylab='Illumination (Average)', col=col, pch=16, ylim=plotRange,
           main=paste0('type = ', type))
    }
  } else {
    plotRange = c(0, max(typeDF[,5])+ulMargin)
    if(points == TRUE) {
      points(typeDF$time/60, typeDF[,5], col=col, pch=16, ylim=plotRange)
    } else {
      plot(typeDF$time/60, typeDF[,5], xlab='Time (minutes)', 
           ylab='Illumination', col=col, pch=16, ylim=plotRange,
           main=paste0('type = ', type))
    }
  }
}

d <- illumDF
#d <- illumDF[illumDF$time > 500*60,]
plot_type(d, type = "WT", average = TRUE, col='black', ulMargin=1000)
plot_type(d, type = "pqsL", average = TRUE, points=TRUE, col='grey', ulMargin=3000)
plot_type(d, type = "pqsE", average = TRUE, points = TRUE, col='red', ulMargin=1000)
plot_type(d, type = "rhlR", average = TRUE, points = TRUE, col='blue', ulMargin=1000)
plot_type(d, type = "lasR1", average = TRUE, points = TRUE, col='green', ulMargin=1000)
plot_type(d, type = "lasR2", average = TRUE, points = TRUE, col='orange', ulMargin=1000)
plot_type(d, type = "pqsH1", average = TRUE, points = TRUE, col='purple', ulMargin=1000)
plot_type(d, type = "pqsH2", average = TRUE, points = TRUE, col='pink2', ulMargin=1000)
plot_type(d, type = "pqsA CTX::PqsE", average = TRUE, points = TRUE, col='lightseagreen', ulMargin=1000)
plot_type(d, type = "mvfR", average = TRUE, points = TRUE, col='yellow', ulMargin=1000)


d <- odDF
plot_type(d, type = "WT", average = TRUE, col='black', ulMargin=.3)
plot_type(d, type = "pqsL", average = TRUE, points = TRUE, col='grey', ulMargin=.3)
plot_type(d, type = "pqsE", average = TRUE, points = TRUE, col='red', ulMargin=.3)
plot_type(d, type = "rhlR", average = TRUE, points = TRUE, col='blue', ulMargin=.3)
plot_type(d, type = "lasR1", average = TRUE, points = TRUE, col='green', ulMargin=.3)
plot_type(d, type = "lasR2", average = TRUE, points = TRUE, col='orange', ulMargin=.3)
plot_type(d, type = "pqsH1", average = TRUE, points = TRUE, col='purple', ulMargin=.3)
plot_type(d, type = "pqsH2", average = TRUE, points = TRUE, col='pink2', ulMargin=.3)
plot_type(d, type = "pqsA CTX::PqsE", average = TRUE, points = TRUE, col='lightseagreen', ulMargin=.3)
plot_type(d, type = "mvfR", average = TRUE, points = TRUE, col='yellow', ulMargin=.3)
plot_type(d, type = "pqsBC", average = TRUE, points = F, col='black', ulMargin=.3)

plot_type(d, type = "LB", average = TRUE, points = F, col='orange', ulMargin=.3)



