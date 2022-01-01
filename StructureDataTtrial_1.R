# Script to structure trial 1 raw data coming from the machine into R DFs
# Dec 12, 2021
#
# Read in Illumination Data 
illumData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/illumData_1.csv", header=F)
odData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/odData_1.csv", header=F)
concMatrix <- as.matrix(read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/concMatrix.csv", header=F))

# Function to return a DF with illum data of all replicates and concentrations at a time point!  
get_1TT <- function(d, concMatrix) {
  illum <- d[-(1:2)]
  design <- matrix(illum, nrow=8, ncol=5, byrow=T)
  design <- design[,-5]
  df <- data.frame(replicate=rep(LETTERS[1:8], each=4), 
                   concentration=c(t(concMatrix)),
                   time = d[1],
                   temperature = d[2],
                   illumination = c(t(design)))
  df
}

# extract and merge illum and od data, one data frame at a time!
illumDF <- NULL
for(i in 2:ncol(illumData)) {
  df <- get_1TT(illumData[,i], concMatrix=concMatrix)
  illumDF <- rbind(illumDF, df)
}
odDF <- NULL
for(i in 2:ncol(odData)) {
  df <- get_1TT(odData[,i], concMatrix=concMatrix)
  odDF <- rbind(odDF, df)
}
names(odDF)[5] <- 'growth'

# Look at data using poly()
# illum data
illumDF$concentration <- factor(illumDF$concentration)[,drop=T]
illumDF.poly <- lm(illumination ~ replicate + concentration + temperature + poly(time, 7), 
                   data = subset(illumDF, concentration!= '0'))
summary(illumDF.poly)

# OD data
odDF$concentration <- factor(odDF$concentration)[,drop=T]
odDF.poly <- lm(growth ~ replicate + concentration + temperature + poly(time, 5), data = odDF)
summary(odDF.poly)

# adjust illumination for growth - GFP replicates and becomes more intense with greater cell density.
# It is imparative to adjust for the part of illumination that comes from growth so that what remains
# is related only to quorum sensing activity.
illumDF$growth <- odDF$growth
illumDF$concentration <- factor(illumDF$concentration)[,drop=T]
illumDF.poly <- lm(illumination ~ poly(growth,3) + replicate + concentration + temperature + poly(time,7), 
                   data = subset(illumDF, concentration!= '0'))
summary(illumDF.poly)

plot_conc <- function(illumDF, conc, average=FALSE, points=FALSE, col='black', ulMargin=0) {
  concDF <- illumDF[illumDF$concentration == conc,]
  if(average == TRUE) {
    x <- tapply(concDF[,5], concDF$time, mean)
    concDFAvg <- data.frame(time = as.numeric(names(x)), illumAvg = x)
    plotRange = c(0, max(concDFAvg$illumAvg)+ulMargin)
    if(points == TRUE) {
      points(concDFAvg$time/60, concDFAvg$illumAvg, col=col, pch=16, ylim=plotRange)
    } else {
      plot(concDFAvg$time/60, concDFAvg$illumAvg, xlab='Time (minutes)', 
           ylab='Illumination (Average)', col=col, pch=16, ylim=plotRange,
           main=paste0('Concentration = ', conc))
    }
  } else {
    plotRange = c(0, max(concDF$illumination)+ulMargin)
    if(points == TRUE) {
      points(concDF$time/60, concDF$illumination, col=col, pch=16, ylim=plotRange)
    } else {
    plot(concDF$time/60, concDF$illumination, xlab='Time (minutes)', 
         ylab='Illumination', col=col, pch=16, ylim=plotRange,
         main=paste0('Concentration = ', conc))
    }
  }
}

d <- illumDF
#d <- illumDF[illumDF$time > 500*60,]
plot_conc(d, conc = 0, average = TRUE, col='grey', ulMargin=1000)
plot_conc(d, conc = 10, average = TRUE, points = TRUE, ulMargin=1000)
plot_conc(d, conc = 20, average = TRUE, points = TRUE, col='red', ulMargin=1000)
plot_conc(d, conc = 50, average = TRUE, points = TRUE, col='blue', ulMargin=1000)

d <- odDF
plot_conc(d, conc = 0, average = TRUE, col='grey', ulMargin=0.2)
plot_conc(d, conc = 10, average = TRUE, points = TRUE)
plot_conc(d, conc = 20, average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 50, average = TRUE, points = TRUE, col='blue')



