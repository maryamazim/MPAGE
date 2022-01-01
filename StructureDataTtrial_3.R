# Script to structure trial 3 raw data coming from the machine into R DFs
# Dec 13, 2021
#
### Not using mvfR data ###

# Read in Illumination Data 
illumData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/illumData_3.csv", header=F)
odData <- read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/odData_3.csv", header=F)
concDelMatrix <- as.matrix(read.csv("/Users/mazim/Documents/RahmeLab/MyProj/Analysis/Data/shared_data/concDelMatrix.csv", header=F))
colnames(concDelMatrix) <- 1:12
rownames(concDelMatrix) <- LETTERS[1:8]

dnames <- illumData[,1]
# Function to return a DF with illum data of all replicates and concentrations at a time point!  
get_1TT <- function(d, dnames, concDelMatrix) {
  illum <- d[-(1:3)]
  dnames <- dnames[-(1:3)]
  df <- data.frame(replicate=substr(dnames, 1, 1), 
                   concentration=0,
                   time = d[2],
                   temperature = d[3],
                   illumination = illum)
  for(i in 1:length(dnames)) df$concentration[i] <- gsub(concDelMatrix[substr(dnames[i], 1, 1), substr(dnames[i], 2, nchar(dnames[i]))], pattern=' ', replacement='')
  df
}

# extract and merge illum and OD data, one data frame at a time!
illumDF.all <- NULL
for(i in 2:ncol(illumData)) {
  df <- get_1TT(illumData[,i], dnames=dnames, concDelMatrix=concDelMatrix)
  illumDF.all <- rbind(illumDF.all, df)
}
odDF.all <- NULL
for(i in 2:ncol(odData)) {
  df <- get_1TT(odData[,i], dnames=dnames, concDelMatrix=concDelMatrix)
  odDF.all <- rbind(odDF.all, df)
}

# concentrations only--
illumDF <- illumDF.all[is.element(illumDF.all$concentration, c('0', '2', '10', '20')),]
odDF <- odDF.all[is.element(odDF.all$concentration, c('0', '2', '10', '20')),]

# gene deletions only--
illumDF <- illumDF.all[!is.element(illumDF.all$concentration, c('0', '2', '10', '20')),]
odDF <- odDF.all[!is.element(odDF.all$concentration, c('0', '2', '10', '20')),]

plot_conc <- function(illumDF, conc, average=FALSE, points=FALSE, col='black', ulMargin=0) {
  concDF <- illumDF[illumDF$concentration == conc,]
  if(average == TRUE) {
    x <- tapply(concDF$illumination, concDF$time, mean)
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

# Look at illumination patterns
d <- illumDF
#d <- illumDF[illumDF$time > 500*60,]
plot_conc(d, conc = 0, average = TRUE, col='grey', ulMargin=1000)
plot_conc(d, conc = 2, average = TRUE, points = TRUE)
plot_conc(d, conc = 10, average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 20, average = TRUE, points = TRUE, col='blue')

plot_conc(d, conc = 2, average = TRUE, ulMargin=1000)
plot_conc(d, conc = 10, average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 20, average = TRUE, points = TRUE, col='blue')

# gene deletions
plot_conc(d, conc = 'lasR', average = TRUE, col='black', ulMargin=3000)
plot_conc(d, conc = 'mvfR', average = TRUE, points = TRUE, col='grey')
plot_conc(d, conc = 'pqsBC', average = TRUE, points = TRUE, col='blue')
plot_conc(d, conc = 'pqsE', average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 'pqsH', average = TRUE, points = TRUE, col='pink2')
plot_conc(d, conc = 'pqsL', average = TRUE, points = TRUE, col='purple')
plot_conc(d, conc = 'rhlR', average = TRUE, points = TRUE, col='green')

# Look at growth patterns
d <- odDF
#d <- illumDF[illumDF$time > 500*60,]
plot_conc(d, conc = 0, average = TRUE, col='grey', ulMargin=.2)
plot_conc(d, conc = 2, average = TRUE, points = TRUE)
plot_conc(d, conc = 10, average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 20, average = TRUE, points = TRUE, col='blue')

# gene deletions
plot_conc(d, conc = 'lasR', average = TRUE, col='black', ulMargin=.3)
plot_conc(d, conc = 'mvfR', average = TRUE, points = TRUE, col='grey')
plot_conc(d, conc = 'pqsBC', average = TRUE, points = TRUE, col='blue')
plot_conc(d, conc = 'pqsE', average = TRUE, points = TRUE, col='red')
plot_conc(d, conc = 'pqsH', average = TRUE, points = TRUE, col='pink2')
plot_conc(d, conc = 'pqsL', average = TRUE, points = TRUE, col='purple')
plot_conc(d, conc = 'rhlR', average = TRUE, points = TRUE, col='green')

# Look at data using poly()
# illum data
illumDF$concentration <- factor(illumDF$concentration, levels=c(0,2,10,20))[,drop=T]
illumDF.poly <- lm(illumination ~ replicate + concentration + temperature + poly(time, 3), 
                   data = subset(illumDF, concentration!= '0'))
summary(illumDF.poly)

# OD data
odDF$concentration <- factor(odDF$concentration, levels=c(0,2,10,20))[,drop=T]
odDF.poly <- lm(illumination ~ replicate + concentration + temperature + poly(time, 3), data = odDF)
summary(odDF.poly)

# adjust illumination for growth
illumDF$growth <- odDF$illumination
illumDF.poly <- lm(illumination ~ growth + replicate + concentration + temperature + poly(time, 3), 
                   data = subset(illumDF, concentration!= '0'))
summary(illumDF.poly)

# Swap conc 2 and 20
illumDF$growth <- odDF$illumination
S_illumDF <- illumDF
S_illumDF$concentration <- factor(S_illumDF$concentration, levels=c(0,2,10,20))[,drop=T]
levels(S_illumDF$concentration) <- c("0", "20", "10", "2")
S_illumDF$concentration <- factor(S_illumDF$concentration, levels=c(0,2,10,20))[,drop=T]
S_illumDF.poly <- lm(illumination ~ growth + replicate + concentration + temperature + poly(time, 3), 
                   data = subset(S_illumDF, concentration!= '0'))
summary(S_illumDF.poly)


