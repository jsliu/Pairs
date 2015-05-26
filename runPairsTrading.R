# test pairs

library(data.table)
library(xts)

rm(list=ls())
setwd("/Users/jason/Documents/code/Pairs")
source("PairsTrading.R")
#source("PairsTradingFuture.R")
options(digits.secs=6)



filename <- "/Users/jason/Documents/code/testPairs/cffex.csv"
data <- fread(filename)
setkey(data,TradingDay,UpdateTime)
data[,DateTime:=strptime(paste(TradingDay,UpdateTime),"%Y%m%d%H:%M:%S")+UpdateMillisec/1000]
setkey(data,InstrumentID)

# ZHONG ZHENG 500
IC.dt <- data["IC1505.CFFEX",] 
IC.xts <- xts(IC.dt[,LastPrice],IC.dt[,DateTime])
names(IC.xts) <- "ZhongZheng500"
# SHANG ZHENG 50
IH.dt <- data["IH1505.CFFEX",]
IH.xts <- xts(IH.dt[,LastPrice],IH.dt[,DateTime])
names(IH.xts) <- "ShangZheng50"
# HU SHENG 300
IF.dt <- data["IF1505.CFFEX",]
IF.xts <- xts(IF.dt[,LastPrice],IF.dt[,DateTime])
names(IF.xts) <- "HuShen300"

# inital check
basket <- cbind(IC.xts,IH.xts,IF.xts)
day.close <- endpoints(basket,on="days")
basket.daily <- basket[day.close[-1]]
tsRainbow <- rainbow(ncol(basket.daily))
plot.zoo(basket.daily,screens=1,col=tsRainbow)

# pairs trading
future.cost <- 0.3e-4
k <- 10;
on <- 'secs'
pfRet <- rep(0,length(day.close)-1)

for (i in 2:length(day.close))
{
  basket.day <- basket[(day.close[i-1]+1):day.close[i],]
  freq <- endpoints(basket.day,on=on,k=k)
  nFreq <- length(freq)
  freq.close <- freq[-c(1,(nFreq-3):nFreq)]
  basket.day.min <- basket.day[freq.close,]
  
  trade.pairs <- basket.day.min[,c(1,3)]
  
  out <- pairs.trading(trade.pairs, future.cost=future.cost, method="KF")
  pfRet[i-1] <- sum(out$pfRet, na.rm=T)
}

# result
x11()
charts.PerformanceSummary(xts(pfRet,time(basket.daily)), wealth.index=T)

#average.length(out$numUnits)
