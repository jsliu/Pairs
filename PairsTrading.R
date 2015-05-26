# pairs trading strategy

library(FKF)
library(xts)
library(urca)
library(lattice)
library(highfrequency)
library(PerformanceAnalytics)

half.life <- function(trade.pairs)
{
  tp <- ca.jo(trade.pairs)
  
  hedge.ratio <- tp@V[,1]
  
  yport <- xts(trade.pairs %*% hedge.ratio, time(trade.pairs))
  lagY <- na.omit(lag(yport,1))
  deltaY <- na.omit(yport-lagY) 
  lambda <- lm(deltaY~lagY)$coeff[2]
  half.life <- -log(2)/lambda
  lookback <- round(half.life)
}

pairs.trading <- function(trade.pairs, future.cost=0.3e-4, lookback=45, 
                          entryZscore=1, exitZscore=0, stoplossZscore=4, 
                          lagged=F, method=c("OLS", "TLS", "KF"))
{
  method <- match.arg(method)
  
  if (method == "OLS")
  {
    # using OLS
    hedge.ratio <- rollapplyr(trade.pairs, width=lookback, by.column=F, function(x) lm(x[,1] ~ x[,2])$coef[2])
    yport <- xts(rowSums(trade.pairs * cbind(1, -hedge.ratio)), time(trade.pairs))
  }
  
  if (method == "TLS")
  {
    # Total least square, to avoid asymmetric in OLS
    hedge.ratio <- rollapplyr(trade.pairs, width=lookback, by.column=F, function(x) {r <- princomp(coredata(x)); r$loadings[1,1]/r$loadings[2,1]})
    yport <- xts(rowSums(trade.pairs * cbind(1, -hedge.ratio)), time(trade.pairs))
  }
  
  
  if (method == "KF")
  {
    # using Kalman Filter
    m <- ncol(trade.pairs)-1  # dimension of state variable, i.e intercept + etf
    n <- nrow(trade.pairs)    # number of observations
    d <- 1                    # dimension of observations, i.e. index future
    
    # apply kalman filter to find out hedge ratio
    delta <- 0.01
    err <- 0.001
    
    a0 <- runif(m, -0.5, 0.5)
    P0 <- diag(100, m, m)
    
    dt <- matrix(0, m, 1)
    Tt <- diag(1, m, m)
    HHt <- diag(delta/(1-delta),m,m)
    ct <- matrix(0, d, 1)
    GGt <- diag(err, d, d)
    yt <- matrix(coredata(trade.pairs[,1]), ncol=n)
    Zt <- array(NA, dim=c(d,m,n))
    Zt[1,,] <- t(coredata(trade.pairs[,2]))
    #Zt[1,,] <- t(cbind(1, coredata(trade.pairs$ETF)))
    
    out <- fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)
    hedge.ratio <- xts(out$at[,2:(n+1)], time(trade.pairs))
    #hedge.ratio <- xts(out$at[2,2:(n+1)], time(trade.pairs))
    yport <- xts(rowSums(trade.pairs * cbind(1, -hedge.ratio)), time(trade.pairs))
    yport[1:lookback] <- NA
  }
  
  # bollinger bands
  zScore <- (yport-rollapplyr(yport, lookback, mean))/rollapplyr(yport, lookback, sd)
  entryZscore <- 1
  exitZscore <- 0
  stoplossZscore <- 4
  longsEntry <- coredata(zScore < -entryZscore)
  longsExit <- coredata(zScore >= -exitZscore)
  longsStoploss <- coredata(zScore < -stoplossZscore)
  shortsEntry <- coredata(zScore > entryZscore)
  shortsExit <- coredata(zScore <= exitZscore)
  shortsStoploss <- coredata(zScore > stoplossZscore)
  numUnitsLong <- rep(NA, length(yport))
  numUnitsShort <- rep(NA, length(yport))
  numUnitsLong[lookback+1] <- 0
  numUnitsLong[longsEntry] <- 1
  numUnitsLong[longsExit] <- 0
  numUnitsLong[longsStoploss] <- 0
  numUnitsShort[lookback+1] <- 0
  numUnitsShort[shortsEntry] <- -1
  numUnitsShort[shortsExit] <- 0
  numUnitsShort[shortsStoploss] <- 0
  numUnits <- na.locf(numUnitsLong, na.rm=F)+na.locf(numUnitsShort, na.rm=F)
  numUnits.xts <- xts(numUnits, order.by=time(yport))
  close.idx <- endpoints(numUnits.xts, on="days", k=1)
  
  # close every day
  numUnits[close.idx] <- 0
  
  #positions <- matrix(rep(numUnits, ncol(trade.pairs)), ncol=ncol(trade.pairs)) * matrix(rep(hedge.ratio, nrow(trade.pairs)), nrow(trade.pairs), byrow=T) * trade.pairs
  positions <- matrix(rep(numUnits, ncol(trade.pairs)), ncol=ncol(trade.pairs)) * cbind(1, -hedge.ratio) * trade.pairs
  
  if (lagged)
    positions <- lag(positions, 1)
  # percentage weight, used to calculate pnl
  #positions <- positions/rowSums(abs(positions))

  N <- nrow(positions)
  ret <- (lag(trade.pairs,-1)-trade.pairs)/trade.pairs
  #ret.net.cost <- cbind(ret[,1]-future.cost, ret[,2]-etf.cost)
  pnl <- rowSums(positions * ret)

  net.pnl <- pnl - rowSums(abs(positions) * c(future.cost, future.cost))
  pf.ret <- net.pnl/rowSums(abs(positions))
  pf.ret[is.na(pf.ret)] <- 0

  list(pfRet=xts(pf.ret, order.by=time(positions)), numUnits=numUnits.xts, positions=positions)
  
}

average.length <- function(numUnits)
{
  numUnits <- coredata(numUnits)
  out <- rle(numUnits[!is.na(numUnits)])
  pos.len <- mean(out$lengths[out$values==1])
  neg.len <- mean(out$lengths[out$values==-1])
  
  list(posLen=pos.len, negLen=neg.len)
}


daily.pnl <- function(pf.ret)
{
  daily.pnl <- aggregate(pf.ret, as.Date(time(pf.ret)), sum)
  barplot(daily.pnl)
  
  daily.pnl
}

information.ratio <- function(pf.ret)
{
  num <- mean(unlist(lapply(split(trade.pairs,f="days"), nrow)))
  ir <- mean(pf.ret, na.rm=T)/sd(pf.ret,na.rm=T)*sqrt(num*260)

  ir
}

read.etf <- function(filename, k=5, on="mins")
{
  etf <- read.csv(filename, header=F)
  etf.xts <- xts(etf[,-1], order.by=as.POSIXct(etf[,1], format="%Y-%m-%d %H:%M:%OS"))
  etf.interval <- do.call(rbind, lapply(split(etf.xts, f="days"), aggregatePrice, on=on, k=k, marketopen="09:15:00", marketclose="15:00:00"))
  
  etf.interval
}

read.index.future <- function(years=14, months=2:6, k=5, on="mins")
{
  # read futures contract
  if.interval <- NULL
  years <- 14
  months <- 2:6
  for (i in years)
  {
    for (j in months)
    {
      y <- i
      m <- ifelse(j < 10, paste0("0", j), j)
      
      file <- paste0("data/hf/IF", y, m, ".csv")
      index.future <- read.csv(file, head=F)
      if.xts <- xts(index.future[,-1], order.by=as.POSIXct(index.future[,1], format="%Y-%m-%d %H:%M:%OS"))
      if.interval <- rbind(if.interval, do.call(rbind, lapply(split(if.xts, f="days"), aggregatePrice, on=on, k=k, marketopen="09:15:00", marketclose="15:00:00")))
    }
  }
  
  if.interval
}


# future.lotsize = 300, etf.lotsize=100
create.pairs <- function(asset1.interval, asset2.interval, asset1.lotsize, asset2.lotsize)
{
  asset2.mid <- xts(rowMeans(asset2.interval[,2:3]), order.by=time(asset2.interval))
  asset1.mid <- xts(rowMeans(asset1.interval[,2:3]), order.by=time(asset1.interval))
  
  trade.pairs <- merge(asset1.mid*asset1.lotsize, asset2.mid*asset2.lotsize)
  
  # time filter
  market.open <- c("09:30:00",  "13:00:00")
  market.close <- c("11:30:00", "15:00:00")
  dates <- as.Date(time(trade.pairs))
  
  time.filter <- rep(FALSE, length(dates))
  for (i in 1:length(market.open))
  {
    time.filter <- time.filter | (time(trade.pairs) >= as.POSIXct(paste(dates, market.open[i])) & time(trade.pairs) <= as.POSIXct(paste(dates, market.close[i])))
  }
  
  trade.pairs <- na.omit(trade.pairs[time.filter, ])
  
  colnames(trade.pairs) <- c("asset1", "asset2")
  
  trade.pairs
}