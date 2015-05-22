pairs.trading.future <- function(trade.pairs, future.cost=0.3e-4, etf.cost=1.5e-4, lookback=45, 
                          entryZscore=1, exitZscore=0, stoplossZscore=4, 
                          lagged=F, method=c("OLS", "TLS", "KF"))
{
  method <- match.arg(method)
  
  index.future <- trade.pairs$asset2
  etf <- trade.pairs$asset1
  
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
    n <- nrow(trade.pairs)  # number of observations
    d <- 1                   # dimension of observations, i.e. index future
    
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
    yt <- matrix(coredata(etf), ncol=n)
    Zt <- array(NA, dim=c(d,m,n))
    Zt[1,,] <- t(coredata(index.future))
    #Zt[1,,] <- t(cbind(1, coredata(trade.pairs$ETF)))
    
    out <- fkf(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt)
    hedge.ratio <- xts(out$at[,2:(n+1)], time(trade.pairs))
    #hedge.ratio <- xts(out$at[2,2:(n+1)], time(trade.pairs))
    yport <- xts(rowSums(trade.pairs * cbind(1, -hedge.ratio)), time(trade.pairs))
    yport[1:lookback] <- NA
  }
  
  # check the change of etf (in the units of volatility)
  etf.change <- (etf-rollapply(etf, 5, mean))/rollapply(etf, 5, sd)
  etf.change[is.na(etf.change)] <- 0
  
  # bollinger bands
  zScore <- (yport-rollapplyr(yport, lookback, mean))/rollapplyr(yport, lookback, sd)
  entryZscore <- 1
  exitZscore <- 0
  stoplossZscore <- 4
  shortsEntry <- coredata(zScore < -entryZscore) & coredata(etf.change > -1)
  shortsExit <- coredata(zScore >= -exitZscore)
  shortsStoploss <- coredata(zScore < -stoplossZscore)
  longsEntry <- coredata(zScore > entryZscore) & coredata(etf.change <= 1)
  longsExit <- coredata(zScore <= exitZscore)
  longsStoploss <- coredata(zScore > stoplossZscore)
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
 
  # trade futures only
  #positions <- matrix(rep(numUnits, ncol(trade.pairs)), ncol=ncol(trade.pairs)) * cbind(1, -hedge.ratio) * trade.pairs
  positions <- numUnits*index.future
  
  if (lagged)
    positions <- lag(positions, 1)
  # percentage weight, used to calculate pnl
  #positions <- positions/rowSums(abs(positions))
  
  N <- nrow(positions)
  ret <- (lag(index.future,-1)-index.future)/index.future
  #ret.net.cost <- cbind(ret[,1]-future.cost, ret[,2]-etf.cost)
  pnl <- rowSums(positions * ret)
 
  net.pnl <- pnl - rowSums(abs(positions) * future.cost)
  pf.ret <- net.pnl/rowSums(abs(positions))
  pf.ret[is.na(pf.ret)] <- 0

  list(pfRet=xts(pf.ret, order.by=time(positions)), numUnits=numUnits.xts, positions=positions)
  
}