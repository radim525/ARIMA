library(stats)
library(nonlinearTseries)
library(gridExtra)
library(data.table)
library(dplyr)
library(xtable)
library(imputeTS)
library(tseries)

setwd("~/")

# adds percentage of NaN values to time series
NAN = function(ts,per) {
set.seed(123)
r <- sample(1:length(ts), size = per*0.01*length(ts), replace = FALSE)
ts[r] = NaN
return(ts)}


# creates plot of the acf of the analysed data with the arima fitting
# and the table with the relevant calculated quantities
# nn - the amount of time series anysed, which should have names L1,..,Ln, maxp,maxq,maxi - maximal values for arima fitting
# the vector names = c("description of L1","description of L2", "description of Ln") should be created in order to appear in the created table describing the generated data

PandT= function(nn,maxp,maxq,maxi){
  
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0("plot.pdf"), height =5*nn, width =5,points=15)
  layout(matrix(c(1:nn),nn,1, byrow=TRUE))
  MAT = matrix( nrow =nn, ncol =10)
  
  for (i1 in 1:nn) {
    L = get(paste0("L",i1))
    x = as.numeric(get(paste0("L",i1)))
    k =is.nan(x)
    i = length(x)
    while( k[i] ==TRUE ){
      i= i-1
    }
    x = x[1:i]
    k =is.nan(x)
    i = 1
    while( k[i] ==TRUE ){
      i= i+1
    }
    x = x[i:length(x)]
    y  = na_interpolation(x, option = "linear")
    final.aic <- Inf
    final.order <- c(0,0,0)
    for (i in 1:maxp) for (j in 1:maxq) for (n in 0:maxi) {
      try( current.aic <- AIC(arima(y, order=c(i, n, j), method="ML")))
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.order <- c(i, n, j)
        try(  final.arma <- arima(y, order=final.order, method="ML"))
        print(final.aic)
        print(final.order)
      }
    }
    MAT[i1,1]  = final.order[1]
    MAT[i1,2]  = final.order[3]
    MAT[i1,3]  = final.order[2]
    MAT[i1,4]  = round(final.aic,2)
    MAT[i1,5]  = length(na.omit(as.numeric(get(paste0("L",i1)))))
    MAT[i1,6]  = length(y)
    
    j1=   try(timeLag(y, technique = "acf", selection.method =  "first.minimum",do.plot = FALSE))
    j2=    try(timeLag(y, technique = "acf",selection.method =  "first.zero" ,do.plot = FALSE))
    j3=    try(timeLag(y, technique = "acf",selection.method =  "first.e.decay",do.plot = FALSE))
    
    if (!is.numeric(j1)) {
      j1=NaN
    }
    if (!is.numeric(j2)) {
      j2=NaN
    }
    if (!is.numeric(j3)) {
      j3=NaN
    }
    #    
    MAT[i1,7]  =j1
    MAT[i1,8]  =j2
    MAT[i1,9]  =j3
    MAT[i1,10]  =adf.test(na.omit(y))$p.value
    
    timeLag(y, technique = "acf",selection.method =  "first.e.decay",main =paste0(names[i1],", p =",final.order[1] ," q =",final.order[3] ," n =",final.order[2]  ))
    legend(x = "topright",          # Position
           legend = c(paste0("acf.first.minimum=",j1), paste0("acf.first.zero=",j2), paste0("acf.first.e.decay=",j3))  )
    
  }
  while (!is.null(dev.list()))  dev.off()
  
  df1 = as.data.frame(MAT)
  colnames(df1) <- c("p","q","i","AIC","len1","len2","acf.first.minimum","acf.first.zero","acf.first.e.decay","adf.p.value")
  
  dd= as.data.frame(names)
  names(dd) = "type"
  DF = cbind(dd,df1)
  names(DF)
  
  png("TAB.png",width = 700, height = 150+25*nn, units = "px", pointsize = 15)
  grid.table(DF)
  dev.off()
  
  write.table( DF, "TAB.dat"  )
}

# same as PandT, but with logarithmic preprocessing

PandTlog= function(nn,maxp,maxq,maxi){
  
  while (!is.null(dev.list()))  dev.off()
  pdf(paste0("plot.pdf"), height =5*nn, width =5,points=15)
  layout(matrix(c(1:nn),nn,1, byrow=TRUE))
  MAT = matrix( nrow =nn, ncol =10)
  
  for (i1 in 1:nn) {
    L = get(paste0("L",i1))
    x = as.numeric(get(paste0("L",i1)))
    k =is.nan(x)
    i = length(x)
    while( k[i] ==TRUE ){
      i= i-1
    }
    x = x[1:i]
    k =is.nan(x)
    i = 1
    while( k[i] ==TRUE ){
      i= i+1
    }
    x = x[i:length(x)]
    y  = na_interpolation(x, option = "linear")
    y = log(y-1.1*min(y))
    final.aic <- Inf
    final.order <- c(0,0,0)
    for (i in 1:maxp) for (j in 1:maxq) for (n in 0:maxi) {
      try( current.aic <- AIC(arima(y, order=c(i, n, j), method="ML")))
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.order <- c(i, n, j)
        try(  final.arma <- arima(y, order=final.order, method="ML"))
        print(final.aic)
        print(final.order)
      }
    }
    MAT[i1,1]  = final.order[1]
    MAT[i1,2]  = final.order[3]
    MAT[i1,3]  = final.order[2]
    MAT[i1,4]  = round(final.aic,2)
    MAT[i1,5]  = length(na.omit(as.numeric(get(paste0("L",i1)))))
    MAT[i1,6]  = length(y)
    
    j1=   try(timeLag(y, technique = "acf", selection.method =  "first.minimum",do.plot = FALSE))
    j2=    try(timeLag(y, technique = "acf",selection.method =  "first.zero" ,do.plot = FALSE))
    j3=    try(timeLag(y, technique = "acf",selection.method =  "first.e.decay",do.plot = FALSE))
    
    if (!is.numeric(j1)) {
      j1=NaN
    }
    if (!is.numeric(j2)) {
      j2=NaN
    }
    if (!is.numeric(j3)) {
      j3=NaN
    }
    #    
    MAT[i1,7]  =j1
    MAT[i1,8]  =j2
    MAT[i1,9]  =j3
    MAT[i1,10]  =adf.test(na.omit(y))$p.value
    
    timeLag(y, technique = "acf",selection.method =  "first.e.decay",main =paste0(names[i1],", p =",final.order[1] ," q =",final.order[3] ," n =",final.order[2]  ))
    legend(x = "topright",          # Position
           legend = c(paste0("acf.first.minimum=",j1), paste0("acf.first.zero=",j2), paste0("acf.first.e.decay=",j3))  )
    
  }
  while (!is.null(dev.list()))  dev.off()
  
  df1 = as.data.frame(MAT)
  colnames(df1) <- c("p","q","i","AIC","len1","len2","acf.first.minimum","acf.first.zero","acf.first.e.decay","adf.p.value")
  
  dd= as.data.frame(names)
  names(dd) = "type"
  DF = cbind(dd,df1)
  names(DF)
  
  png("TABlog.png",width = 700, height = 150+25*nn, units = "px", pointsize = 15)
  grid.table(DF)
  dev.off()
  
  write.table( DF, "TABlog.dat"  )
}


# # # Example # # #


# create set of analysed ARIMA time series, they should have names L1,L2..,Ln
# their names should appear in the names vector

names = c("arma44len100","arma44len500","arima212len100","arima212len2500")

L1 <- arima.sim(list(order = c(4,0,4), ar=c(0.45,0.11,0.240,0.11), ma = c(0.270,0.11,0.11,0.38)), n = 100)
ts.plot(L1)

L1= NAN(L1,10)
ts.plot(L1)

L2 <- arima.sim(list(order = c(4,0,4), ar=c(0.45,0.11,0.240,0.11), ma = c(0.270,0.11,0.11,0.38)), n = 500)
ts.plot(L1)

L2= NAN(L2,10)
ts.plot(L2)

L3 <- arima.sim(list(order = c(2,1,2), ar=c(0.45,0.11), ma = c(0.270,0.22)), n = 100)
ts.plot(L3)

L3= NAN(L3,10)
ts.plot(L3)

L4 <- arima.sim(list(order = c(2,1,2), ar=c(0.45,0.11), ma = c(0.270,0.22)), n = 2500)
ts.plot(L4)

L4= NAN(L4,10)
ts.plot(L4)

# Calls function for calculation of the arima fitting for the generated data, produces table and plot
PandT(4,3,1,1)


# Same as above with the applciation of logairthms as preprocessing of the data
PandTlog(4,3,1,1)

# The results appear in main working directory

