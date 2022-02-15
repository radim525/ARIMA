# ARIMA
testing the quality of ARIMA fitting
This code in R serves for checking of the quality of the ARIMA fitting, it contains 3 functions made composed from various functions
from selected libraries.

NAN(ts,per) is the function for disturbing given time series, it adds percentage of NaN values to the time series 


PandT(nn,maxp,maxq,maxi) this function calculates the p,q,i, values of ARIMA fitting for the time series, which should have names
L1,..,Ln, while nn, stands for the amount n of analysed data, maxp, maxq, maxi, are the maximal values of the parameter space for which the best AIC is searched as the choosing criteria.
The function does produce pdf plot whith the ACF plots of analysed data, where the name of time sereis and the fitts appear above the plot, the legend inside the plot shows 3 automatic choices of the ACF values according to the function acf from the library NonlinearTseries.
The table is saved in dat and as png file showing calculated quantities
name of time series denoted as type column, p, q, i , values of ARIMA fitting, len1 length of the input time series, len2 length of the time sereis with the interpolated missing values, the values of the ACF values of automatic choice, the p value of
the Augemented Dickey Fuller test from the library tseries, where the alternative hypotesis is stationary

