R-scripts to test for correlation between traits under directional evolution

They require the packages ape, ade4, caper, lmtest, tseries

All the functions take as arguments a phylogeny "phy" which has to be provided 
as an object of type "phylo" and one or two numeric vectors X (and Y) storing 
the values of the trait(s) at the tips of "phy" (their indices have to be 
consistent with those of the tips of "phy").

Trend-detection tests

SR_trend.test takes as argument a phylogeny phy and a vector X and returns 
Pearson's correlation test between the tip values and their times in order to 
detect a trend in the trait evolution and the results of the diagnostic tests 
for least squares regression validity conditions (Durbin-Watson's, 
Harrison-McCabe's and Jarque-Bera's tests).

hR_trend.test takes as argument phy and X and returns the test of the nullity of the slope
when regressing the independent contrasts on the "h_k", in order to detect a trend 
in the trait evolution. This function also provides the results of the diagnostic tests 
for least squares regression validity conditions (Durbin-Watson's, Harrison-McCabe's and 
Jarque-Bera's tests).

Correlation tests

SR_corr.test, IC_corr.test, PGLS_corr.test, DC_corr.test, PGLSt_corr.test and MR_corr.test 
take as argument a phylogeny phy and two vectors X and Y and implement the SR, IC, DC and MR tests 
as described in the manuscript. They return the corresponding test of the nullity of the slope in the
associated regression in order to test for correlation. These functions also provide the results 
of the diagnostic tests for least squares regression validity conditions (Durbin-Watson's, 
Harrison-McCabe's and Jarque-Bera's tests).


The script in the file "simulations.R" compute all the figures of the manuscript. 
