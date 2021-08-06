# BBS_seasonal_shift_ALHU
modified BBS model that accounts for changing seasonal patterns in Allen's Hummingbird.

Using the package bbsBayes, this code creates a modified version of the GAMYE model that also estimates a non-linear smooth for the effect of season (day of the season). 
The smooths are allowwed to vary by decade, and the parameters of the smooths are fit with a first-difference time-series, so that the shape of the smooth in each decade is similar to the shape in the decade before. 

