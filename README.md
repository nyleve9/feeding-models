# feeding-models
QB6 - elk plant intake: bite model vs. density model

The bite model, in which intake rate is regulated by bite mass, offered a better approximation of the Elk feeding data than the density model alternative, in which intake rate is regulated by the interaction between plant density and plant size.  

Table 1 shows the parameter estimates for both models along with the 95% confidence intervals for the parameter estimates for the bite model:

Parameter
Parameter estimate
95% Confidence interval
Bite model
Rmax
37.956
34.124 – 42.377

h
0.0111
0.00790 – 0.0151

standard deviation
8.583
7.586 – 9.817
Density model
Vmax
1,000,000


h
0.0101


standard deviation
16.068



Table 2 shows the corresponding AIC (Akaike Information Criterion), ΔAIC (difference in AIC from the best-fitting model), AICc (AIC values with a correction for small sample size), Akaike weights ("probability" of the model), and negative log likelihoods for each model.  The lowest negative log likelihood was found to correspond with the bite model, and since ΔAIC for the density model is >10 and wr = 1 for the bite model, our analysis shows strong support for the bite model:

Model
-ln(L(θ))
# pars
AIC
ΔAIC
AICc
wr
Bite
413.99
3
833.926
0
834.140
1
Density
486.70
3
979.406
145.480
979.620
0


In addition, plots of observations vs. predictions revealed that the bite model showed a consistent 1:1 relationship between observations and predictions with little bias, whereas the density model consistently underestimated low intake rates and overestimated high intake rates (figure 1):

Figure 2 shows the best-fit model plotted against the data:

Figures 3 and 4 show the negative log likelihood profiles for the 3 parameters included in the bite model, along with corresponding visualizations of the 95% confidence intervals found for each parameter.




Figure 3:

Figure 4:

Model diagnostics show that choosing to use a Normal distribution of residuals (differences between observed intake values and predicted values) was appropriate.  We found this choice reasonable based on the continuous and non-negative nature of the response variable (intake rate) and also based on support provided by the figures 5 - 7:
Figure 5 is a histogram showing that the residuals are normally distributed:

Figure 6 is a plot of the residuals vs. the fitted intake values; it indicates that there is no pattern to the distribution of residuals that would indicate that our model is a poor fit, and there are no outstanding outliers: 

Figure 7 is a Q-Q plot, which compares the distribution of the data to the distribution of predicted values from the bite model.  In general, this Q-Q plot follows the line y=x, which supports that the distributions are similar:

