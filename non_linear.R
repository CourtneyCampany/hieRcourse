

# nls examples ------------------------------------------------------------


chapm <- function(x, Asym, b,c) Asym*(1-exp(-b*x))^c

curve(chapm(x, Asym=100, b=.15, c=3), from =0, to=50)
curve(chapm(x, Asym=100, b=.05, c=3), add=TRUE)


data("Loblolly")

with(Loblolly, plot(age, height))

#model
nls_lob <- nls(height~ chapm(age, Asym, b, c), 
               data=Loblolly,
               start=list(Asym=100, b=.1, c=2.5))

library(nlstools)
overview(nls_lob)

plot_nls(nls_lob, ylim=c(0,80), xlim=c(0,30))

#cannot extract coefs with overview, need to use broom
library(broom)
tidy(nls_lob, conf.int=TRUE) #can get a few relevant stats from the model
glance(nls_lob)
coef(nls_lob)


nlsfitSS <- nls(height ~ SSgompertz(age, Asym, b2,b3),
                data=Loblolly)
tidy(nlsfitSS)

#diagnostics (similar model assumptions, equal variance)

library(car)
qqPlot(residuals(nlsfitSS))

plot(fitted(nlsfitSS), residuals(nlsfitSS))
abline(h=0) #not good 

plot(Loblolly$height, fitted(nls_lob), pch=16, col="dimgrey")
abline(0,1)

#goodness of fit, cannto get r2
glance(nls_lob)$sigma #residual standard error, se of residuals


library(nlshelper)
plot_nls(nlsfitSS)

#predictiosn from non linear model
predict(nlsfitSS, newdata = data.frame(age=12))

newdat <- data.frame(age=c(18, 30))
# Predict from the fitted model for the new dataframe:
predict(nls_lob, newdata=newdat)

#if you want to evaluate a parameter you can have two models and use AIC to check
#drop the C from the first model for example

fit1 <- nls(height~ Asym*(1-exp(-b*x))^c, 
               data=Loblolly,
               start=list(Asym=100, b=.1, c=2.5))

fit2 <- nls(height~ Asym*(1-exp(-b*x)), 
            data=Loblolly,
            start=list(Asym=100, b=.1, c=2.5))

AIC(fit1, fit2)
anova(fit1, fit2)
#we need C which gives us the inflection point

#visualize the correlation between the parameters using 'which'
library(ellipse)

plot(ellipse(nls_lob, which=c("Asym","b")), type='l')

# Optionally, add the estimated coefficients as a point
p <- coef(nls_lob)
points(p["Asym"], p["b"], pch=19)
#possible to plot all confidence ellipses with the plotcorr function, 
#for a quick visualizationof correlation and confidence regions of the parameters,:
plotcorr(summary(nls_lob, correlation=TRUE)$correlation)
#in a matrix format, positive or negative and width of cirlce is how correlated


# foot lenth example ------------------------------------------------------

foot <- read.csv("anthropometry.csv")
  foot <- foot[complete.cases(foot),]

#fit grouping variables to fit different variables using nlsList
#only use selfstarting functions!!!!!

# Note the grouping operator '|' to specify to fit the curve by gender
fit_foot_gender <- nlsList(foot_length ~
                             SSgompertz(age, Asym, b2, b3) | gender,
                           data=foot)

# Load nlshelper to allow a tidy statement of the fit
library(nlshelper) #load this before broom
broom::tidy(fit_foot_gender, conf.int=TRUE)
#CI do not overlap, probably different

palette(c("forestgreen", "cornflowerblue"))
plot_nls(fit_foot_gender, pch=16, cex=.2)


fit_foot_0 <- nls(foot_length ~
                             SSgompertz(age, Asym, b2, b3),
                           data=foot)

anova_nlslist(fit_foot_gender, fit_foot_0)


tidy(fit_foot_0, conf.int = TRUE)
confint(fit_foot_0) #same method they are the same

library(nlstools)
confint2(fit_foot_0)
#this will work when others do not converge

#table of stats
cbind(Estimate=coef(fit_foot_0), confint2(fit_foot_0))



# chickweight exercise ----------------------------------------------------

chickweight <- as.data.frame(ChickWeight)

chickss <- chickweight[chickweight$Chick==6 & chickweight$Time >0,]
with(chickss, plot(weight~Time))

chickfit1 <- nls(weight ~SSgompertz(Time, Asym, b2, b3),data=chickss)

chickfit2 <- nls(weight ~ SSweibull(Time, Asym, Drop,lrc,pwr), data=chickss)

chickfit3 <- nls(weight ~SSlogis(Time, Asym, xmid, scal), data=chickss)

AIC(chickfit1, chickfit2, chickfit3)
#wiebul fits best

plot_nls(chickfit1, pch=16, cex=1)
plot_nls(chickfit2, add=TRUE, col="blue")
plot_nls(chickfit3, add=TRUE, col="red")

#what is the weight of the chick at t=25?
predict(chickfit2, newdata = data.frame(Time=25))


# enzyme exercise ---------------------------------------------------------

data('Puromycin')

with(Puromycin, plot(rate~conc, pch=16,col=state))
with(Puromycin, plot(log(rate)~log(conc), pch=16,col=state)) #log does not make linear, need nls

enz1 <- nls(rate ~SSmicmen(conc, Vm,k),data=Puromycin)
plot_nls(enz1, pch=16, cex=1)

library(nlshelper) #load this before broom
broom::tidy(enz1, conf.int=TRUE)


