

# local regression --------------------------------------------------------

hobart <- read.csv("sydney_hobart_times.csv")

with(hobart, plot(Time~Year))


loess1 <- loess(Time ~ Year, data=hobart, span=.8)
loess2 <- loess(Time ~ Year, data=hobart, span=.4)


library(nlshelper)
plot_loess(loess1)
plot_loess(loess2)

#predict
hobart_pred <- predict(loess2, se=TRUE)
#make a CI for the predictions

hobart_pred <- transform(hobart_pred, 
                           lci = fit - 2*se.fit,
                           uci = fit + 2*se.fit, 
                           Year=hobart$Year)


with(hobart_pred, plot(fit~Year))
with(hobart_pred, lines(lci~Year))
with(hobart_pred, lines(uci~Year))



# GAM ---------------------------------------------------------------------


howell <- read.csv('howell.csv')

with(howell, plot(height~age, col=sex, pch=16))

#is sex differ? or where are they different

library(mgcv)
g <- gam(height ~ s(age, by=sex, k=8), data=howell)

summary(g)


library(visreg)
visreg(g, 'age', by="sex", overlay=TRUE)

#can use with mixed effect model or with different distributions
#good for different times series with irregular data



# quantile regression -----------------------------------------------------

choat <- read.csv('Choat_precipP50.csv')

with(choat, plot(P50~annualprecip))

library(quantreg)
fit1 <- rq(P50~ annualprecip, data=choat, tau=.1) #tau=50 is the median regression
fit2 <- rq(P50~ annualprecip, data=choat, tau=.9)

summary(fit1)
summary(fit2)

abline(fit1) #bottom quantile decreaes with less rain
abline(fit2, col='blue') #90% quanitle doesnt change with rain
#careful when choosing quantile, you need significant data at the level


fit1 <- rq(P50~ annualprecip, data=choat, tau=.1) #tau=50 is the median regression
predat <- data.frame(annualprecip = seq(375, 3100, by=10))
pred1 <- as.data.frame(predict(fit1, newdata=predat, interval='confidence'))

#add confidence to the bottom quantile by predicting (above)
lines(predat$annualprecip, pred1$lower, lty=3)
lines(predat$annualprecip, pred1$higher, lty=3)


# non linear quanitle regression ------------------------------------------
library(quantreg)

foot <- read.csv("anthropometry.csv")
foot <- foot[complete.cases(foot),]
foot_male <- foot[foot$gender=="male",]

foot_nlrq <- nlrq(foot_length ~ SSgompertz(age, Asym, b2, b3),
                  tau=.9, data=foot_male)
#90% foot length at a given age

coef(foot_nlrq)


with(foot_male,
     plot(age, foot_length, pch=16, cex=0.2, col="grey",
          xlim=c(0,22), ylim=c(100,320)))

#examine multiple quantiles

# Loop through desired quantiles ...
for(TAU in c(0.1, 0.5, 0.9)){
  # ... use them in a non-linear quantile regression
  fit_mod <- nlrq(foot_length ~
                    SSgompertz(age, Asym, b2, b3),
                  tau=TAU,
                  data=foot_male)
  # Set up x-values to predict over
  xpred <- seq(min(foot$age), max(foot$age), length=101)
  # And add a line to the plot
  lines(xpred, predict(fit_mod, data.frame(age=xpred)))
  # Label the curve on the right (y-value is predicted from the curve)
  text(21, predict(fit_mod, data.frame(age=max(foot_male$age))),
       labels=as.character(TAU))
}