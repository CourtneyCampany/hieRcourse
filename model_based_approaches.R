# using adonis ------------------------------------------------------------
tibplat <- read.csv("tibplat.csv")
#analysis of variance using distance matrices (do not use complex designs or,
# continuous variables, can use with unbalanced designs)
adonis(tibplat[,3:10]~fertilization * enclosure, data=tibplat)

#can use mvabund to generate effect sizes and predict
library(mvabund)

vars <- tibplat[,1:2]
mod <- mvabund(tibplat[,3:10])

plot(mod)
m1 <- manyglm(mod~fertilization*enclosure, data=vars, family='poisson')
plot(m1)
summary(m1)
anova(m1, nBoot=100, p.uni='adjusted')
