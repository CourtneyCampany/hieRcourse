#pca response variables that are linear and continutous (pca and rda)
#if not linear can you tranform and then use pca?
#with count data that have a minimal about of zeros then cca

library(vegan)
data(varechem)
str(varechem)

#pca analysis (good for continuous data)

#variance of variables
var(varechem$Ca)
var(varechem$K)
#variance way higher in Ca compared to K, so need to scale variables (scale =T in rda())

chem1 <- rda(varechem, scale=T) #mean is 0 and SD=1 with scale=T

windows()
plot(chem1)

data("varespec")
str(varespec)

spec1 <- rda(varespec, scale=T)
windows()
plot(spec1)
summary(spec1)#eigenvalues explain amnt of variation for each axis

#seems to violate assumptions of normal residuals (far right variables)
#so we will use correspondace analysis which is better for count data (chi2)
spec2 <- cca(varespec)
windows()
plot(spec2)

#when trying to decide which analysis to use you can look at the gradient length of the 
#data using the decorana() function, axis length <3 choose pca

#principal coordinate analysis----------------------------------

spec3 <- wcmdscale(vegdist(varespec, method='bray'), eig=T)
windows()
plot(spec3)
#look at eigenvalues
spec3$eig
eigens <- spec3$eig[spec3$eig >=0]
eigens/sum(eigens)

#may need to look at ?decostand for ecological standardizations
sptrans <- decostand(varespec, "normalize") #applys across datasetcourtc

spec3.scrs <- scores(spec3,display='sites') 
#sites refers to rows, #species refers to columns but no colums in pcoa
str(spec3.scrs)
windows()
plot(spec3.scrs[,1:2], pch=18, col="blue")


#exercis
allom <- read.csv("Allometry.csv")
allom1 <- rda(allom[,-1])
windows()
plot(allom1)
var(allom$branchmass)
var(allom$leafarea)
allom2 <- rda(allom[,-1], scale=T)
windows()
plot(allom2)

allomtrans <- decostand(allom[,-1], method="range")
allom3 <- rda(allomtrans, scale=T)
windows()
plot(allom3)

#row 41 seems to be bad
#look at NMDS for ranks for possible outliers (section 1.2.4 in manual)



