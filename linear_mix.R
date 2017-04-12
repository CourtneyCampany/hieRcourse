library(lme4)
library(visreg)
library(lmerTest)
library(car)

# pref data example -------------------------------------------------------

pref <- read.csv("prefdata.csv")

mod1 <- lm(LMA ~ dfromtop * species, data=pref)

visreg(mod1, xvar="dfromtop", by="species", overlay=TRUE)

#relationship within each tree
lmlis1 <- lmList(LMA ~ dfromtop | ID, data=pref)

#plot all fits

#better idea to fit mixed effect model

# random intercept
me1 <- lmer(LMA ~ dfromtop * species + (1|ID), data=pref)

#random intercept and slope
me2 <- lmer(LMA ~ dfromtop * species + (dfromtop|ID), data=pref)
summary(me2) #can see that effect from dfomtop in random is small (test to drop)
anova(me2) #only F statistic, no P

Anova(me2, test="F")

anova(me1, me2) #likelyhood ratio test 
AIC(me1, me2)
#aic is lower (more than 2) for me1 so do not need a random slope

visreg(me1, "dfromtop", by="species", overlay=T)

#model diagnostics
plot(me1)
 #first plot model better predcits smaller values than large values
# so violates equal variation assumpution
qqPlot(resid(me1))
#both diagnostics suggest that a transformation is necessary


me3 <- lmer(log(LMA) ~ dfromtop * species + (dfromtop|ID), data=pref)
plot(me3) #better
qqPlot(resid(me3))
plot(resid(me3)~fitted(me3))
Anova(me3)
VarCorr(me3)

# mouse -------

mouse <- read.csv("wildmousemetabolism.csv")
  mouse$id <- as.factor(mouse$id)
    
plot(jitter(rmr)~temp, pch=21, bg="gold", ylim=c(0, 0.6),data=mouse)

mouse15 <- mouse[mouse$temp==15,]
  mouse15$id <- with(mouse15, reorder(id, rmr, median, na.rm=TRUE))  
  boxplot(rmr~id, data=mouse15)

#ramdom effect is run nested within individual
mouse_mod0 <- lmer(rmr~temp + (1|id/run), data=mouse)
VarCorr(mouse_mod0)
Anova(mouse_mod0)

#center the intercept to a biological meaningful value (31)
mouse$temp31 <- mouse$temp -31

mouse_mod31 <- lmer(rmr~temp31 + (1|id/run), data=mouse)
summary(mouse_mod31)
AIC(mouse_mod0, mouse_mod31) #should be the same

mouse_mod31bm <- lmer(rmr~temp31*bm + (1|id/run), data=mouse)
library(pbkrtest)
KRmodcomp(mouse_mod31, mouse_mod31bm)
#body mass is very sig
Anova(mouse_mod31bm)
VarCorr(mouse_mod31bm) #random effect of id is smaller


#now add random slope 
mouse_mod3 <- lmer(rmr~bm*temp3 + (temp31|id/run), data=mouse)
anova(mouse_mod31bm,mouse_mod3) #model allowing different slopes is better
#so individuals differ in their response to temp
visreg(mouse_mod3, "bm", by='temp31', overlay=TRUE, partial=FALSE)

#test more fixed effects
mouse_mod4 <- lmer(rmr~bm*temp31 + sex+ (temp31|id/run), data=mouse)
mouse_mod5 <- lmer(rmr~bm*temp31 + wheel + (temp31|id/run), data=mouse)
#no effect of sex but there is one for wheel
KRmodcomp(mouse_mod3, mouse_mod4)
KRmodcomp(mouse_mod3, mouse_mod5)
Anova(mouse_mod5, test='F')

#final plot that will predict using random slope and intecepts for each id

#dataframe with each combo of temp and id for the first run
pred_mouse <- expand.grid(temp31=c(-16,-11,0), id=levels(mouse$id), run=1)

library(doBy)
bm_agg <- summaryBy(bm~id, data=mouse, FUN=mean, keep.names = TRUE)

#know predict
pred_mouse <- merge(pred_mouse, bm_agg)
pred_mouse$rmr_pred <- predict(mouse_mod3, pred_mouse)
library(scales)

plot(rmr~jitter(temp), pch=21, bg=alpha("gold", .5), ylim=c(0, 0.6),
     data=mouse[mouse$run==1,])
#add predlines
invisible(lapply(split(pred_mouse, pred_mouse$id), function(x)lines(x$temp31+31,
                       x$rmr_pred, col=alpha("gold", .8))))


# block design ------------------------------------------------------------

litter <- read.csv("litter.csv")
  litter$plot <- as.factor(litter$plot)
  litter$block <- as.factor(litter$block)

library(lubridate)
  litter$date<- mdy(litter$date)
  litter$date2 <- litter$date - as.Date('2006-05-23')

library(lattice)
bwplot(masslost ~ factor(date)| profile:herbicide,data=litter)

#see that some litter bags are lost so the design is unbalanced
ftable(xtabs(~date2 + profile+herbicide, data=litter))

#see manual for examples of Anova, because data is unbalance the choice of 
#test to get type 1 or type 11 test will give diff result (if balanced they
#would be the same)

litter_mod1 <- lmer(masslost ~ date2 +herbicide*profile+(1|block/plot),
                    data=litter) #start with  heirarcy:block then plot
Anova(litter_mod1, test='F')       
VarCorr(litter_mod1)
#check the effect of plot
litter_mod2 <- lmer(masslost ~ date2 +herbicide*profile+(1|block),data=litter)
Anova(litter_mod2, test='F')
VarCorr(litter_mod2)
visreg(litter_mod1, "profile", by="herbicide", overlay=TRUE, 
       cond=list(date2=100))
visreg(litter_mod1, "date2", by='profile', overlay=TRUE)

anova(litter_mod1, litter_mod2)
#zero effect of plot....but is it built into the experimental design? so keep

#likelyhood ratio tests

#remove interaction
litter_mod_noint <- lmer(masslost ~ date2 +herbicide+profile+(1|block/plot),
                    data=litter) 
#1. anova test
anova(litter_mod1,litter_mod_noint)
#Maximum likelyhood choosen auto because effects (fixed or random) are differ

#2. KRmodcomp
pbkrtest::KRmodcomp(litter_mod1, litter_mod_noint)

#3. AIC
AIC(litter_mod1, litter_mod_noint)

#in all three the model with the interaction is a better model
#the interaction in the surface treatment with different herbicides (visreg)

#investigate interaction by creating a unique vairable id
litter$combotrt <- paste(litter$herbicide, litter$profile, sep='-')
  litter$combotrt[litter$profile=='buried'] <- 'buried'
  #this drops the herbicide trt from the buried combo trt as it was the same
  #in visreg (could do pairwise tests for all combos if I wanted)
  litter$combotrt<- as.factor(litter$combotrt)

#no compare two models with the new 3 factor level and one without herbicide
  litter_mod3 <- lmer(masslost ~ date2 +combotrt+(1|block/plot),
                           data=litter) 
  litter_mod4 <- lmer(masslost ~ date2 +profile+(1|block/plot),
                           data=litter) 
  
  AIC(litter_mod2, litter_mod3, litter_mod4)
  #combo trt (interaction) is the best model
  
#repeated measurements--------------
