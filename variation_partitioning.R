library(vegan)

data(varechem)
data(varespec)

str(varechem)
#lets split variation between elemental data and the soil physical properties
df.chem <- scale(varechem[, c("Al", "K", "P")])
df.soil <- scale(varechem[,c("Baresoil", "Humdepth")])
#up to 5 tables below
spec.var <- varpart(varespec, df.chem, df.soil)
plot(spec.var, Xnames=c('chem', 'soil'))
#includes covariation in venn diagram
spec.var #extract r2
spec.chem <- rda(varespec ~ df.chem)
anova(spec.chem) #includes overlap in venn
spec.chem1 <- rda(varespec ~ df.chem + Condition(df.soil))
anova(spec.chem1) #soil chem variables only in plant community

#more comparisons
spec.var2 <- varpart(varespec, df.chem, df.soil, scale(varechem$N))
windows()
plot(spec.var2, Xnames=c('chem', 'soil','N'))


#more complicated example with spatial processess

data(mite)
str(mite)
data("mite.env")
summary(mite.env)
#coords
data(mite.xy)

plot(mite.xy)

#principal coodinates of spatial matrices

mite.pcnm <- as.data.frame(scores(pcnm(dist(mite.xy))))#extract scores
str(mite.pcnm)

mite.var<- varpart(mite, ~ +WatrCont + SubsDens,
                   ~ +Substrate + Shrub, 
                   mite.pcnm,
                   data=mite.env)
windows()
plot(mite.var, Xnames=c('water+density', 'substrate+shrub', "space"))

temp <- cbind(mite.env, mite.pcnm[, c('PCNM1','PCNM2')])
anova(rda(mite~WatrCont+SubsDens + Condition(Substrate) + Condition(Shrub)
          +Condition(PCNM1)+Condition(PCNM2),
          data=temp)) #significance of watercontent and density only



# exercises ---------------------------------------------------------------

endo <- read.csv("endophytes.csv")
endo_env <- read.csv("endophytes_env.csv")
endo_dist <- read.csv("endophytes_dist.csv")

str(endo)
str(endo_env)

endo_chem<- scale(endo_env[,3:4])

endo.var <- varpart(endo, ~species, ~type, endo_chem, data=endo_env)
windows()
plot(endo.var, Xnames=c('species', 'type','chemistry'))
endo.chem1 <- rda(endo ~ species + Condition(type)
                  + Condition(endo_chem), data=endo_env)
anova(endo.chem1)
#species explains the most variation and is significant

endo.pcnm <- as.data.frame(scores(pcnm(dist(endo_dist))))#extract coords

endo.var2<- varpart(endo, ~species, ~type, endo_chem,
                    endo.pcnm,
                   data=endo_env)
windows()
plot(endo.var2,Xnames=c('species', 'type','chemistry', "space")) 
#space not greater than variables

temp_endo <- cbind(endo_env, endo.pcnm[, c('PCNM1','PCNM2')])

anova(rda(endo~PCNM1+PCNM2 + Condition(species) + Condition(type)
          +Condition(endo_chem),
          data=temp_endo)) #space not significant
anova(rda(endo~endo_chem + Condition(species) + Condition(type)
          +Condition(PCNM1)+Condition(PCNM2),
          data=temp_endo)) #chem is sig
anova(rda(endo~species + Condition(endo_chem) + Condition(type)
          +Condition(PCNM1)+Condition(PCNM2),
          data=temp_endo)) #species is sig
anova(rda(endo~type + Condition(species) + Condition(endo_chem)
          +Condition(PCNM1)+Condition(PCNM2),
          data=temp_endo)) #type is sig


library(ade4)
#generate dummy variabels with dudi.hillsmith

# mite data ---------------------------------------------------------------

str(mite)
#capscale performs constrianed pca
mite.hel <- decostand(mite, method='hellinger')
mite.cap1 <-capscale(mite.hel ~ ., data=mite.env, distance='bray')
windows()
plot(mite.cap1)
#arrows associated with continutous
#no arrows for categorical variables

mite.cap0 <- capscale(mite.hel ~ 1, data=mite.env, distance='bray')
step.env <- ordistep(mite.cap0,formula(mite.cap1))
step.env$anova

library(ade4)
mite.env.mod <- dudi.mix(mite.env)
head(mite.env.mod$tab)
mite.env.mod2 <- mite.env.mod$tab

#new columns for different levels of everything categorical
#removes problems with categorical variables
mite.cap0 <- capscale(mite.hel ~ 1, data=mite.env.mod2, distance='bray')
mite.cap1 <-capscale(mite.hel ~ ., data=mite.env.mod2, distance='bray')
plot(mite.cap1)
step.env <- ordistep(mite.cap0,formula(mite.cap1))
step.env$anova


