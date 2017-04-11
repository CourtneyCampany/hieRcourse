#two table analyses

#RDA same data assumptions as pca

library(vegan)
data(varechem)
data(varespec)
rda1 <- rda(varespec~scale(varechem), scale=T) 
#scale works only for response spec data
windows()
plot(rda1)

summary(rda1)
#contraned and unconstrained components
#65% of inertia explained in 2nd table
#look at importance of components
#23% explained by first two axises

#14 variables explained the 65% in the species data, then the pca axises explain the rest
#biplot explaines the arrows on the plot

#is that explained variation significant (each of the 14 variables is a df)

anova(rda1)
#lots of residual df, meaning several variables are explained the same variation or
#or not explaining enough to validate keeping

chem.scl <- as.data.frame(scale(varechem))
rda0 <- cca(varespec~ ., data=chem.scl, scale=T) 
rda0a <- cca(varespec~ 1, data=chem.scl, scale=T) #intercept only
step.chem <- ordistep(rda0a, scope=formula(rda0), permutations=how(nperm=999))
step.chem$anova
plot(step.chem)
anova(step.chem)
summary(step.chem) #31% of variation by the 3 variables that we kept


# exercises ---------------------------------------------------------------

endo <- read.csv("endophytes.csv")
endo_env <- read.csv("endophytes_env.csv")
?capscale

endo_envss <- endo_env[,3:5]

endo1<- dbrda(endo~ ., data=end_envss) 
windows()
plot(endo1)

#not sure how this function below works, adds lines to plots?
fit <- envfit(endo, endo_envss)
plot(fit)

endo_full <- dbrda(endo~ ., data=endo_envss, scale=T) 
endo0 <- dbrda(endo~ 1, data=endo_envss, scale=T) #intercept only
step.endo <- ordistep(endo0, scope=formula(endo_full), permutations=how(nperm=999))
step.endo$anova
plot(step.endo)
anova(step.endo)

