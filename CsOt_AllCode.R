

##last modified - 1/26/2023
###############################################
## The following is the complete analysis for Dylan Keel's master's thesis
## complete Spring 2023
## The analysis will first fit and select models for C.shasta, then Chinook
## Summary tables, statistics, and figures are at the end of the code


#library required packages

library(doSNOW)
library(doParallel)
library(doMPI)
library(tidyverse)
library(tictoc)
library(fields)
library(ggplot2)
library(sp)
library(lattice)
library(INLA)
library(plotly)
library(geometry)
library(viridis)
#library(devtools) #might need devtools to install some packaages
library(tictoc)
library(kableExtra)

library(rgdal)

#sometimes this happens...
#install.packages("gstat",force=TRUE)
library(gstat)
library(remotes)
library(INLAOutputs)


library(inlamesh3d)
library(inlatools)
library(INLAutils)
library(ggregplot)

#getwd()



###############Chinook GLMM analysis ####

## my working directory
setwd("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis")
# source useful scripts
source("HighstatLibV13.R")
source("INLA_plotting_functions.R")


## Read in all the data
CSD2=read.csv("CSD.csv")


## Find the upper lower and mean concentrations per liter
CSD2$C.shasta.vol.cor.conc=replace_na(round(CSD2$C.shasta.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc=replace_na(round(CSD2$Chinook.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.upper=replace_na(round(CSD2$C.shasta.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.upper=replace_na(round(CSD2$Chinook.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.lower=replace_na(round(CSD2$C.shasta.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.lower=replace_na(round(CSD2$Chinook.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


## look at the data
par(mar=c(2,2,2,2))
hist(CSD2$C.shasta.vol.cor.conc,breaks=50)
hist(CSD2$Chinook.vol.cor.conc,breaks=50)

## check for multicolinearity
MyVar <- c("D.thal.m","DistanceFromShore.m.z","D.M.")

Mypairs(CSD2[,MyVar])

# procede when finding none


#scale the covariates

CSD2.scale <- dplyr::select(CSD2, all_of(MyVar)) 
names(CSD2.scale) <- paste(MyVar, '.z', sep='')

scale.par = do.call(data.frame, 
                    list(mean = apply(CSD2.scale, 2, mean),
                         sd = apply(CSD2.scale, 2, sd),
                         min = apply(CSD2.scale, 2, min),
                         max = apply(CSD2.scale, 2, max)))

CSD2.scale <- as.data.frame(scale(CSD2.scale))
CSD2 <- cbind(CSD2, CSD2.scale)


#now check all the covariates for multicolinearity

CSD2$Distance.Left.Bank.m[CSD2$Distance.Left.Bank.m==0]<-0.01
MyVar <- c("DistanceFromShore.m.z", "D.M.","V.M.","Type01","SiteN","D.thal.m")
Mypairs(CSD2[,MyVar])

#colinear variables are depth and type, and distance from shore and velocity




# plot them against the response variable
MyVar<-c("D.M..z","DistanceFromShore.m.z.z","V.M..z","Type01","SiteN","D.thal.m.z")
#standard
Cs.Copies=CSD2$Chinook.vol.cor.conc
MyMultipanel.ggp2(Z = CSD2, 
                  varx = MyVar, 
                  vary = "Chinook.vol.cor.conc", 
                  ylab = "Response variable",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


#make sure data is in the right format
CSD2$Site<-as.factor(CSD2$Site)
CSD2$SiteN<-as.factor(CSD2$Site)
CSD2$SiteN<-recode_factor(CSD2$SiteN,KLM=1,I5=2,TOH=3,BVR=4,KMN=5,SV=6)
Site=as.numeric(CSD2$SiteN)
CSD2=subset(CSD2,select=-Site)
D.M..z=as.numeric(CSD2$D.M..z)
DistanceFromShore.m.z.z=as.numeric(CSD2$DistanceFromShore.m.z.z)
Type01=as.factor(CSD2$Type01)
CSD2$Depth2=as.numeric(CSD2$D.M..z)*as.numeric(CSD2$D.M..z)
Depth2=CSD2$Depth2

#make base model forms
base.form <- Cs.Copies ~ D.M..z + DistanceFromShore.m.z.z +  f(Site,model="iid")
base.form.v <- Cs.Copies ~ D.M..z + V.M..z +  f(Site,model="iid")
base.form.t <-Cs.Copies ~ Type01 + DistanceFromShore.m.z.z +  f(Site,model="iid")
base.form.t.v <- Cs.Copies ~ Type01 + V.M..z +  f(Site,model="iid")
base.form.thal.d<-Cs.Copies ~ D.M..z + D.thal.m.z +  f(Site,model="iid")
base.form.thal.t<-Cs.Copies ~ Type01 + D.thal.m.z +  f(Site,model="iid")
base.form.thal<-Cs.Copies ~ D.thal.m.z +  f(Site,model="iid")

pcprior <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))

base.form.d2<-Cs.Copies~Depth2+D.M..z+DistanceFromShore.m.z.z+  f(Site,model="iid",hyper=pcprior)
base.form.d<-Cs.Copies~D.M..z+  f(Site,model="iid",hyper=pcprior)
base.form.t<-Cs.Copies~Type01+f(Site,model="iid",hyper=pcprior)

hist(log(CSD2$Chinook.vol.cor.conc),breaks=100)

zeros=which(CSD2$Chinook.vol.cor.conc==0)
length(zeros)

######here we fit some models without spatial effects
####normal non-spatial with distance from shore


nonspat.d2 <- inla(base.form.d2,
                
                control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                family = "poisson", 
                
                control.predictor = list(link = 1, compute = TRUE),
                data = CSD2)
summary(nonspat.d2)



nonspat.d2.nbn <- inla(base.form.d2,
                   
                   control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                   family = "nbinomial", 
                   
                   control.predictor = list(link = 1, compute = TRUE),
                   data = CSD2)
summary(nonspat.d2.nbn)


### Check to see if zero inflation is warranted

nonspat.d2.zip <- inla(base.form.d2,
                   
                   control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                   family = "zeroinflatedpoisson1", 
                   
                   control.predictor = list(link = 1, compute = TRUE),
                   data = CSD2)
summary(nonspat.d2.zip)



nonspat.d2.zinbn <- inla(base.form.d2,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "zeroinflatednbinomial1", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = CSD2)
summary(nonspat.d2.zinbn)


###########################Check for Overdispersion

source("Modified Distribution Check inlatools.R")
source("Modified Dispersion Check inlatools.R")
# 
# bind_rows(
#   nonspat.d2.nbn$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2.zinbn$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = TRUE
#     ),
#   nonspat.d2.zip$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = TRUE
#     )
# ) %>%
#   select(
#     distribution, zeroinflation, parameter, 
#     mean, lcl = `0.025quant`, ucl = `0.975quant`
#   ) -> summary_zip_fixed
# bind_rows(
#   nonspat.d2.nbn$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2.zinbn$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = TRUE
#     ),
#   nonspat.d2.zip$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = TRUE
#     )
# ) %>%
#   select(
#     distribution, zeroinflation, parameter, 
#     mean, lcl = `0.025quant`, ucl = `0.975quant`
#   ) %>%
#   arrange(parameter, distribution, zeroinflation) -> summary_zip_hyper
# 
# # 
# summary_zip_hyper
# summary_zip_fixed
# 
# nonspat.d2_dc <- dispersion_check(nonspat.d2)
# nonspat.d2.nbn_dc <- dispersion_check(nonspat.d2.nbn)
# nonspat.d2.zip_dc <- dispersion_check(nonspat.d2.zip)
# nonspat.d2.zinbn_dc <- dispersion_check(nonspat.d2.zinbn)
# 
# plot(nonspat.d2_dc)
# plot(nonspat.d2.nbn_dc)
# plot(nonspat.d2.zip)
# plot(nonspat.d2.zinbn)

#Cs.Copies

# 
# 
# list(
#   poisson = nonspat.d2,
#   negbin = nonspat.d2.nbn,
#   zipoisson = nonspat.d2.zip,
#   zinegbin = nonspat.d2.zinbn
# ) %>%
#   fast_distribution_check() -> zip_fdc
# 
# plot(zip_fdc, scales = "free")
# 



# 
# 
# nonspat <- inla(base.form,
#                 
#                 control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
#                 family = "poisson", 
#              
#                 control.predictor = list(link = 1, compute = TRUE),
#                 data = CSD2)
# summary(nonspat)
# 
# sum(nonspat$cpo$cpo)

# 
# 
# #### normal non-spatial with velocity
# nonspat.v <- inla(base.form.v,
#                   
#                   control.compute = list(waic = TRUE),
#                   family = "poisson", 
#                   Ntrial = 1,
#                   
#                   data = CSD2)
# summary(nonspat.v)
# 
# nospat.thal<- inla(base.form.thal.d,
#                    
#                    control.compute = list(waic = TRUE),
#                    family = "poisson", 
#                    Ntrial = 1,
#                    
#                    data = CSD2)
# summary(nospat.thal)
# 
# # 
# # #with type
# nonspat.t <- inla(base.form.t,
#                   
#                   control.compute = list(waic = TRUE),
#                   family = "poisson",
#                   Ntrial = 1,
#                   
#                   data = CSD2)
# summary(nonspat.t)

# check out the models
par(mar=c(3,3,3,3))

observed<-Cs.Copies

D<-INLAutils::plot_inla_residuals(nonspat.d2.nbn, observed)


P<-autoplot(nonspat.d2.zinbn,which=(c(1:5)))
P

summary(nonspat.d2.zinbn)
summary(nonspat.d2.nbn)

#average distance between sites
summary(CSD2$Meters)
(300533-212602)/5

## convert latitude and longitude to UTM coordinates in km

xy <- with(CSD2, data.frame(Sample.ID, Meters/1000/111.111, Distance.Left.Bank.m/1000/111.111))

coordinates(xy) <- c("Meters.1000.111.111", "Distance.Left.Bank.m.1000.111.111")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
utm<-NULL
utm <- data.frame(spTransform(xy, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")))
names(utm) <- c('Sample.ID', 'Xkm', 'Ykm', 'drop')

#View(utm)
CSD2<- merge(CSD2, utm)
#CSD2<-CSD2[1:211,]
par(mfrow=c(1,1))
#plot(utm$Ykm~utm$Xkm)

#nonspat w/out velocity adj
Pi1 <- nonspat.d2.nbn$summary.fitted.values[, "mean"]

D1  <- (Cs.Copies - Pi1) / sqrt(Pi1)

summary(D1)

D1
# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink x 2 orders of magnitude
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(utm$Xkm),
                     Y = as.numeric(utm$Ykm))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 5 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)




## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 0.05,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(Km)") + ylab("Sample variogram")
p1 <- p + theme(text = element_text(size=15))+ylim(0,225)
p1 

#with a sill of ~165 and a nugget of ~110 we find:

Sill.Nug.R=110/165
1-Sill.Nug.R

#approximately 1/3 of the variance arises from spatial structure



#we see that spatial autocorrelation increases from 0 to 60 meters
# in our data

Loc<-NULL
Loc <- cbind(utm$Xkm,utm$Ykm)
Loc
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (km)",
     ylab = "Frequency")
#View(Loc)


RangeGuess <-0.03

#max.edge = diff(range(Loc[1]))/15


#tic()
#mesh1 = inla.mesh.2d(loc=Loc,
#  max.edge = max.edge)
#toc()

#plot(mesh1, main="1st attempt"); points(Loc, col="blue")
#names(Loc)


###make the mesh

#plot(Loc[,2]~Loc[,1])


# There is still a spatial correlation in the residuals of spat.inla.  The first thing that you should try is to
# use a mesh with more triangles. Use a smaller MaxEdge and also a smaller
# cutoff. By doing that we allow for smaller-scale correlation.
# Right now we have:

require(splancs)


# 
# Hull <- inla.nonconvex.hull(Loc, convex = -0.085)
# MaxEdge  <- RangeGuess /6
# mesh2d     <- inla.mesh.2d(boundary = Hull,
#                            max.edge = c(0.01,10000) * MaxEdge,
#                            cutoff = MaxEdge ,
#                            max.n=4000)  #control max.n to be <3000
# 

# Hull <- inla.nonconvex.hull(Loc, convex = -0.083)
# MaxEdge  <- RangeGuess 
# mesh2d     <- inla.mesh.2d(boundary = Hull,
#                            max.edge = c(0.0001,100000) * MaxEdge,
#                            cutoff = MaxEdge/1000 ,
#                            max.n=4000)  #control max.n to be <3000



Hull <- inla.nonconvex.hull(Loc, convex = -0.085)
MaxEdge  <- RangeGuess /6
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(0.01,10000) * MaxEdge,
                           cutoff = MaxEdge ,
                           max.n=1500)  #control max.n to be <3000
# 
mesh2d$n
# 
# Hull <- inla.nonconvex.hull(Loc, convex = -0.085)
# MaxEdge  <- RangeGuess 
# mesh2d     <- inla.mesh.2d(boundary = Hull,
#                            max.edge = c(0.001,4000) * MaxEdge,
#                            cutoff = MaxEdge / 5,
#                            max.n=1000)  #control max.n to be <3000
# 
# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1)



# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess/10, 0.005), 
                            prior.sigma = c(.1, .05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
#w.index




all.form= Cs.Copies ~ Depth2+ D.M..z + DistanceFromShore.m.z.z + Type01 + V.M..z + 
  D.thal.m.z + f(Site, model = "iid")


stackform<-as.formula(all.form)


modrando<-"Site"

terms <- all.form[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

# put together a dataframe of values for the stack
Xm <- model.matrix(stackform, data = CSD2)
Xm <- as.data.frame(Xm)
N <- nrow(CSD2)



StackFit <- inla.stack(
  remove.unused=TRUE,
  tag = "Fit",
  data = list(Copies=Cs.Copies), 
  A = list(1, 1, A),                  
  effects = list(   
    Intercept = rep(1, N),
    Xm        = Xm[,-1],    #Covariates without the intercept
    w         = w.index))

inla.stack.data(StackFit)

sf.d2.d.dm<-Copies~-1+Depth2+D.M..z+DistanceFromShore.m.z.z+ 
  f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)

sf.d<-Copies~-1+D.M..z+ f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)

sf.t<-Copies~-1+Type01+ f(w, model = spde)
#0
spat.inla.d<- inla(formula = update(sf.d, . ~ . + f(Site, model ="iid")),
                   family = "nbinomial",
                   data = inla.stack.data(StackFit),
                   control.compute = list(waic = TRUE, config=TRUE),
                   control.predictor = list(link = 1, compute = TRUE,
                                            A = inla.stack.A(StackFit)))




spat.inla.d2.d.m<- inla(formula = update(sf.d2.d.dm, . ~ . + f(Site, model ="iid")),
                        family = "nbinomial",
                        data = inla.stack.data(StackFit),
                        control.compute = list(waic = TRUE, config=TRUE),
                        control.predictor = list(link = 1, compute = TRUE,
                                                 A = inla.stack.A(StackFit)))



# 
GRSpred <- inlaPredict(spat.inla.d2.d.m,  mdl.stack = StackFit, spatialMod=T,orig.dat=CSD2,data.pts=192)

dc.siddm=dispersion_check(spat.inla.d2.d.m)
plot(dc.siddm)

dist.siddm=distribution_check(spat.inla.d2.d.m)
plot(dist.siddm)



# spat.inla.d<- inla(formula = update(sf.d2.d.dm, . ~ . + f(Site, model ="iid")),
#                    family = "nbinomial",
#                    data = inla.stack.data(StackFit),
#                    control.compute = list(waic = TRUE, config=TRUE),
#                    control.predictor = list(link = 1, compute = TRUE,
#                                             A = inla.stack.A(StackFit)))
# 
# 
# 
# spat.inla.t<- inla(formula = update(sf.t, . ~ . + f(Site, model ="iid")),
#                    family = "nbinomial",
#                    data = inla.stack.data(StackFit),
#                    control.compute = list(waic = TRUE, config=TRUE),
#                    control.predictor = list(link = 1, compute = TRUE,
#                                             A = inla.stack.A(StackFit)))

#summary(nonspat.d2.nbn)
#2


#dc2=distribution_check(spat.inla.d2.d.m)

#plot(dc2)
#mean(dc2$median)

autoplot(spat.inla.d2.d.m)

INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)

par(mar=c(4,4,4,4))
D<-INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)


## convert latitude and longitude to UTM coordinates in km
xy <- with(CSD2, data.frame(Sample.ID, Meters/1000/111.111, Distance.Left.Bank.m/1000/111.111))

coordinates(xy) <- c("Meters.1000.111.111", "Distance.Left.Bank.m.1000.111.111")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
utm<-NULL
utm <- data.frame(spTransform(xy, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")))
names(utm) <- c('Sample.ID', 'Xkm', 'Ykm', 'drop')

#View(utm)
CSD2<- merge(CSD2, utm)
#CSD2<-CSD2[1:211,]
par(mfrow=c(1,1))
plot(utm$Ykm~utm$Xkm)


length(spat.inla.d2.d.m$summary.fitted.values[,"mean"])
#nonspat w/out velocity adj
Pi1 <- spat.inla.d2.d.m$summary.fitted.values[, "mean"]


fitIndex <- inla.stack.index(StackFit, tag='Fit')$data
fitted <- spat.inla.d2.d.m$summary.fitted.values[fitIndex,]

Pi1=fitted$mean

# Pi1
# 
# observed
# 

D1  <- (Cs.Copies - Pi1) / sqrt(Pi1)

summary(D1)

# e <- Cs.Copies-Pi1
# 
# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1

MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(utm$Xkm),
                     Y = as.numeric(utm$Ykm))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 5 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 0.05,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(Km)") + ylab("Sample variogram")
p <- p + theme(text = element_text(size=15))+ylim(0,225)
p 

p1

1-110/150

# 
# 
# # Step 5.	Make a stack.
# 
# # 
# r2_train_spde_a=r2_val_spde_a=rmse_train_spde_a=rmse_val_spde_a=c() # Prepare empty array to store goodness-of-fit statistics
# # 
# waic_index1=matrix(nrow=10,ncol=100)
# waic_index2=matrix(nrow=10,ncol=100)
# waic_index3=matrix(nrow=10,ncol=100)
# waic_index4=matrix(nrow=10,ncol=100)
# waic_index5=matrix(nrow=10,ncol=100)
# waic_index6=matrix(nrow=10,ncol=100)
# waic_index7=matrix(nrow=10,ncol=100)
# waic_index8=matrix(nrow=10,ncol=100)
# waic_index9=matrix(nrow=10,ncol=100)
# waic_index10=matrix(nrow=10,ncol=100)
# waic_index11=matrix(nrow=10,ncol=100)
# waic_index12=matrix(nrow=10,ncol=100)
# waic_index13=matrix(nrow=10,ncol=100)
# waic_index14=matrix(nrow=10,ncol=100)
# waic_index15=matrix(nrow=10,ncol=100)
# waic_index16=matrix(nrow=10,ncol=100)
# 
# waic_array_Cs.inla.d2.d.dm=c()
# waic_array_Cs.inla.d2.d.dt=c()
# waic_array_Cs.inla.d2.d.v=c()
# waic_array_Cs.inla.d2.d=c()
# waic_array_Cs.inla.d=c()
# waic_array_Cs.inla.d.v=c()
# waic_array_Cs.inla.d.dt=c()
# waic_array_Cs.inla.d.dm=c()
# waic_array_Cs.inla.dm=c()
# waic_array_Cs.inla.dt=c()
# waic_array_Cs.inla.v=c()
# waic_array_Cs.inla.v.t=c()
# waic_array_Cs.inla.dm.t=c()
# waic_array_Cs.inla.dt.t=c()
# waic_array_Cs.inla.t=c()
# waic_array_Cs.inla.null=c()
# 
# # # 
# # 
# # all.form= Cs.Copies ~ D.M..z + DistanceFromShore.m.z.z + Type01 + V.M..z +D.thal.m.z +
# #   f(Site,model="iid")
# # set.seed(NULL)
# # 
# 
# RMSE=function(set,outcome,data,fit){
#   res = data[set,outcome]-fit[set]
#   RMSE_val <- sqrt(mean(res^2,na.rm=T)) 
#   return(RMSE_val)  
# }
# 
# pseudo_r2=function(set,outcome,data,fit){
#   res =  data[set,outcome]-fit[set]
#   RRes=sum((res)^2,na.rm = T)
#   RRtot=sum((data[set,outcome]-mean(fit[set],na.rm=T))^2,na.rm = T)
#   pseudo_r2_val=1-RRes/RRtot
#   return(pseudo_r2_val)  
# }
# 
# 
# 
# 
# library(foreach)
# library(doParallel)

#setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-6) #not to overload your computer
#registerDoParallel(cl)
# 
#
# tic()



#Make the Model Formulas
#1

sf.d2.d.dm<-Copies~-1+Depth2+D.M..z+DistanceFromShore.m.z.z+ f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)

#2
sf.d2.d.dt<- Copies~ -1+Depth2+D.M..z+D.thal.m.z+ f(w, model = spde)

#3
sf.d2.d.v<- Copies~ -1 + Depth2+D.M..z+V.M..z+ f(w, model = spde)

#4
sf.d2.d<- Copies~ -1+ Depth2+D.M..z+f(w, model = spde)

#5
sf.d<-  Copies~ -1+ D.M..z+f(w, model = spde)

#6
sf.d.v<- Copies~ -1 + D.M..z+V.M..z+ f(w, model = spde)

#7
sf.d.dt<- Copies~ -1 + D.M..z+D.thal.m.z+ f(w, model = spde)

#8
sf.d.dm<- Copies~ -1 + DistanceFromShore.m.z.z+D.M..z +f(w, model = spde)

#9
sf.dm<- Copies~ -1 + DistanceFromShore.m.z.z + f(w, model = spde)

#10
sf.dt<- Copies~ -1 + D.thal.m.z +f(w, model = spde)

#11
sf.v<-Copies~ -1 + V.M..z+ f(w, model = spde)

#12
sf.v.t<-Copies~ -1 + V.M..z+ Type01+ f(w, model = spde)

#13
sf.dm.t<-Copies~ -1 + Type01+ DistanceFromShore.m.z.z+ f(w, model = spde)

#14
sf.dt.t<-Copies~ -1 + Type01+ D.thal.m.z+ f(w, model = spde)

#15
sf.t<-Copies~ -1 + Type01+  f(w, model = spde)

#14
sf.null<-Copies~ -1 + f(w, model = spde)

# 
# 
# for (p in 1:100){
#   set.seed(NULL)
#   spec = c(val1 = .1, val2 = .1, val3 = .1,
#            val4 = .1, val5 = .1, val6 = .1,
#            val7 = .1, val8 = .1, val9 = .1,
#            val10 = .1)
# 
#   g = sample(cut(
#     seq(nrow(CSD2)),
#     nrow(CSD2)*cumsum(c(0,spec)),
#     labels = names(spec)
#   ))
# #
# #   #
# #   # foreach(k =1:10, .combine = cbind, .multicombine = T,
# #   #          .packages =c('tidyverse','fields','sp','stats',
# #   #                       'ggplot2','lattice','INLA','plotly','geometry','viridis','tictoc','kableExtra',
# #   #                       'rgdal','gstat','remotes','inlamesh3d','inlatools','INLAutils','ggregplot')
# #   # ) %dopar% {
# #   #
# #
#   for(k in 1:10){
#     # Define the index
#     index_val=which(g==paste0("val",k))
#     index_train=which(g!=paste0("train",k))
# 
#     # Define the weights
#     A.train <- inla.spde.make.A(mesh=mesh2d, loc=Loc[index_train,])
#     A.val <- inla.spde.make.A(mesh=mesh2d, loc=Loc[index_val,])
# #     #
# #
# #
# #
# #
# #     # Build the stack
# #
# #
# #
#     stackform<-as.formula(all.form)
# #
#     modrando<-"Site"
#     terms <- all.form[3]
#     terms <- gsub(" ", "", as.character(terms))
#     terms <- unlist(strsplit(terms, '\\+'))
#     terms <- terms[-grep('iid', terms)]
#     if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
#     terms <- c(terms, modrando)
# #
# #
#     ## make the formula for the stack
#     stackform <- formula(paste('~', paste(terms, collapse='+')))
# 
#     # put together a dataframe of values for the for validation stack
#     Xm <- model.matrix(stackform, data = CSD2)
#     Xm <- as.data.frame(Xm)
#     Xm.val <- Xm[index_val,]
#     N <- length(index_val)
#     #nrow(CSD2)
#     
#     stack.val <- inla.stack(
#       remove.unused=TRUE,
#       tag = "val",
#       data = list(Copies=Cs.Copies[index_val]),
#       A = list(1, 1, A.val),
#       effects = list(
#         Intercept = rep(1, N),
#         Xm.val        = Xm.val[,-1],    #Covariates without the intercept
#         w         = w.index))
# 
#     # put together a dataframe of values for the stack
#     Xm.train <- Xm[index_train,]
#     N <- #nrow(CSD2)
#       length(index_train)
# 
#     stack.train <- inla.stack(
#       remove.unused=TRUE,
#       tag = "train",
#       data = list(Copies=Cs.Copies[index_train]),
#       A = list(1, 1, A.train),
#       effects = list(
#         Intercept = rep(1, N),
#         Xm.train        = Xm.train[,-1],    #Covariates without the intercept
#         w         = w.index))
# 
#     join.stack <- inla.stack(stack.train,stack.val)
#     #, stack.val)
# 
#      # Write the test formula
# #
# 
# 
# 
# 
# 
# 
#     # Fit the models
#     #1
# 
#     Chin.inla.d2.d.dm <- inla(formula = update(sf.d2.d.dm, . ~ . + f(Site, model ="iid")),
#                          family = "nbinomial",
#                          data = inla.stack.data(join.stack),
#                          control.compute = list( waic = TRUE, config=TRUE),
#                          control.predictor = list(link = 1, compute = TRUE,
#                                                   A = inla.stack.A(join.stack)))
#     #2
#     Chin.inla.d2.d.dt <- inla(formula = update(sf.d2.d.dt, . ~ . + f(Site, model ="iid")),
#                           family = "nbinomial",
#                           data = inla.stack.data(join.stack),
#                           control.compute = list( waic = TRUE,config=TRUE),
#                           control.predictor = list(link = 1, compute = TRUE,
#                                                    A = inla.stack.A(join.stack)))
#     #3
#     Chin.inla.d2.d.v <- inla(formula = update(sf.d2.d.v, . ~ . + f(Site, model ="iid")),
#                           family = "nbinomial",
#                           data = inla.stack.data(join.stack),
#                           control.compute = list(waic = TRUE,config=TRUE),
#                           control.predictor = list(link = 1, compute = TRUE,
#                                                    A = inla.stack.A(join.stack)))
#     #4
#     Chin.inla.d2.d <- inla(formula = update(sf.d2.d, . ~ . + f(Site, model ="iid")),
#                           family = "nbinomial",
#                           data = inla.stack.data(join.stack),
#                           control.compute = list(waic = TRUE,config=TRUE),
#                           control.predictor = list(link = 1, compute = TRUE,
#                                                    A = inla.stack.A(join.stack)))
#     #5
#     Chin.inla.d <- inla(formula = update(sf.d, . ~ . + f(Site, model ="iid")),
#                              family = "nbinomial",
#                              data = inla.stack.data(join.stack),
#                              control.compute = list( waic = TRUE,config=TRUE),
#                              control.predictor = list(link = 1, compute = TRUE,
#                                                       A = inla.stack.A(join.stack)))
#     #6
#     Chin.inla.d.v <- inla(formula = update(sf.d.v, . ~ . + f(Site, model ="iid")),
#                             family = "nbinomial",
#                             data = inla.stack.data(join.stack),
#                             control.compute = list(waic = TRUE, config=TRUE),
#                             control.predictor = list(link = 1, compute = TRUE,
#                                                      A = inla.stack.A(join.stack)))
# 
#     #7
#     Chin.inla.d.dt <- inla(formula = update(sf.d.dt, . ~ . + f(Site, model ="iid")),
#                                family = "nbinomial",
#                                data = inla.stack.data(join.stack),
#                                control.compute = list(waic = TRUE,config=TRUE),
#                                control.predictor = list(link = 1, compute = TRUE,
#                                                         A = inla.stack.A(join.stack)))
#     #8
#     Chin.inla.d.dm <- inla(formula = update(sf.d.dm, . ~ . + f(Site, model ="iid")),
#                             family = "nbinomial",
#                             data = inla.stack.data(join.stack),
#                             control.compute = list(waic = TRUE, config=TRUE),
#                             control.predictor = list(link = 1, compute = TRUE,
#                                                      A = inla.stack.A(join.stack)))
# 
#     #9
# 
#     Chin.inla.dm<- inla(formula = update(sf.dm, . ~ . + f(Site, model ="iid")),
#                            family = "nbinomial",
#                            data = inla.stack.data(join.stack),
#                            control.compute = list(waic = TRUE, config=TRUE),
#                            control.predictor = list(link = 1, compute = TRUE,
#                                                     A = inla.stack.A(join.stack)))
# 
#     #10
#     Chin.inla.dt <- inla(formula = update(sf.dt, . ~ . + f(Site, model ="iid")),
#                                family = "nbinomial",
#                                data = inla.stack.data(join.stack),
#                                control.compute = list(waic = TRUE,config=TRUE),
#                                control.predictor = list(link = 1, compute = TRUE,
#                                                         A = inla.stack.A(join.stack)))
#     #11
#     Chin.inla.v <- inla(formula = update(sf.v, . ~ . + f(Site, model ="iid")),
#                              family = "nbinomial",
#                              data = inla.stack.data(join.stack),
#                              control.compute = list(waic = TRUE, config=TRUE),
#                              control.predictor = list(link = 1, compute = TRUE,
#                                                       A = inla.stack.A(join.stack)))
# 
#     #12
#     Chin.inla.v.t<- inla(formula = update(sf.v.t, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE, config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
#     #13
#     Chin.inla.dm.t<- inla(formula = update(sf.dm.t, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE, config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
# 
# 
#     #14
#     Chin.inla.dt.t<- inla(formula = update(sf.dt.t, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE, config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
#     #15
#     Chin.inla.t<- inla(formula = update(sf.t, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE, config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
#     #16
#     Chin.inla.null<- inla(formula = update(sf.null, . ~ . + f(Site, model ="iid")),
#                           family = "nbinomial",
#                           data = inla.stack.data(join.stack),
#                           control.compute = list(waic = TRUE, config=TRUE),
#                           control.predictor = list(link = 1, compute = TRUE,
#                                                    A = inla.stack.A(join.stack)))
# 
#     
# 
# end_time=Sys.time()
# #
#     # Extract the fitted values
#     index_inla_train = inla.stack.index(join.stack,"train")$data
#     index_inla_val = inla.stack.index(join.stack,"val")$data
# 
#     model.list=list(
#       Chin.inla.d2.d.dm,
#       Chin.inla.d2.d.dt,
#       Chin.inla.d2.d.v,
#       Chin.inla.d2.d,
#       Chin.inla.d,
#       Chin.inla.d.v,#6
#       Chin.inla.d.dt,
#       Chin.inla.d.dm,
#       Chin.inla.dm,#9
#       Chin.inla.dt,
#       Chin.inla.v,
#       Chin.inla.v.t,
#       Chin.inla.dm.t,
#       Chin.inla.dt.t,
#       Chin.inla.t,
#       Chin.inla.null
#     )
# 
# #
# #
#     r2_train_spde=c()
#     r2_val_spde=c()
#     rmse_train_spde=c()
#     rmse_val_spde=c()
# #
#     for (i in 1:length(model.list)){
# #
# #
#       model=model.list[[i]]
# 
# 
#       results.train=model$summary.fitted.values$mean[index_inla_train]
#       results.val=model$summary.fitted.values$mean[index_inla_val]
# #
#       M_fit_spde=array(NA,length(model$summary.fitted.values[,"mean"]))
#       M_fit_spde[index_train]=results.train
#       M_fit_spde[index_val]=results.val
# 
#       # Compute goodness-of-fit statistics
#       r2_train_spde[i]=pseudo_r2(index_train,"Chinook.vol.cor.conc",CSD2,M_fit_spde)
#       r2_val_spde[i]=pseudo_r2(index_val,"Chinook.vol.cor.conc",CSD2,M_fit_spde)
#       rmse_train_spde[i]=RMSE(index_train,"Chinook.vol.cor.conc",CSD2,M_fit_spde)
#       rmse_val_spde[i]=RMSE(index_val,"Chinook.vol.cor.conc",CSD2,M_fit_spde)
# 
# #
# #
# #
# #
#       waic_index1[k,p]=waic_array_Chin.inla.d2.d.dm=r2_val_spde[1]
#       waic_index2[k,p]=waic_array_Chin.inla.d2.d.dt=r2_val_spde[2]
#       waic_index3[k,p]=waic_array_Chin.inla.d2.d.v=r2_val_spde[3]
#       waic_index4[k,p]= waic_array_Chin.inla.d2.d=r2_val_spde[4]
#       waic_index5[k,p]=waic_array_Chin.inla.d=r2_val_spde[5]
#       waic_index6[k,p]=waic_array_Chin.inla.d.v=r2_val_spde[6]
#       waic_index7[k,p]=waic_array_Chin.inla.d.dt=r2_val_spde[7]
#       waic_index8[k,p]=waic_array_Chin.inla.d.dm=r2_val_spde[8]
#       waic_index9[k,p]=waic_array_Chin.inla.dm=r2_val_spde[9]
#       waic_index10[k,p]=waic_array_Chin.inla.dt=r2_val_spde[10]
#       waic_index11[k,p]=waic_array_Chin.inla.v=r2_val_spde[11]
#       waic_index12[k,p]=waic_array_Chin.inla.v.t=r2_val_spde[12]
#       waic_index13[k,p]=waic_array_Chin.inla.dm.t=r2_val_spde[13]
#       waic_index14[k,p]=waic_array_Chin.inla.dt.t=r2_val_spde[14]
#       waic_index15[k,p]=waic_array_Chin.inla.t=r2_val_spde[15]
#       waic_index16[k,p]=waic_array_Chin.inla.null=r2_val_spde[16]
#    
# #
#     }
# 
#   }
#   #print(paste("Average R2 in the validation set:", round(mean(r2_val_spde_a)*100),"%",
#   #            "\nAverage R2 in the training set:", round(mean(r2_train_spde_a)*100),"%"))
#   print(paste("another one bites the dust"))
# 
# }
# toc()
#
#
# Chinook.metawaic=list(waic_index1,waic_index2,waic_index3,waic_index4,waic_index5,waic_index6,
#                       waic_index7,waic_index8,waic_index9,waic_index10,waic_index11,waic_index12,
#                       waic_index13,waic_index14,waic_index15,waic_index16)

#save(Chinook.metawaic,file="Chinook.metar.new.rda")


#### END OF KFCV ####
# 
# summary(Chin.inla.null)
# summary(Chin.inla.dm.t)
# summary(Chin.inla.dt.t)
# summary(Chin.inla.t)
# summary(Chin.inla.d)




load("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis/Chinook.metar.new.rda")
msr=matrix(nrow=10,ncol=100)
msr2=matrix(nrow=10,ncol=100)

the.data=Chinook.metawaic
# 
# mean(unlist(the.data[16]))
# mean(unlist(the.data[21])-unlist(the.data[22]))
# unlist(the.data[21])
# unlist(the.data[22])
# 
# mean(unlist(the.data[21])-unlist(the.data[11]))

#compare each position
#for(k in 1:100){
  #for(i in 1:10){
for(k in 1:100){
for(i in 1:10){
    vector= as.numeric((lapply(the.data, function(x) return(x[i,k]))))
    value=max(vector)
    msr[i,k]=  ifelse(value==vector[1],1,
                      ifelse(value==vector[2],2,
                             ifelse(value==vector[3],3,
                                    ifelse(value==vector[4],4,
                                           ifelse(value==vector[5],5,
                                                  ifelse(value==vector[6],6,
                                                         ifelse(value==vector[7],7,
                                                                ifelse(value==vector[8],8,
                                                                       ifelse(value==vector[9],9,
                                                                              ifelse(value==vector[10],10,
                                                                                     ifelse(value==vector[11],11, 
                                                                                            
                                                                                            ifelse(value==vector[12],12,
                                                                                                   ifelse(value==vector[13],13,
                                                                                                          ifelse(value==vector[14],14,
                                                                                                                 ifelse(value==vector[15],15,
                                                                                                                        ifelse(value==vector[16],16,
                                                                                                                               ifelse(value==vector[17],17,
                                                                                                                                      ifelse(value==vector[18],18,
                                                                                                                                             ifelse(value==vector[19],19,
                                                                                                                                                    ifelse(value==vector[20],20,
                                                                                                                                                           ifelse(value==vector[21],21,
                                                                                                                                                                  ifelse(value==vector[22],22
                                                                                                                                                                         
                                                                                                                                                                  ))))))))))))))))))))))                                                                                                  
    
    
    
    
    
  }
}












msr3=matrix(nrow=10,ncol=100)

#compare each position
for(k in 1:100){
  for(i in 1:10){
    vector= as.numeric((lapply(the.data, function(x) return(x[i,k]))))
    value=max(vector)                                                                                                    
    msr3[i,k]=  ifelse(value==vector[1],value,
                       ifelse(value==vector[2],value,
                              ifelse(value==vector[3],value,
                                     ifelse(value==vector[4],value,
                                            ifelse(value==vector[5],value,
                                                   ifelse(value==vector[6],value,
                                                          ifelse(value==vector[7],value,
                                                                 ifelse(value==vector[8],value,
                                                                        ifelse(value==vector[9],value,
                                                                               ifelse(value==vector[10],value,
                                                                                      ifelse(value==vector[11],value, 
                                                                                             
                                                                                             ifelse(value==vector[12],value,
                                                                                                    ifelse(value==vector[13],value,
                                                                                                           ifelse(value==vector[14],value,
                                                                                                                  ifelse(value==vector[15],value,
                                                                                                                         ifelse(value==vector[16],value,
                                                                                                                                ifelse(value==vector[17],value,
                                                                                                                                       ifelse(value==vector[18],value,
                                                                                                                                              ifelse(value==vector[19],value,
                                                                                                                                                     ifelse(value==vector[20],value,
                                                                                                                                                            ifelse(value==vector[21],value,
                                                                                                                                                                   ifelse(value==vector[22],value
                                                                                                                                                                   )
                                                                                                                                                            )
                                                                                                                                                     )
                                                                                                                                              )
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                       )
                                                                                                                                )
                                                                                                                                
                                                                                                                         ))))))))))))))))
  }
}

mean(msr3)





msr2=as.factor(msr)
msr
soomeri=order(summary(msr2),decreasing = F)
soomeri
soomeri2=summary(msr2)


soomeri2=as.data.frame(soomeri2)
soomeri2

soomeri2.mod=c(1:16)

#soomeri2=c(3,3,35,8,18,7,22,10,25,0,6,84,75,17,64,57,18,154,71,97,219,7)

soomeri2=c(43,59,7,41,25,12,42,38,77,46,48,38,238,161,83,42)


soomeri2=data.frame(soomeri2,soomeri2.mod)

sort(soomeri2$soomeri2, decreasing = T)

# 
# Model.Num=soomeri2$mod
# times.selected=sort(soomeri2$soomeri2, decreasing = T)
soomeri2$weight=soomeri2$soomeri2/1000*100
transferability=c()

for (i in soomeri2$soomeri2.mod){
  transferability[i]= mean(unlist(the.data[i]))
  
}
#transferability=c(transferability[1:9],transferability[11:22])

mod.sel=data.frame(soomeri2,transferability)
#mod.sel.sort=sort(mod.sel$weight,decreasing=T)

mod.sel
#mod.sel.sort
# 
# model.list=list(nospat.inla.m,nospat.inla.d,nospat.inla.v,nospat.inla.t,
#                 nospat.inla.thal,nospat.inla.d.m,nospat.inla.thal.d,
#                 nospat.inla.t.m,nospat.inla.d.v,nospat.inla.thal.d,
#                 nospat.inla.null,spat.inla.d,spat.inla.m,spat.inla.t,
#                 spat.inla.v,spat.inla.thal,spat.inla.d.v,spat.inla.d.m,
#                 spat.inla.t.v,spat.inla.t.m,spat.inla.thal.d,spat.inla.null)
# mod.rsq=list()


model.name=c("Depth-Squared + Depth + Distance From Shore",
             "Depth-Squared + Depth + Distance From Thalweg",
             "Depth-Squared + Depth + Velocity",
             "Depth-Squared + Depth",
             "Depth",
             "Depth + Velocity",
             "Depth + Distance From Thalweg",
             "Depth + Distance From Shore",
             "Distance From Shore",
             "Distance From Thalweg",
             "Velocity",
             "Velocity + Sample Type",
             "Distance From Shore + Sample Type",
             "Distance From Thalweg + Sample Type",
             "Sample Type",
             "NULL"
            
             
             
)



mod.sel.sort=data.frame(mod.sel,model.name, row.names = NULL)
mod.sel.sort


#colnames(mod.sel.sort)=c("weight","model.name")

#mod.sel.sort=left_join(mod.sel.sort,mod.sel,by="weight")


####R2####

all.form2=Cs.Copies ~ D.M..z + Type01 + DistanceFromShore.m.z.z+ f(Site, model = "iid")

stackform<-as.formula(all.form2)


modrando<-"Site"
terms <- all.form2[3]
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

# put together a dataframe of values for the stack
Xm <- model.matrix(stackform, data = CSD2)
Xm <- as.data.frame(Xm)
N <- nrow(CSD2)



StackFit2 <- inla.stack(
  remove.unused=TRUE,
  tag = "Fit",
  data = list(Copies=Cs.Copies), 
  A = list(1, 1, A),                  
  effects = list(   
    Intercept = rep(1, N),
    Xm        = Xm[,-1],    #Covariates without the intercept
    w         = w.index))

#5
Chin.inla.d <- inla(formula = update(sf.d, . ~ . + f(Site, model ="iid")),
                    family = "nbinomial",
                    data = inla.stack.data(StackFit2),
                    control.compute = list( waic = TRUE,config=TRUE),
                    control.predictor = list(link = 1, compute = TRUE,
                                             A = inla.stack.A(StackFit2)))

#15
Chin.inla.t<- inla(formula = update(sf.t, . ~ . + f(Site, model ="iid")),
                   family = "nbinomial",
                   data = inla.stack.data(StackFit2),
                   control.compute = list(waic = TRUE, config=TRUE),
                   control.predictor = list(link = 1, compute = TRUE,
                                            A = inla.stack.A(StackFit2)))


#16
Chin.inla.null<- inla(formula = update(sf.null, . ~ . + f(Site, model ="iid")),
                      family = "nbinomial",
                      data = inla.stack.data(StackFit2),
                      control.compute = list(waic = TRUE, config=TRUE),
                      control.predictor = list(link = 1, compute = TRUE,
                                               A = inla.stack.A(StackFit2)))

Chin.inla.top<- inla(formula = update(sf.dm.t, . ~ . + f(Site, model ="iid")),
                      family = "nbinomial",
                      data = inla.stack.data(StackFit2),
                      control.compute = list(waic = TRUE, config=TRUE),
                      control.predictor = list(link = 1, compute = TRUE,
                                               A = inla.stack.A(StackFit2)))


# 

#empirical cumulative distribution

# dc.cs=distribution_check(Cs.inla.null)
# pl
hist(Cs.Copies,breaks=100)

summary(Chin.inla.t)
#summary(Chin.inla.t)
mod.sel.sort

colnames(mod.sel.sort)=c("Times Selected","Model Number","Weight","Transferability","Model Covariates")
#mod.sel.sort$mod=as.character((mod.sel.sort$mod))

#mod.sel.sort=left_join(mod.sel.sort,mod.rsq,by="mod")
mod.sel.sort=mod.sel.sort%>%mutate_if(is.numeric, format, digits=2,nsmall = 0)

mod.sel.sort



mod.sel.sorted=data.frame(mod.sel.sort[,5],mod.sel.sort[,3],mod.sel.sort[,4])

mod.sel.sorted
mod.sel.sorted <- mod.sel.sorted%>% arrange(desc(mod.sel.sort...4.))

mod.sel.sorted
n.covariates=c(2,3,3,2,2,2,2,1,3,2,1,2,1,1,1,0)




mod.sel.sorted=data.frame(mod.sel.sorted,n.covariates)
mod.sel.sorted
mod.sel.sorted$mod.sel.sort...4.=as.numeric(mod.sel.sorted$mod.sel.sort...4.)
mod.sel.sorted$n.covariates=as.numeric(mod.sel.sorted$n.covariates)

#ttt=expression("Transferability"~r^{"2"})



p_text <- function(x,y) {
  case_when(
    x >= max(mod.sel.sorted$mod.sel.sort...4.)- 0.05 & y==1|0 ~ "*"
    
  )
}

stars=p_text(mod.sel.sorted$mod.sel.sort...4.,mod.sel.sorted$n.covariates)
stars=str_replace_na(stars,"")


mod.sel.sorted$mod.sel.sort...4.=str_c(mod.sel.sorted$mod.sel.sort...4.,stars)

mod.sel.sorted <- mod.sel.sorted%>% arrange(desc(mod.sel.sort...4.))

mod.sel.sorted=mod.sel.sorted[,1:3]

kbl(mod.sel.sorted, caption = "<b>Table 2: <i>O. tshawytscha</i> models 
    with model weights (percent of times selected) and cross-validation transferability (mean k-fold R<sup>2</sup>).<b/>",
  col.names = c("Model", "Weight(%)", "Transferability"))%>%  
  kable_classic(full_width = F, html_font = "Cambria")%>%
  add_footnote("Indicates most parsimonious models with transferability within 0.05 of of top model.", 
               notation="symbol",threeparttable = T, escape=T) 





disc.siddm=dispersion_check(Chin.inla.d)
plot(disc.siddm)

ds.siddm=distribution_check(Chin.inla.d)
plot(ds.siddm)





summary(Chin.inla.t)
summary(Chin.inla.d)




SpatField.w <- inla.spde2.result(inla = Chin.inla.d,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)
Kappa <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x),
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.range.nominal[[1]] )

Sigma_u
Range

result = inla.spde2.result(Chin.inla.d, "w", spde)

par(mar=c(6,6,6,6))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density",xlab="Nominal Range", ylab="Frequency")

plot(result[["marginals.kappa"]][[1]], type = "l",
     main = "Kappa, posterior density", xlab="Kappa", ylab="Frequency")
Kappa


LocMesh <- mesh2d$loc[,1:2]
D<-dist(LocMesh)
D <- as.matrix(D)
#D
# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, #max(D)
             .04, length = 100)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1

Cor.M

#.338, .076, .003
cor.plot.data=data.frame(d.vec,Cor.M)
Range
cor.plot.data

cor.plot=ggplot()+geom_line(data=cor.plot.data,aes(x=d.vec*1000,y=Cor.M))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(0,50)+
  ylim(0,1)+
 geom_abline(aes(intercept = 0, slope = 0),linetype=3)+
  xlab("Distance (m)")+
  ylab("Matern Correlation Values")+
  theme(axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))#+
 # geom_text(x = 9, y = 0.08, aes(label = "10% modeled correlation"))

cor.plot
Kappa


fitIndex <- inla.stack.index(StackFit2, tag='Fit')$data
fitted.b<- Chin.inla.d$summary.fitted.values[fitIndex,]
fitted.nonspat<-Chin.inla.null$summary.fitted.values[fitIndex,]
Pred   <- fitted.b[, "mean"]
VarY <- Pred
resid   <- (Cs.Copies - Pred) / sqrt(VarY)
MyData3 <- data.frame(resid = resid,
                      Xkm = utm$Xkm,
                      Ykm = utm$Ykm)


Chin.inla.d$summary.hyperpar


INLAutils::plot_inla_residuals(Chin.inla.d, observed)

autoplot(Chin.inla.d)

fitIndex <- inla.stack.index(StackFit2, tag='Fit')$data
fitted <- Chin.inla.d$summary.fitted.values[fitIndex,]

ggplot_inla_residuals(Chin.inla.d,CSD2$Chinook.vol.cor.conc,CI = TRUE,binwidth = NULL)
ggplot_inla_residuals2(Chin.inla.d,CSD2$Chinook.vol.cor.conc, CI=TRUE,method = NA)
MyData3$MySize <- 2 * abs(MyData3$resid) / max(MyData3$resid)
MyData3$MyCol <- ifelse(MyData3$resid> 0, 1, 2)

#View(MyData3)
lattice::xyplot(MyData3$Ykm ~ MyData3$Xkm,
                data = MyData3,
                cex = MyData3$MySize,
                col = MyData3$MyCol,
                pch = 1)
par(mfrow=c(1,1))

hist(MyData3$resid, breaks = 20)

sd(CSD2$D.M.)
1-exp(-0.267)
1-exp(-0.084)

summary(Chin.inla.d)
#After accounting for the random effects of Site and Space 
#for every 48cm increase in depth, we can expect 
#between a 8% and a 23.4% decrease in the number of copies of Chinook DNA 
#per liter in the Klamath River.

summary(Chin.inla.t)
1-exp(-0.499)
1-exp(-0.153)

#After accounting for the random effects of Site and Space,
# we expect between a 14.2 and a 39.3% decrease in the number of copies of Chinook
#DNA per liter in the Klamath River when a sample is collected with a isokinetic sampler
# rather than a whirlpak bag.


#rm(mdl)

library(stats)


summary(Chin.inla.d)

source("INLA_plotting_functions.R")



# GRSpred <- inlaPredict(mdl=Chin.inla.d,mdl.stack = StackFit2, spatialMod=T,
#                        orig.dat=CSD2,data.pts=45)
# 
predPlot4 <- inlaPlot(inlaPrediction = GRSpred, plotVar = 'D.M..z', orig.dat = CSD2,
                      unscale=F, scale.par=scale.par,
                      xlabl='Scaled Depth', ylabl='Multiplicative Effect',printplotobject = T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),axis.text.x = element_text(size = 16))+
  theme(axis.title = element_text(face="bold"))

predPlot4



predPlot5<- inlaPlot(inlaPrediction = GRSpred, plotVar = 'DistanceFromShore.m.z.z', orig.dat = CSD2,
                     unscale=F, scale.par=scale.par,
                     xlabl='Scaled Distance to Shore', ylabl='Multiplicative Effect',printplotobject = T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),axis.text.x = element_text(size = 16))+
  theme(axis.title = element_text(face="bold"))

predPlot5



sd(CSD2$Distance.From.Shore)
summary(Chin.inla.top)


1-exp(-0.170)
exp(0.022)-1



#Not selected in best model, but for every 18 meter increase in distance to shore
#we expect between a 15.6 % decrease to a 2.2% increase in O. tshawytscha eDNA
#concentration. There is a 92.5% chance that increasing distance to shore has a
# negative effect on O. tshawytscha eDNA conc.



## chose a marginal and compare the with the results computed by the
## inla-program
r = Chin.inla.d$summary.fixed["x",]
m = Chin.inla.d$marginals.fixed$D.M..z

## compute the 95% HPD interval
inla.hpdmarginal(0.95, m)

x = seq(-6, 6, len = 1000)
y = dnorm(x)
inla.hpdmarginal(0.95, list(x=x, y=y))

## compute the the density for exp(r), version 1
r.exp = inla.tmarginal(exp, m)
## or version 2
r.exp = inla.tmarginal(function(x) exp(x), m)

## to plot the marginal, we use the inla.smarginal, which interpolates (in
## log-scale). Compare with some samples.
plot(inla.smarginal(m), type="l",xlab="marginal",ylab="probability")
abline(v=0,col="red")

s = inla.rmarginal(1000, m)
hist(inla.rmarginal(1000, m), add=TRUE, prob=TRUE)
lines(density(s), lty=2)

m1 = inla.emarginal(function(x) x^1, m)
m2 = inla.emarginal(function(x) x^2, m)
stdev = sqrt(m2 - m1^2)
q = inla.qmarginal(c(0.025,0.925), m)

## inla-program results
print(r)

## inla.marginal-results (they shouldn't be perfect!)
print(c(mean=m1, sd=stdev, "0.025quant" = q[1], "0.925quant" = q[2]))
## using the buildt-in function
inla.zmarginal(m)


#92.5% chance that distance to shore has a negative effect on Chinook eDNA conc.


fitted.d<- Chin.inla.d$summary.fitted.values[fitIndex,]
RRtot=sum((Cs.Copies-mean(Cs.Copies))^2)
RRes=sum((Cs.Copies-fitted.d$mean)^2)
pseudo_r2_val.d=1-RRes/RRtot
pseudo_r2_val.d
#fitted.d

txt="R^{2} == 0.45"

fitted <- Chin.inla.t$summary.fitted.values[fitIndex,]

RRtot=sum((Cs.Copies-mean(Cs.Copies))^2)
RRes=sum((Cs.Copies-fitted$mean)^2)
pseudo_r2_val.srf=1-RRes/RRtot
pseudo_r2_val.srf


fitted.n <-Chin.inla.null$summary.fitted.values[fitIndex,]

RRtot=sum((Cs.Copies-mean(Cs.Copies))^2)
RRes=sum((Cs.Copies-fitted.n$mean)^2)
pseudo_r2_val.null=1-RRes/RRtot
pseudo_r2_val.null


(pseudo_r2_val.srf-pseudo_r2_val.null)/pseudo_r2_val.null



(pseudo_r2_val.d-pseudo_r2_val.null)/pseudo_r2_val.null

#type explains 18.7% of the variance

#depth explains 17.2% of the variance


45-17.2

45-18.7

CSD2$fitted=fitted.d$mean
CSD2$fitted.null=fitted.n$mean

txt="R^{2} == 0.45"
fitted.d

plot=ggplot(CSD2,aes(x=log(fitted),y=log(Chinook.vol.cor.conc)))+
  geom_point()+
  theme_bw()+
  xlim(c(2.8,8.5))+
 ylim(c(2.8,8.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  geom_abline(aes(intercept = 0, slope = 1))+
  xlab("log(Fitted Copies/L)")+
  ylab("log(Observed Copies/L)")+
  theme(axis.title = element_text(face="bold"))

plot=plot+geom_text(x = 5.5, y = 8, aes(label = txt),parse=T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))

plot

library(ggpubr)
ggarrange(plot,cor.plot,predPlot4, 
          #labels=c("Goodness of Fit", "Imposed Matern Correlation Values"),
          ncol=1,nrow=3, legend = "top")+
  theme(plot.margin = margin(1,1,1,1, "cm")) 


# 
# + 
#   scale_x_discrete(breaks=c(0,100000,200000),
#                    labels=c("0e+00", "1e+05", "2e+05"))
#   #+geom_text(data=R2data,
# mapping = aes(
#  x = 50000,
#  y = 200000
#  ))
#

# ggplot(data = CSD2)+geom_point(aes(x=CSD2$fitted.null,y=CSD2$C.shasta.vol.cor.conc))+
#   facet_wrap(facets=CSD2$SiteN)+theme_classic()+geom_abline(aes(intercept = 0, slope = 1))
# 
# ggplot(data = CSD2)+geom_point(aes(x=CSD2$fitted.nonspat,y=CSD2$C.shasta.vol.cor.conc))+
#   facet_wrap(facets=CSD2$SiteN)+theme_classic()+geom_abline(aes(intercept = 0, slope = 1))




#summary(spat.inla.null)


# 
# 
# 
# par(mfrow = c(1, 1))
# plot(spat.inla.d$marginals.fix[[1]], type = "l", xlim=c(0,1),ylim=c(0,.2),
#      xlab = "Intercept", ylab = "Density")

#### end ####

rm(list = ls(all.names = TRUE))
#####################################
################C.shasta GLMM analysis #####
#analysis follows identical steps to Chinook

library(doSNOW)
library(doParallel)
library(doMPI)
library(tidyverse)
library(tictoc)
library(fields)
library(ggplot2)
library(sp)
library(lattice)
library(INLA)
library(plotly)
library(geometry)
library(viridis)
#library(devtools) #might need devtools to install some packaages
library(tictoc)
library(kableExtra)

library(rgdal)

#sometimes this happens...
#install.packages("gstat",force=TRUE)
library(gstat)
library(remotes)
library(INLAOutputs)


library(inlamesh3d)
library(inlatools)
library(INLAutils)
library(ggregplot)

#getwd()
###############C. shasta GLMM analysis ####





## my working directory
setwd("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis")
# source useful scripts
source("HighstatLibV13.R")
source("INLA_plotting_functions.R")


## Read in all the data
CSD2=read.csv("CSD.csv")


## Find the upper lower and mean concentrations per liter
CSD2$C.shasta.vol.cor.conc=replace_na(round(CSD2$C.shasta.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc=replace_na(round(CSD2$Chinook.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.upper=replace_na(round(CSD2$C.shasta.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.upper=replace_na(round(CSD2$Chinook.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.lower=replace_na(round(CSD2$C.shasta.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.lower=replace_na(round(CSD2$Chinook.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


## look at the data
par(mar=c(2,2,2,2))
hist(CSD2$C.shasta.vol.cor.conc,breaks=50)
hist(CSD2$Chinook.vol.cor.conc,breaks=50)

## check for multicolinearity
MyVar <- c("D.thal.m","DistanceFromShore.m.z","D.M.")

Mypairs(CSD2[,MyVar])

# procede when finding none


#scale the covariates

CSD2.scale <- dplyr::select(CSD2, all_of(MyVar)) 
names(CSD2.scale) <- paste(MyVar, '.z', sep='')

scale.par = do.call(data.frame, 
                    list(mean = apply(CSD2.scale, 2, mean),
                         sd = apply(CSD2.scale, 2, sd),
                         min = apply(CSD2.scale, 2, min),
                         max = apply(CSD2.scale, 2, max)))

CSD2.scale <- as.data.frame(scale(CSD2.scale))
CSD2 <- cbind(CSD2, CSD2.scale)


#now check all the covariates for multicolinearity

CSD2$Distance.Left.Bank.m[CSD2$Distance.Left.Bank.m==0]<-0.01
MyVar <- c("DistanceFromShore.m.z", "D.M.","V.M.","Type01","SiteN","D.thal.m")
Mypairs(CSD2[,MyVar])

#colinear variables are depth and type, and distance from shore and velocity




# plot them against the response variable
MyVar<-c("D.M..z","DistanceFromShore.m.z.z","V.M..z","Type01","SiteN","D.thal.m.z")
#standard
Cs.Copies=CSD2$C.shasta.vol.cor.conc
MyMultipanel.ggp2(Z = CSD2, 
                  varx = MyVar, 
                  vary = "C.shasta.vol.cor.conc", 
                  ylab = "Response variable",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


#make sure data is in the right format
CSD2$Site<-as.factor(CSD2$Site)
CSD2$SiteN<-as.factor(CSD2$Site)
CSD2$SiteN<-recode_factor(CSD2$SiteN,KLM=1,I5=2,TOH=3,BVR=4,KMN=5,SV=6)
Site=as.numeric(CSD2$SiteN)
CSD2=subset(CSD2,select=-Site)
D.M..z=as.numeric(CSD2$D.M..z)
DistanceFromShore.m.z.z=as.numeric(CSD2$DistanceFromShore.m.z.z)
Type01=as.factor(CSD2$Type01)
Depth2=as.numeric(CSD2$D.M..z)*as.numeric(CSD2$D.M..z)

#make base model forms
base.form <- Cs.Copies ~ D.M..z + DistanceFromShore.m.z.z +  f(Site,model="iid")
base.form.v <- Cs.Copies ~ D.M..z + V.M..z +  f(Site,model="iid")
base.form.t <-Cs.Copies ~ Type01 + DistanceFromShore.m.z.z +  f(Site,model="iid")
base.form.t.v <- Cs.Copies ~ Type01 + V.M..z +  f(Site,model="iid")
base.form.thal.d<-Cs.Copies ~ D.M..z + D.thal.m.z +  f(Site,model="iid")
base.form.thal.t<-Cs.Copies ~ Type01 + D.thal.m.z +  f(Site,model="iid")
base.form.thal<-Cs.Copies ~ D.thal.m.z +  f(Site,model="iid")

#pcprior <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))


# we'll fit a relatively full model to evaluate the proper PDF
# for the analysis

base.form.d2<-Cs.Copies~Depth2+D.M..z+DistanceFromShore.m.z.z+  f(Site,model="iid")

hist(Cs.Copies,breaks=100)

zeros=which(CSD2$C.shasta.vol.cor.conc==0)
length(zeros)

######here we fit some models without spatial effects
####normal non-spatial with distance from shore


nonspat.d2 <- inla(base.form.d2,
                   
                   control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                   family = "poisson", 
                   
                   control.predictor = list(link = 1, compute = TRUE),
                   data = CSD2)
summary(nonspat.d2)



nonspat.d2.nbn <- inla(base.form.d2,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "nbinomial", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = CSD2)
summary(nonspat.d2.nbn)

#we'll check to see if zero inflation is warranted

nonspat.d2.zip <- inla(base.form.d2,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "zeroinflatedpoisson1", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = CSD2)
summary(nonspat.d2.zip)



nonspat.d2.zinbn <- inla(base.form.d2,
                         
                         control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                         family = "zeroinflatednbinomial1", 
                         
                         control.predictor = list(link = 1, compute = TRUE),
                         data = CSD2)
summary(nonspat.d2.zinbn)


###########################Check for Overdispersion

source("Modified Distribution Check inlatools.R")
source("Modified Dispersion Check inlatools.R")

# bind_rows(
#   nonspat.d2.nbn$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2.zinbn$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = TRUE
#     ),
#   nonspat.d2.zip$summary.fixed %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = TRUE
#     )
# ) %>%
#   select(
#     distribution, zeroinflation, parameter, 
#     mean, lcl = `0.025quant`, ucl = `0.975quant`
#   ) -> summary_zip_fixed
# bind_rows(
#   nonspat.d2.nbn$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = FALSE
#     ),
#   nonspat.d2.zinbn$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "negative binomial",
#       zeroinflation = TRUE
#     ),
#   nonspat.d2.zip$summary.hyperpar %>%
#     rownames_to_column("parameter") %>%
#     mutate(
#       distribution = "Poisson",
#       zeroinflation = TRUE
#     ))
# # ) %>%
#   select(
#     distribution, zeroinflation, parameter, 
#     mean, lcl = `0.025quant`, ucl = `0.975quant`
#   ) %>%
#   arrange(parameter, distribution, zeroinflation) -> summary_zip_hyper
# 

#summary_zip_hyper
#summary_zip_fixed

# these functions check to see if the observed dispersion
# is probable given the data and the model

nonspat.d2_dc <- dispersion_check(nonspat.d2)
nonspat.d2.nbn_dc <- dispersion_check(nonspat.d2.nbn)
nonspat.d2.zip_dc <- dispersion_check(nonspat.d2.zip)
nonspat.d2.zinbn_dc <- dispersion_check(nonspat.d2.zinbn)


# Here's what we want to test our models for overdisperssion
plot(nonspat.d2_dc)

mean(nonspat.d2_dc$data)
mean(nonspat.d2_dc$model)
#the poisson model is clearly overdispersed we want the 
#  0.1 < P(D|data>D|model) < 0.9

plot(nonspat.d2.nbn_dc)

# Ok, sot the negative-binomial may be underdispersed, so we can
# check what the underdispersion parameter actually is

mean(nonspat.d2.nbn_dc$data)
mean(nonspat.d2.nbn_dc$model)

# that's low, but underdispersion can arise from an autocorrelation
# structure. Let's not rule out the NB distribution until we check
# for spatial autocorrelation



# plot(nonspat.d2.zip)
# plot(nonspat.d2.zinbn)

# check out the model
par(mar=c(3,3,3,3))

observed<-Cs.Copies

D<-INLAutils::plot_inla_residuals(nonspat.d2.nbn, observed)


#P<-autoplot(nonspat.d2.zinbn,which=(c(1:5)))
#P

#summary(nonspat.d2.zinbn)
summary(nonspat.d2.nbn)


# Let's make a simplified world to view our data in 2D to check for
# spatial-autocorrelation. Let's assume that the river is a straight
# line and that sinuosity is not important.


#average distance between sites
summary(CSD2$Meters)
(300533-212602)/5

## convert latitude and longitude to UTM coordinates in km

xy <- with(CSD2, data.frame(Sample.ID, Meters/1000/111.111, Distance.Left.Bank.m/1000/111.111))

coordinates(xy) <- c("Meters.1000.111.111", "Distance.Left.Bank.m.1000.111.111")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
utm<-NULL
utm <- data.frame(spTransform(xy, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")))
names(utm) <- c('Sample.ID', 'Xkm', 'Ykm', 'drop')

#View(utm)
CSD2<- merge(CSD2, utm)
#CSD2<-CSD2[1:211,]
par(mfrow=c(1,1))
plot(utm$Ykm~utm$Xkm)

#nonspat 
Pi1 <- nonspat.d2.nbn$summary.fitted.values[, "mean"]

# calculate the residuals
D1  <- (Cs.Copies - Pi1) / sqrt(Pi1)

summary(D1)

# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink the x coordinate by 2 orders of magnitude to decrease 
# processing time
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(utm$Xkm)/50,
                     Y = as.numeric(utm$Ykm))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 5 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3, 4, and 6 have clear spatial patterns
# in residuals, but now we need to formally evaluate the on the site scale

#View(MyData)




## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 0.06, # 60m cutoff to see autocorrelation at the site level
                        cressie = TRUE)



p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(Km)") + ylab("Semivariance")
p1 <- p + theme(text = element_text(size=15))+ylim(0,30000)
p1 +geom_hline(yintercept = 9600,color="red")+
annotate("text", x=0.04, y=9100, label="Nugget", size=7, color="red")+
  geom_hline(yintercept = 24200,color="darkgreen")+
  annotate("text", x=0.04, y=25400, label="Sill", size=7, color="darkgreen")+
  geom_vline(xintercept=0.05, color  = "brown")+
  annotate("text",x=0.052, y=10000, color="brown",size=7,angle=-90,label="Estimated Range")

#we see that spatial autocorrelation increases from 0 to 50 meters
# in our data. That's about the average width of the cross-sections

#with a nugget of 10,000 and a sill of ~25,000 we say

1-9600/24200

#~60.3% of the observed variance (at the site level) is spatial
# that's a lot.

#Now we calculate distances between our sample points and make a mesh
Loc<-NULL
Loc <- cbind(utm$Xkm/50, utm$Ykm)
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (km)",
     ylab = "Frequency")
#View(Loc)


# how far does the autocorrelation extend?

RangeGuess <-0.05 #50 meters



require(splancs)


# methods for building the mesh described in appendix A

Hull <- inla.nonconvex.hull(Loc, convex = -0.085)
MaxEdge  <- RangeGuess /6
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(0.01,10000) * MaxEdge,
                           cutoff = MaxEdge ,
                           max.n=1500)  #control max.n to be <3000

# 
mesh2d$n
# 
# Hull <- inla.nonconvex.hull(Loc, convex = -0.085)
# MaxEdge  <- RangeGuess 
# mesh2d     <- inla.mesh.2d(boundary = Hull,
#                            max.edge = c(0.001,4000) * MaxEdge,
#                            cutoff = MaxEdge / 5,
#                            max.n=1000)  #control max.n to be <3000
# 
# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1)

#note the x coordinate is artificially shrunk here to minimize model
# run time


# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


#here we use an informative prior derived from the semivariogram

spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess, 0.05), 
                            prior.sigma = c(1, 0.05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)





# make a broad model formula
all.form= Cs.Copies ~ Depth2+ D.M..z + DistanceFromShore.m.z.z + Type01 + V.M..z + 
  D.thal.m.z + f(Site, model = "iid")


#let's make the stack

stackform<-as.formula(all.form)

# random effect
modrando<-"Site"

terms <- all.form[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

# put together a dataframe of values for the stack
Xm <- model.matrix(stackform, data = CSD2)
Xm <- as.data.frame(Xm)
N <- nrow(CSD2)




StackFit <- inla.stack(
  remove.unused=TRUE,
  tag = "Fit",
  data = list(Copies=Cs.Copies), 
  A = list(1, 1, A),                  
  effects = list(   
    Intercept = rep(1, N),
    Xm        = Xm[,-1],    #Covariates without the intercept
    w         = w.index))

sf.d2.d.dm<-Copies~-1+Depth2+D.M..z+DistanceFromShore.m.z.z+ f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)

#0
spat.inla.d2.d.m<- inla(formula = update(sf.d2.d.dm, . ~ . + f(Site, model ="iid")),
                        family = "nbinomial",
                        data = inla.stack.data(StackFit),
                        control.compute = list(waic = TRUE, config=TRUE),
                        control.predictor = list(link = 1, compute = TRUE,
                                                 A = inla.stack.A(StackFit)))

summary(spat.inla.d2.d.m)


# now we check for over/under dispersion
disc.siddm=dispersion_check(spat.inla.d2.d.m)
disc.siddm$data
# There is still a low probability of the observed dispersion
# given the modelbut 0.61 is getting better
plot(disc.siddm)




dc.siddm=distribution_check(spat.inla.d2.d.m)
plot(dc.siddm)

#https://inlatools.netlify.app/articles/distribution.html here is a good
# tutorial in how these tools work


#dc2=distribution_check(spat.inla.d2.d.m)

#plot(dc2)
#mean(dc2$median)

autoplot(spat.inla.d2.d.m)

#INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)

par(mar=c(4,4,4,4))
D<-INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)


## convert latitude and longitude to UTM coordinates in km
xy <- with(CSD2, data.frame(Sample.ID, Meters/1000/111.111, Distance.Left.Bank.m/1000/111.111))

coordinates(xy) <- c("Meters.1000.111.111", "Distance.Left.Bank.m.1000.111.111")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
utm<-NULL
utm <- data.frame(spTransform(xy, CRS("+proj=utm +zone=10 +datum=WGS84 +units=km")))
names(utm) <- c('Sample.ID', 'Xkm', 'Ykm', 'drop')

#View(utm)
CSD2<- merge(CSD2, utm)
#CSD2<-CSD2[1:211,]
par(mfrow=c(1,1))
plot(utm$Ykm~utm$Xkm)


length(spat.inla.d2.d.m$summary.fitted.values[,"mean"])
#nonspat w/out velocity adj
Pi1 <- spat.inla.d2.d.m$summary.fitted.values[, "mean"]


fitIndex <- inla.stack.index(StackFit, tag='Fit')$data
fitted <- spat.inla.d2.d.m$summary.fitted.values[fitIndex,]

Pi1=fitted$mean


D1  <- (Cs.Copies - Pi1) / sqrt(Pi1)

summary(D1)

# e <- Cs.Copies-Pi1
# 
# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1

MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(utm$Xkm),
                     Y = as.numeric(utm$Ykm))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 5 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 0.06,
                        cressie = TRUE)

#hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))

p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))

p <- p + xlab("Distance(Km)") + ylab("Semivariance")
p <- p + theme(text = element_text(size=15))+ylim(0,12500)

#here is the semivariogram with the spatial effect
p 


#here is the semivariogram without the spatial effect
p1


# the spatial effect addresses the spatial structure of the data well



##########Following is the K.fold model fitting that is commented out



# Step 5.	Make a stack.

# # 
# r2_train_spde_a=r2_val_spde_a=rmse_train_spde_a=rmse_val_spde_a=c() # Prepare empty array to store goodness-of-fit statistics
# # 
# waic_index1=matrix(nrow=10,ncol=100)
# waic_index2=matrix(nrow=10,ncol=100)
# waic_index3=matrix(nrow=10,ncol=100)
# waic_index4=matrix(nrow=10,ncol=100)
# waic_index5=matrix(nrow=10,ncol=100)
# waic_index6=matrix(nrow=10,ncol=100)
# waic_index7=matrix(nrow=10,ncol=100)
# waic_index8=matrix(nrow=10,ncol=100)
# waic_index9=matrix(nrow=10,ncol=100)
# waic_index10=matrix(nrow=10,ncol=100)
# waic_index11=matrix(nrow=10,ncol=100)
# waic_index12=matrix(nrow=10,ncol=100)
# waic_index13=matrix(nrow=10,ncol=100)
# waic_index14=matrix(nrow=10,ncol=100)
# waic_index15=matrix(nrow=10,ncol=100)
# waic_index16=matrix(nrow=10,ncol=100)
# 
# waic_array_Cs.inla.d2.d.dm=c()
# waic_array_Cs.inla.d2.d.dt=c()
# waic_array_Cs.inla.d2.d.v=c()
# waic_array_Cs.inla.d2.d=c()
# waic_array_Cs.inla.d=c()
# waic_array_Cs.inla.d.v=c()
# waic_array_Cs.inla.d.dt=c()
# waic_array_Cs.inla.d.dm=c()
# waic_array_Cs.inla.dm=c()
# waic_array_Cs.inla.dt=c()
# waic_array_Cs.inla.v=c()
# waic_array_Cs.inla.v.t=c()
# waic_array_Cs.inla.dm.t=c()
# waic_array_Cs.inla.dt.t=c()
# waic_array_Cs.inla.t=c()
# waic_array_Cs.inla.null=c()
# 
# # # 
# # 
# # all.form= Cs.Copies ~ D.M..z + DistanceFromShore.m.z.z + Type01 + V.M..z +D.thal.m.z +
# #   f(Site,model="iid")
# # set.seed(NULL)
# # 
# 
# RMSE=function(set,outcome,data,fit){
#   res = data[set,outcome]-fit[set]
#   RMSE_val <- sqrt(mean(res^2,na.rm=T)) 
#   return(RMSE_val)  
# }
# 
# pseudo_r2=function(set,outcome,data,fit){
#   res =  data[set,outcome]-fit[set]
#   RRes=sum((res)^2,na.rm = T)
#   RRtot=sum((data[set,outcome]-mean(fit[set],na.rm=T))^2,na.rm = T)
#   pseudo_r2_val=1-RRes/RRtot
#   return(pseudo_r2_val)  
# }
# 
# 
# 
# 
# library(foreach)
# library(doParallel)

#setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-6) #not to overload your computer
#registerDoParallel(cl)
# 
#
# tic()



#Make the Model Formulas
#1

sf.d2.d.dm<-Copies~-1+Depth2+D.M..z+DistanceFromShore.m.z.z+ f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)

#2
sf.d2.d.dt<- Copies~ -1+Depth2+D.M..z+D.thal.m.z+ f(w, model = spde)

#3
sf.d2.d.v<- Copies~ -1 + Depth2+D.M..z+V.M..z+ f(w, model = spde)

#4
sf.d2.d<- Copies~ -1+ Depth2+D.M..z+f(w, model = spde)

#5
sf.d<-  Copies~ -1+ D.M..z+f(w, model = spde)

#6
sf.d.v<- Copies~ -1 + D.M..z+V.M..z+ f(w, model = spde)

#7
sf.d.dt<- Copies~ -1 + D.M..z+D.thal.m.z+ f(w, model = spde)

#8
sf.d.dm<- Copies~ -1 + DistanceFromShore.m.z.z+D.M..z +f(w, model = spde)

#9
sf.dm<- Copies~ -1 + DistanceFromShore.m.z.z + f(w, model = spde)

#10
sf.dt<- Copies~ -1 + D.thal.m.z +f(w, model = spde)

#11
sf.v<-Copies~ -1 + V.M..z+ f(w, model = spde)

#12
sf.v.t<-Copies~ -1 + V.M..z+ Type01+ f(w, model = spde)

#13
sf.dm.t<-Copies~ -1 + Type01+ DistanceFromShore.m.z.z+ f(w, model = spde)

#14
sf.dt.t<-Copies~ -1 + Type01+ D.thal.m.z+ f(w, model = spde)

#15
sf.t<-Copies~ -1 + Type01+  f(w, model = spde)

#14
sf.null<-Copies~ -1 + f(w, model = spde)

# 
# 
# for (p in 1:100){
#   set.seed(NULL)
#   spec = c(val1 = .1, val2 = .1, val3 = .1,
#            val4 = .1, val5 = .1, val6 = .1,
#            val7 = .1, val8 = .1, val9 = .1,
#            val10 = .1)
#   
#   g = sample(cut(
#     seq(nrow(CSD2)),
#     nrow(CSD2)*cumsum(c(0,spec)),
#     labels = names(spec)
#   ))
#   #
#   #   #
#   #   # foreach(k =1:10, .combine = cbind, .multicombine = T,
#   #   #          .packages =c('tidyverse','fields','sp','stats',
#   #   #                       'ggplot2','lattice','INLA','plotly','geometry','viridis','tictoc','kableExtra',
#   #   #                       'rgdal','gstat','remotes','inlamesh3d','inlatools','INLAutils','ggregplot')
#   #   # ) %dopar% {
#   #   #
#   #
#   for(k in 1:10){
#     # Define the index
#     index_val=which(g==paste0("val",k))
#     index_train=which(g!=paste0("train",k))
#     
#     # Define the weights
#     A.train <- inla.spde.make.A(mesh=mesh2d, loc=Loc[index_train,])
#     A.val <- inla.spde.make.A(mesh=mesh2d, loc=Loc[index_val,])
#     #     #
#     #
#     #
#     #
#     #
#     #     # Build the stack
#     #
#     #
#     #
#     stackform<-as.formula(all.form)
#     #
#     modrando<-"Site"
#     terms <- all.form[3]
#     terms <- gsub(" ", "", as.character(terms))
#     terms <- unlist(strsplit(terms, '\\+'))
#     terms <- terms[-grep('iid', terms)]
#     if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
#     terms <- c(terms, modrando)
#     #
#     #
#     ## make the formula for the stack
#     stackform <- formula(paste('~', paste(terms, collapse='+')))
#     
#     # put together a dataframe of values for the for validation stack
#     Xm <- model.matrix(stackform, data = CSD2)
#     Xm <- as.data.frame(Xm)
#     Xm.val <- Xm[index_val,]
#     N <- length(index_val)
#     #nrow(CSD2)
#     
#     stack.val <- inla.stack(
#       remove.unused=TRUE,
#       tag = "val",
#       data = list(Copies=Cs.Copies[index_val]),
#       A = list(1, 1, A.val),
#       effects = list(
#         Intercept = rep(1, N),
#         Xm.val        = Xm.val[,-1],    #Covariates without the intercept
#         w         = w.index))
#     
#     # put together a dataframe of values for the stack
#     Xm.train <- Xm[index_train,]
#     N <- #nrow(CSD2)
#       length(index_train)
#     
#     stack.train <- inla.stack(
#       remove.unused=TRUE,
#       tag = "train",
#       data = list(Copies=Cs.Copies[index_train]),
#       A = list(1, 1, A.train),
#       effects = list(
#         Intercept = rep(1, N),
#         Xm.train        = Xm.train[,-1],    #Covariates without the intercept
#         w         = w.index))
#     
#     join.stack <- inla.stack(stack.train,stack.val)
#     #, stack.val)
#     
#     # Write the test formula
#     #
#     
#     
#     
#     
#     
#     
#     # Fit the models
#     #1
#     
#     Cs.inla.d2.d.dm <- inla(formula = update(sf.d2.d.dm, . ~ . + f(Site, model ="iid")),
#                             family = "nbinomial",
#                             data = inla.stack.data(join.stack),
#                             control.compute = list( waic = TRUE, config=TRUE),
#                             control.predictor = list(link = 1, compute = TRUE,
#                                                      A = inla.stack.A(join.stack)))
#     #2
#     Cs.inla.d2.d.dt <- inla(formula = update(sf.d2.d.dt, . ~ . + f(Site, model ="iid")),
#                             family = "nbinomial",
#                             data = inla.stack.data(join.stack),
#                             control.compute = list( waic = TRUE,config=TRUE),
#                             control.predictor = list(link = 1, compute = TRUE,
#                                                      A = inla.stack.A(join.stack)))
#     #3
#     Cs.inla.d2.d.v <- inla(formula = update(sf.d2.d.v, . ~ . + f(Site, model ="iid")),
#                            family = "nbinomial",
#                            data = inla.stack.data(join.stack),
#                            control.compute = list(waic = TRUE,config=TRUE),
#                            control.predictor = list(link = 1, compute = TRUE,
#                                                     A = inla.stack.A(join.stack)))
#     #4
#     Cs.inla.d2.d <- inla(formula = update(sf.d2.d, . ~ . + f(Site, model ="iid")),
#                          family = "nbinomial",
#                          data = inla.stack.data(join.stack),
#                          control.compute = list(waic = TRUE,config=TRUE),
#                          control.predictor = list(link = 1, compute = TRUE,
#                                                   A = inla.stack.A(join.stack)))
#     #5
#     Cs.inla.d <- inla(formula = update(sf.d, . ~ . + f(Site, model ="iid")),
#                       family = "nbinomial",
#                       data = inla.stack.data(join.stack),
#                       control.compute = list( waic = TRUE,config=TRUE),
#                       control.predictor = list(link = 1, compute = TRUE,
#                                                A = inla.stack.A(join.stack)))
#     #6
#     Cs.inla.d.v <- inla(formula = update(sf.d.v, . ~ . + f(Site, model ="iid")),
#                         family = "nbinomial",
#                         data = inla.stack.data(join.stack),
#                         control.compute = list(waic = TRUE, config=TRUE),
#                         control.predictor = list(link = 1, compute = TRUE,
#                                                  A = inla.stack.A(join.stack)))
#     
#     #7
#     Cs.inla.d.dt <- inla(formula = update(sf.d.dt, . ~ . + f(Site, model ="iid")),
#                          family = "nbinomial",
#                          data = inla.stack.data(join.stack),
#                          control.compute = list(waic = TRUE,config=TRUE),
#                          control.predictor = list(link = 1, compute = TRUE,
#                                                   A = inla.stack.A(join.stack)))
#     #8
#     Cs.inla.d.dm <- inla(formula = update(sf.d.dm, . ~ . + f(Site, model ="iid")),
#                          family = "nbinomial",
#                          data = inla.stack.data(join.stack),
#                          control.compute = list(waic = TRUE, config=TRUE),
#                          control.predictor = list(link = 1, compute = TRUE,
#                                                   A = inla.stack.A(join.stack)))
#     
#     #9
#     
#     Cs.inla.dm<- inla(formula = update(sf.dm, . ~ . + f(Site, model ="iid")),
#                       family = "nbinomial",
#                       data = inla.stack.data(join.stack),
#                       control.compute = list(waic = TRUE, config=TRUE),
#                       control.predictor = list(link = 1, compute = TRUE,
#                                                A = inla.stack.A(join.stack)))
#     
#     #10
#     Cs.inla.dt <- inla(formula = update(sf.dt, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE,config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
#     #11
#     Cs.inla.v <- inla(formula = update(sf.v, . ~ . + f(Site, model ="iid")),
#                       family = "nbinomial",
#                       data = inla.stack.data(join.stack),
#                       control.compute = list(waic = TRUE, config=TRUE),
#                       control.predictor = list(link = 1, compute = TRUE,
#                                                A = inla.stack.A(join.stack)))
#     
#     #12
#     Cs.inla.v.t<- inla(formula = update(sf.v.t, . ~ . + f(Site, model ="iid")),
#                        family = "nbinomial",
#                        data = inla.stack.data(join.stack),
#                        control.compute = list(waic = TRUE, config=TRUE),
#                        control.predictor = list(link = 1, compute = TRUE,
#                                                 A = inla.stack.A(join.stack)))
#     #13
#     Cs.inla.dm.t<- inla(formula = update(sf.dm.t, . ~ . + f(Site, model ="iid")),
#                         family = "nbinomial",
#                         data = inla.stack.data(join.stack),
#                         control.compute = list(waic = TRUE, config=TRUE),
#                         control.predictor = list(link = 1, compute = TRUE,
#                                                  A = inla.stack.A(join.stack)))
#     
#     
#     #14
#     Cs.inla.dt.t<- inla(formula = update(sf.dt.t, . ~ . + f(Site, model ="iid")),
#                         family = "nbinomial",
#                         data = inla.stack.data(join.stack),
#                         control.compute = list(waic = TRUE, config=TRUE),
#                         control.predictor = list(link = 1, compute = TRUE,
#                                                  A = inla.stack.A(join.stack)))
#     #15
#     Cs.inla.t<- inla(formula = update(sf.t, . ~ . + f(Site, model ="iid")),
#                      family = "nbinomial",
#                      data = inla.stack.data(join.stack),
#                      control.compute = list(waic = TRUE, config=TRUE),
#                      control.predictor = list(link = 1, compute = TRUE,
#                                               A = inla.stack.A(join.stack)))
#     #16
#     Cs.inla.null<- inla(formula = update(sf.null, . ~ . + f(Site, model ="iid")),
#                         family = "nbinomial",
#                         data = inla.stack.data(join.stack),
#                         control.compute = list(waic = TRUE, config=TRUE),
#                         control.predictor = list(link = 1, compute = TRUE,
#                                                  A = inla.stack.A(join.stack)))
#     
#     
#     
#     end_time=Sys.time()
#     #
#     # Extract the fitted values
#     index_inla_train = inla.stack.index(join.stack,"train")$data
#     index_inla_val = inla.stack.index(join.stack,"val")$data
#     
#     model.list=list(
#       Cs.inla.d2.d.dm,
#       Cs.inla.d2.d.dt,
#       Cs.inla.d2.d.v,
#       Cs.inla.d2.d,
#       Cs.inla.d,
#       Cs.inla.d.v,#6
#       Cs.inla.d.dt,
#       Cs.inla.d.dm,
#       Cs.inla.dm,#9
#       Cs.inla.dt,
#       Cs.inla.v,
#       Cs.inla.v.t,
#       Cs.inla.dm.t,
#       Cs.inla.dt.t,
#       Cs.inla.t,
#       Cs.inla.null
#     )
#     
#     #
#     #
#     r2_train_spde=c()
#     r2_val_spde=c()
#     rmse_train_spde=c()
#     rmse_val_spde=c()
#     #
#     for (i in 1:length(model.list)){
#       #
#       #
#       model=model.list[[i]]
#       
#       
#       results.train=model$summary.fitted.values$mean[index_inla_train]
#       results.val=model$summary.fitted.values$mean[index_inla_val]
#       #
#       M_fit_spde=array(NA,length(model$summary.fitted.values[,"mean"]))
#       M_fit_spde[index_train]=results.train
#       M_fit_spde[index_val]=results.val
#       
#       # Compute goodness-of-fit statistics
#       r2_train_spde[i]=pseudo_r2(index_train,"C.shasta.vol.cor.conc",CSD2,M_fit_spde)
#       r2_val_spde[i]=pseudo_r2(index_val,"C.shasta.vol.cor.conc",CSD2,M_fit_spde)
#       rmse_train_spde[i]=RMSE(index_train,"C.shasta.vol.cor.conc",CSD2,M_fit_spde)
#       rmse_val_spde[i]=RMSE(index_val,"C.shasta.vol.cor.conc",CSD2,M_fit_spde)
#       z
#       #
#       #
#       #
#       #
#       waic_index1[k,p]=waic_array_Cs.inla.d2.d.dm=r2_val_spde[1]
#       waic_index2[k,p]=waic_array_Cs.inla.d2.d.dt=r2_val_spde[2]
#       waic_index3[k,p]=waic_array_Cs.inla.d2.d.v=r2_val_spde[3]
#       waic_index4[k,p]= waic_array_Cs.inla.d2.d=r2_val_spde[4]
#       waic_index5[k,p]=waic_array_Cs.inla.d=r2_val_spde[5]
#       waic_index6[k,p]=waic_array_Cs.inla.d.v=r2_val_spde[6]
#       wazzzic_index7[k,p]=waic_array_Cs.inla.d.dt=r2_val_spde[7]
#       waic_index8[k,p]=waic_array_Cs.inla.d.dm=r2_val_spde[8]
#       waic_index9[k,p]=waic_array_Cs.inla.dm=r2_val_spde[9]
#       waic_index10[k,p]=waic_array_Cs.inla.dt=r2_val_spde[10]
#       waic_index11[k,p]=waic_array_Cs.inla.v=r2_val_spde[11]
#       waic_index12[k,p]=waic_array_Cs.inla.v.t=r2_val_spde[12]
#       waic_index13[k,p]=waic_array_Cs.inla.dm.t=r2_val_spde[13]
#       waic_index14[k,p]=waic_array_Cs.inla.dt.t=r2_val_spde[14]
#       waic_index15[k,p]=waic_array_Cs.inla.t=r2_val_spde[15]
#       waic_index16[k,p]=waic_array_Cs.inla.null=r2_val_spde[16]
#       
#       #
#     }
#     
#   }
#   #print(paste("Average R2 in the validation set:", round(mean(r2_val_spde_a)*100),"%",
#   #            "\nAverage R2 in the training set:", round(mean(r2_train_spde_a)*100),"%"))
#   print(paste("another one bites the dust"))
#   
# }
# toc()
# #
# #
# C.shasta.metawaic=list(waic_index1,waic_index2,waic_index3,waic_index4,waic_index5,waic_index6,
#                        waic_index7,waic_index8,waic_index9,waic_index10,waic_index11,waic_index12,
#                        waic_index13,waic_index14,waic_index15,waic_index16)
# 
# save(C.shasta.metawaic,file="C.shasta.metar.new.rda")
# 


#### end ####


#we load the file from the K.fold cross-validation
# select the most parsimonious model with mean transferability within 0.05 of the top model

load("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis/C.shasta.metar.new.rda")
msr=matrix(nrow=10,ncol=100)
msr2=matrix(nrow=10,ncol=100)

the.data=C.shasta.metawaic

for(k in 1:100){
  for(i in 1:10){
    vector= as.numeric((lapply(the.data, function(x) return(x[i,k]))))
    value=max(vector)
    msr[i,k]=  ifelse(value==vector[1],1,
                      ifelse(value==vector[2],2,
                             ifelse(value==vector[3],3,
                                    ifelse(value==vector[4],4,
                                           ifelse(value==vector[5],5,
                                                  ifelse(value==vector[6],6,
                                                         ifelse(value==vector[7],7,
                                                                ifelse(value==vector[8],8,
                                                                       ifelse(value==vector[9],9,
                                                                              ifelse(value==vector[10],10,
                                                                                     ifelse(value==vector[11],11, 
                                                                                            
                                                                                            ifelse(value==vector[12],12,
                                                                                                   ifelse(value==vector[13],13,
                                                                                                          ifelse(value==vector[14],14,
                                                                                                                 ifelse(value==vector[15],15,
                                                                                                                        ifelse(value==vector[16],16,
                                                                                                                               ifelse(value==vector[17],17,
                                                                                                                                      ifelse(value==vector[18],18,
                                                                                                                                             ifelse(value==vector[19],19,
                                                                                                                                                    ifelse(value==vector[20],20,
                                                                                                                                                           ifelse(value==vector[21],21,
                                                                                                                                                                  ifelse(value==vector[22],22
                                                                                                                                                                         
                                                                                                                                                                  ))))))))))))))))))))))                                                                                                  
    
    
    
    
    
  }
}












msr3=matrix(nrow=10,ncol=100)

#compare each position
for(k in 1:100){
  for(i in 1:10){
    vector= as.numeric((lapply(the.data, function(x) return(x[i,k]))))
    value=max(vector)                                                                                                    
    msr3[i,k]=  ifelse(value==vector[1],value,
                       ifelse(value==vector[2],value,
                              ifelse(value==vector[3],value,
                                     ifelse(value==vector[4],value,
                                            ifelse(value==vector[5],value,
                                                   ifelse(value==vector[6],value,
                                                          ifelse(value==vector[7],value,
                                                                 ifelse(value==vector[8],value,
                                                                        ifelse(value==vector[9],value,
                                                                               ifelse(value==vector[10],value,
                                                                                      ifelse(value==vector[11],value, 
                                                                                             
                                                                                             ifelse(value==vector[12],value,
                                                                                                    ifelse(value==vector[13],value,
                                                                                                           ifelse(value==vector[14],value,
                                                                                                                  ifelse(value==vector[15],value,
                                                                                                                         ifelse(value==vector[16],value,
                                                                                                                                ifelse(value==vector[17],value,
                                                                                                                                       ifelse(value==vector[18],value,
                                                                                                                                              ifelse(value==vector[19],value,
                                                                                                                                                     ifelse(value==vector[20],value,
                                                                                                                                                            ifelse(value==vector[21],value,
                                                                                                                                                                   ifelse(value==vector[22],value
                                                                                                                                                                   )
                                                                                                                                                            )
                                                                                                                                                     )
                                                                                                                                              )
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                              
                                                                                                                                       )
                                                                                                                                )
                                                                                                                                
                                                                                                                         ))))))))))))))))
  }
}

msr2=as.factor(msr)
soomeri=order(summary(msr2),decreasing = F)
soomeri2=summary(msr2)
soomeri2=as.data.frame(soomeri2)
soomeri2.mod=c(1:16)
soomeri2=c(43,129,123,10,8,106,101,45,21,58,81,110,32,46,11,9)
soomeri2=data.frame(soomeri2,soomeri2.mod)
sort(soomeri2$soomeri2, decreasing = T)
soomeri2$weight=soomeri2$soomeri2/1000*100
transferability=c()

for (i in soomeri2$soomeri2.mod){
  transferability[i]= mean(unlist(the.data[i]),na.rm=T)
  
}

mod.sel=data.frame(soomeri2,transferability)



model.name=c("Depth-Squared + Depth + Distance From Shore",
             "Depth-Squared + Depth + Distance From Thalweg",
             "Depth-Squared + Depth + Velocity",
             "Depth-Squared + Depth",
             "Depth",
             "Depth + Velocity",
             "Depth + Distance From Thalweg",
             "Depth + Distance From Shore",
             "Distance From Shore",
             "Distance From Thalweg",
             "Velocity",
             "Velocity + Sample Type",
             "Distance From Shore + Sample Type",
             "Distance From Thalweg + Sample Type",
             "Sample Type",
             "NULL"
             
             
             
)


###make the table

mod.sel.sort=data.frame(mod.sel,model.name, row.names = NULL)
colnames(mod.sel.sort)=c("Times Selected","Model Number","Weight","Transferability","Model Covariates")
mod.sel.sort=mod.sel.sort%>%mutate_if(is.numeric, format, digits=2,nsmall = 0)
mod.sel.sorted=data.frame(mod.sel.sort[,5],mod.sel.sort[,3],mod.sel.sort[,4])
mod.sel.sorted <- mod.sel.sorted%>% arrange(desc(mod.sel.sort...4.))
n.covariates=c(3,3,2,2,3,2,1,2,2,2,1,1,2,1,0,1)
mod.sel.sorted=data.frame(mod.sel.sorted,n.covariates)
mod.sel.sorted$mod.sel.sort...4.=as.numeric(mod.sel.sorted$mod.sel.sort...4.)
mod.sel.sorted$n.covariates=as.numeric(mod.sel.sorted$n.covariates)

p_text <- function(x,y) {
  case_when(
    x >= max(mod.sel.sorted$mod.sel.sort...4.)- 0.05 & y==0 ~ "*"
    
  )
}

stars=p_text(mod.sel.sorted$mod.sel.sort...4.,mod.sel.sorted$n.covariates)
stars=str_replace_na(stars,"")
mod.sel.sorted$mod.sel.sort...4.=str_c(mod.sel.sorted$mod.sel.sort...4.,stars)
mod.sel.sorted <- mod.sel.sorted%>% arrange(desc(mod.sel.sort...4.))
mod.sel.sorted=mod.sel.sorted[,1:3]

kbl(mod.sel.sorted, caption = "<b>Table 1: <i>C. shasta</i> models 
    with model weights (percent of times selected) and cross-validation transferability (mean k-fold R<sup>2</sup>).<b/>",
    col.names = c("Model", "Weight(%)", "Transferability"))%>%  
  kable_classic(full_width = F, html_font = "Cambria")%>%
  add_footnote("Indicates most parsimonious models with transferability within 0.05 of of top model.", 
               notation="symbol",threeparttable = T, escape=T) 




# we select the null model




# Fit the model with all data for summary figures etc.
# same as before, but with all data points


stackform<-as.formula(all.form)


modrando<-"Site"
terms <- all.form[3]
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

# put together a dataframe of values for the stack
Xm <- model.matrix(stackform, data = CSD2)
Xm <- as.data.frame(Xm)
N <- nrow(CSD2)



StackFit <- inla.stack(
  remove.unused=TRUE,
  tag = "Fit",
  data = list(Copies=Cs.Copies), 
  A = list(1, 1, A),                  
  effects = list(   
    Intercept = rep(1, N),
    Xm        = Xm[,-1],    #Covariates without the intercept
    w         = w.index))

#5
Cs.inla.null <- inla(formula = update(sf.null, . ~ . + f(Site, model ="iid")),
                    family = "nbinomial",
                    data = inla.stack.data(StackFit),
                    control.compute = list( waic = TRUE,config=TRUE),
                    control.predictor = list(link = 1, compute = TRUE,
                                             A = inla.stack.A(StackFit)))



ds.inla.null=distribution_check(Cs.inla.null)

plot(ds.inla.null)


########Seems we still have trouble with the lowest values

SpatField.w <- inla.spde2.result(inla = Cs.inla.null,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)
Kappa <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x),
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.range.nominal[[1]] )


result = inla.spde2.result(Cs.inla.null, "w", spde)

Range

par(mar=c(4,4,4,4))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")


LocMesh <- mesh2d$loc[,1:2]
D<-dist(LocMesh)
D <- as.matrix(D)
#D
# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, #max(D)
             2, length = 10000)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1

Sigma_u
Range

cor.plot.data=data.frame(d.vec,Cor.M)
cor.plot.data

cor.plot=ggplot()+geom_line(data=cor.plot.data,aes(x=d.vec,y=Cor.M))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(0,2)+
  ylim(0,1.15)+
  geom_abline(aes(intercept = 0.05, slope = 0),linetype=3)+
  xlab("Distance (km)")+
  ylab("Matern Correlation Values")+
  theme(axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))+
  geom_text(x = 1.5, y = 0.065, aes(label = "5% autocorrelation"))

cor.plot

par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec,
     y = Cor.M,
     pch = 16,
     type = "l",
     cex.lab = 1.5,
     xlab = "Distance",
     ylab = "Correlation",
     xlim = c(0, 100))
abline(h = 0.1, lty = 2)


cor.plot.data
Kappa

fitIndex <- inla.stack.index(StackFit, tag='Fit')$data


fitted.b<- Cs.inla.null$summary.fitted.values[fitIndex,]

#fitted.nonspat<-Chin.inla.null$summary.fitted.values[fitIndex,]



Pred   <- fitted.b[, "mean"]
VarY <- Pred
resid   <- (Cs.Copies - Pred) / sqrt(VarY)

MyData3 <- data.frame(resid = resid,
                      Xkm = utm$Xkm,
                      Ykm = utm$Ykm)


INLAutils::plot_inla_residuals(Cs.inla.null, observed)

# let's look at the residuals again
ggplot_inla_residuals(Cs.inla.null,CSD2$C.shasta.vol.cor.conc,CI = TRUE,binwidth = NULL)
ggplot_inla_residuals2(Cs.inla.null,CSD2$C.shasta.vol.cor.conc, CI=TRUE,method = NA)


MyData3$MySize <- 2 * abs(MyData3$resid) / max(MyData3$resid)
MyData3$MyCol <- ifelse(MyData3$resid> 0, 1, 2)

#View(MyData3)
lattice::xyplot(MyData3$Ykm ~ MyData3$Xkm,
                data = MyData3,
                cex = MyData3$MySize,
                col = MyData3$MyCol,
                pch = 1)
par(mfrow=c(1,1))

hist(MyData3$resid, breaks = 20)

summary(Cs.inla.null)
#After accounting for the random effects of Site and Space 
#for no tested effects were included in our best model.


#fitted.d=nospat.inla.t.m$summary.fitted.values[fitIndex,]
fitted.d<- Cs.inla.null$summary.fitted.values[fitIndex,]
RRtot=sum((Cs.Copies-mean(Cs.Copies))^2)
RRes=sum((Cs.Copies-fitted.d$mean)^2)
pseudo_r2_val.d=1-RRes/RRtot
pseudo_r2_val.d


txt="R^{2} == 0.49"

CSD2=data.frame(CSD2,fitted.d)

# 49%% of the variance is accounted for by random effects


plot=ggplot(CSD2,aes(x=mean/1000,y=C.shasta.vol.cor.conc/1000))+
  geom_point()+
  theme_bw()+
xlim(c(0,200))+
  ylim(c(0,200))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  geom_abline(aes(intercept = 0, slope = 1))+
  xlab("Fitted (Copies/L) x 1000")+
  ylab("Observed (Copies/L) x 1000")+
  theme(axis.title = element_text(face="bold"))

plot=plot+geom_text(x =50, y = 150, aes(label = txt),parse=T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))

#plot

library(ggpubr)
ggarrange(plot,cor.plot,
          #labels=c("Goodness of Fit", "Imposed Matern Correlation Values"),
          ncol=1,nrow=2, legend = "top")+
  theme(plot.margin = margin(1,1,1,1, "cm")) 

summary(Cs.inla.null)

autoplot(Cs.inla.null)
#### end ####




dccs2=dispersion_check(Cs.inla.null)
plot(dccs2)

par(mfrow = c(1, 3))
plot(spat.inla.d$marginals.fix[[1]], type = "l", xlim=c(0,20),ylim=c(0,.2),
     xlab = "Intercept", ylab = "Density")

summary(Cs.inla.null)


##### Summary Statistics and ANCOVA ####
#Final Summary Stats

library(tidyverse)
library(ggplot2)
setwd("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis")
data=read.csv("MASTER.FINAL.ddPCRdata.Cshasta2021.csv")
#view(data)
# # of runs - # of NTCs
length(unique(data$Sample.ID))

# # of FB samples
length(grep("FB",data$Sample.ID))

# Chinook Field blank cocns.
data$Chinook.c.well[(grep("FB",data$Sample.ID))]

#C.shasta Field blank concs
data$C.shasta.c.well[(grep("FB",data$Sample.ID))]

# #3 was rerun twice (#2 and #1) and showed negative. Error likely associated 
# with lab work

CSD2=read.csv("CSD.csv")
#CSD2=cbind(CSD2,Depth)


CSD2$C.shasta.vol.cor.conc=replace_na(round(CSD2$C.shasta.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc=replace_na(round(CSD2$Chinook.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.upper=replace_na(round(CSD2$C.shasta.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.upper=replace_na(round(CSD2$Chinook.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.lower=replace_na(round(CSD2$C.shasta.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.lower=replace_na(round(CSD2$Chinook.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)

setwd("C:/Users/dylan/OneDrive/Desktop/Thesis/GS/GradSchool/GRAD SCHOOL/Project/Thesis")

library(tidyverse)

#write.csv(CSD2,file="CSD2.csv")
CSD2=read.csv("CSD.csv")
#CSD2=cbind(CSD2,Depth)


CSD2$C.shasta.vol.cor.conc=replace_na(round(CSD2$C.shasta.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc=replace_na(round(CSD2$Chinook.c.well*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.upper=replace_na(round(CSD2$C.shasta.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.upper=replace_na(round(CSD2$Chinook.c.well.upper*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)


CSD2$C.shasta.vol.cor.conc.lower=replace_na(round(CSD2$C.shasta.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)
CSD2$Chinook.vol.cor.conc.lower=replace_na(round(CSD2$Chinook.c.well.lower*CSD2$Elution.volum..uL./15*CSD2$Volume..ml./1000,0),0)

KLMData = CSD2
KLMData.KLM=subset(KLMData,SiteN==1)


KLMData$C.shasta.Copies.Second=KLMData$C.shasta.vol.cor.conc*KLMData$SubQ..L.s.
KLMData$Chinook.Copies.Second=KLMData$Chinook.vol.cor.conc*KLMData$SubQ..L.s.

KLMData$C.shasta.Copies.Second.low=KLMData$C.shasta.vol.cor.conc.lower*KLMData$SubQ..L.s.
KLMData$C.shasta.Copies.Second.high=KLMData$C.shasta.vol.cor.conc.upper*KLMData$SubQ..L.s.
KLMData$Chinook.Copies.Second.low=KLMData$Chinook.vol.cor.conc.lower*KLMData$SubQ..L.s.
KLMData$Chinook.Copies.Second.high=KLMData$Chinook.vol.cor.conc.upper*KLMData$SubQ..L.s.
# MAKE TOTALS FOR C.shasta

Total.Copies.Second<-KLMData%>%
  group_by(SiteN,Type,Distance.Left.Bank.m) %>%
  summarise(Sums = mean(C.shasta.Copies.Second),Q=mean(Q.),Depth=max(Point.Depth),
            D.Thal=mean(D.thal.m),low=mean(C.shasta.Copies.Second.low),high=mean(C.shasta.Copies.Second.high))

Total.Copies.Second=subset(Total.Copies.Second, Type != "DI")

Total.Copies.Second<-KLMData%>%
  group_by(SiteN,Distance.Left.Bank.m) %>%
  summarise(Sums = mean(C.shasta.Copies.Second),Q=mean(Q.),Depth=max(Point.Depth),
            D.Thal=mean(D.thal.m),low=mean(C.shasta.Copies.Second.low),high=mean(C.shasta.Copies.Second.high))

#View(Total.Copies.Second)

Total.Copies.Second.KLM=subset(Total.Copies.Second,SiteN==1)
Total.Copies.Second.I5=subset(Total.Copies.Second,SiteN==2)
Total.Copies.Second.TOH=subset(Total.Copies.Second,SiteN==3)
Total.Copies.Second.BVR=subset(Total.Copies.Second,SiteN==4)
Total.Copies.Second.KMN=subset(Total.Copies.Second,SiteN==5)
Total.Copies.Second.SV=subset(Total.Copies.Second,SiteN==6)



#MAKE TOTALS FOR CHINOOK


Total.Copies.Second.Cs<-KLMData%>%
  group_by(SiteN,Type,Distance.Left.Bank.m) %>%
  summarise(Sums = mean(Chinook.Copies.Second),Q=mean(Q.),Depth=max(Point.Depth),
            D.Thal=mean(D.thal.m),low=mean(Chinook.Copies.Second.low),high=mean(Chinook.Copies.Second.high))

Total.Copies.Second.Cs=subset(Total.Copies.Second.Cs, Type != "DI")

Total.Copies.Second.Cs<-KLMData%>%
  group_by(SiteN,Distance.Left.Bank.m) %>%
  summarise(Sums = mean(Chinook.Copies.Second),Q=mean(Q.),Depth=max(Point.Depth),
            D.Thal=mean(D.thal.m),low=mean(Chinook.Copies.Second.low),high=mean(Chinook.Copies.Second.high))

#View(Total.Copies.Second.Cs)

Total.Copies.Second.Cs.KLM=subset(Total.Copies.Second.Cs,SiteN==1)
Total.Copies.Second.Cs.I5=subset(Total.Copies.Second.Cs,SiteN==2)
Total.Copies.Second.Cs.TOH=subset(Total.Copies.Second.Cs,SiteN==3)
Total.Copies.Second.Cs.BVR=subset(Total.Copies.Second.Cs,SiteN==4)
Total.Copies.Second.Cs.KMN=subset(Total.Copies.Second.Cs,SiteN==5)
Total.Copies.Second.Cs.SV=subset(Total.Copies.Second.Cs,SiteN==6)





boxplot.data=read.csv("Boxplot.data.csv")

boxplot.data$SiteN=as.factor(boxplot.data$SiteN)
boxplot.data$Species=as.factor(boxplot.data$Species)
boxplot.data$Method=as.factor(boxplot.data$Method)
coef=.000000021



boxplot.data.C.shasta=filter(boxplot.data, Species=="C.shasta")
boxplot.data.Chinook=filter(boxplot.data, Species=="Chinook")
rm(SiteN)


dist=read.csv("2021FIELD_DATA_Klamath_Discharge.csv")



dist$Q.<-dist$Q.internal*0.028316846592
dist$Width.Interval.m<-dist$Width.Interval*0.3048

dist2<- dist %>% group_by(Site) %>%
  arrange(Site) %>% 
  mutate(w = cumsum(Width.Interval.m)) %>%
  mutate(wm=w-Width.Interval.m)%>%
  mutate(Q.total=sum(Q.))
#View(dist2)

dist2$wt <- with(dist2, wm + (w - wm)/2)

Site<-unique(dist2$Site)
Site<-paste(Site, sep=" ", collapse=NULL)
cmsv<-unique(dist2$Q.total)
cmsv<-round(cmsv,2)
cmsv<-paste(cmsv, sep=" ", collapse=NULL)
cmsv<-paste(cmsv,"m^3/s", collapse=NULL)
Distance.Left.Bank.Meters<-c(1,1,1,1,1,1)
Discharge.per.cell<-c(1,1,1,1,1,1)
Q.total.L=dist2$Q.total*1000



cms<-cbind(Site,cmsv,Distance.Left.Bank.Meters,Discharge.per.cell,Q.total.L)
cms<-as.data.frame(cms)
#cms

dist2$Site = factor(dist2$Site, levels =c("KLM","I5","TOH","BVR","KMN","SV"))


dist2.KLM=subset(dist2,Site=="KLM")
#dist2.KLM=dist2.KLM[2:20,]
#View(dist2.KLM)
dist2.KLM$Width.Interval.m[1:3]=0
#dist2.KLM$Width.Interval.m[19]=4


KLMData.KLM.G<-subset(KLMData.KLM,Type=="G")
KLMData.KLM.G["Distance.Left.Bank.m"][KLMData.KLM.G["Distance.Left.Bank.m"]== 0] <-0.1
KLMData.KLM.G=arrange(KLMData.KLM.G,Distance.Left.Bank.m)
KLMData.KLM.P<-subset(KLMData.KLM,Type=="P")







par(mfrow=c(7,1))
op = par(font = 2)
opar = par(no.readonly = TRUE)
par(mar = c(1,20,1,1))
name.1="Depth (m)"
name.2=c(0.8,0.6,0.2)
name.3 = c("","At Depth", "At Surface")
name.4 = "C.shasta"
name.5 = "Chinook"
#name.6=c(0,10,20,30,40)

# 
# 
#### Additional XS plots ####
# ticks=seq(0,40,by=5)
# 
# plot(x=rev(dist2.KLM$Distance.Left.Bank.meters),y=dist2.KLM$Depth*0.3048,xlim=c(0,40),
#      ylim = rev(range(dist2.KLM$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=rev(dist2.KLM$Distance.Left.Bank.meters),y=dist2.KLM$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,40),ylim=c(0,1.75)
# )+
#   lines(x=rev(dist2.KLM$Distance.Left.Bank.meters),y=dist2.KLM$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=rev(dist2.KLM$Distance.Left.Bank.meters),y=as.numeric(dist2.KLM$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.7,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# text(x=3,y=1.2,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.KLM.P$Distance.Left.Bank.m-5,y=KLMData.KLM.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,35),xlim=c(0,40))+
#   lines(x=KLMData.KLM.G$Distance.Left.Bank.m-5,y=KLMData.KLM.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.7,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# 
# text(x=37,y=30,labels=xprs2, cex=1.8)
# 
# # legend(x="left", legend = name.3,
# #        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.KLM.P$Distance.Left.Bank.m-5,y=KLMData.KLM.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,4000),xlim=c(0,40),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.KLM.G$Distance.Left.Bank.m-5,y=KLMData.KLM.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.7,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# text(x=37,y=3500,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("40.82" ~ (m^3/s)))
# 
# barplot(width=rev(dist2.KLM$Width.Interval.m),height=rev(dist2.KLM$Q.),xlim=c(0,40),space=0,
#         bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n",col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# text(x=4,y=4,labels=xp,cex=1.2)
# 
# 
# sum(Total.Copies.Second.KLM$Sums)
# 
# xprs4=expression(bold( atop("",atop("139,029,674 Copies/Second","Suspended eDNA Load"))))
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=rev(dist2.KLM$Width.Interval.m),height = Total.Copies.Second.KLM$Sums,space=0, xlim=c(0,40),
#         ylab="",xlab="",bty="n",xaxt="n",yaxt="n", col = "black")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.8,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=25,y=20000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.KLM$Sums)
# 
# xprs5=expression(bold( atop("",atop(" 48,845,897 Copies/Second","Suspended eDNA Load"))))
# 
# barplot(width=rev(dist2.KLM$Width.Interval.m),height=Total.Copies.Second.Cs.KLM$Sums,space=0,
#         xlim=c(0,40),
#         ylab="",bty="n",col="white",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=25,y=6500000,labels=xprs5, cex=1.8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##############I5
# 
# 
# dist2.I5=subset(dist2,Site=="I5")
# #dist2.KLM=dist2.KLM[2:20,]
# #View(dist2.I5)
# #dist2.KLM$Width.Interval.m[1:3]=0
# #dist2.KLM$Width.Interval.m[19]=4
# 
# KLMData.I5=subset(KLMData,SiteN==2)
# 
# KLMData.I5.G<-subset(KLMData.I5,Type=="G")
# KLMData.I5.G["Distance.Left.Bank.m"][KLMData.I5.G["Distance.Left.Bank.m"]== 0] <-0.1
# KLMData.I5.G=arrange(KLMData.I5.G,Distance.Left.Bank.m)
# KLMData.I5.P<-subset(KLMData.I5,Type=="P")
# 
# 
# 
# 
# 
# 
# 
# par(mfrow=c(7,1))
# op = par(font = 2)
# opar = par(no.readonly = TRUE)
# par(mar = c(1,20,1,1))
# name.1="Depth (m)"
# name.2=c(0.8,0.6,0.2)
# name.3 = c("","At Depth", "At Surface")
# name.4 = "C.shasta"
# name.5 = "Chinook"
# #name.6=c(0,10,20,30,40)
# 
# 
# ticks=seq(0,35,by=5)
# 
# plot(x=dist2.I5$Distance.Left.Bank.meters-3,y=dist2.I5$Depth*0.3048,xlim=c(0,37),
#      ylim = rev(range(dist2.I5$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.16, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.3, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=dist2.I5$Distance.Left.Bank.meters-3,y=dist2.I5$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,37),ylim=c(0,1.75)
# )+
#   lines(x=dist2.I5$Distance.Left.Bank.meters-3,y=dist2.I5$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=dist2.I5$Distance.Left.Bank.meters-3,y=as.numeric(dist2.I5$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.7,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.165, 0),xpd = TRUE)
# 
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# 
# text(x=3,y=1.4,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.I5.P$Distance.Left.Bank.m-5,y=KLMData.I5.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,150),xlim=c(0,37))+
#   lines(x=KLMData.I5.G$Distance.Left.Bank.m-5,y=KLMData.I5.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.9,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.18, 0),xpd = TRUE)
# #
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# 
# text(x=3,y=125,labels=xprs2, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.I5.P$Distance.Left.Bank.m-5,y=KLMData.I5.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,2000),xlim=c(0,37),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.I5.G$Distance.Left.Bank.m-5,y=KLMData.I5.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.9,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.18, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# text(x=3,y=1500,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# 
# #View(dist2.I5)
# dist2.I5$Width.Interval.m[1]=2
# dist2.I5$Width.Interval.m[16]=2
# cmsv
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("35.08" ~ (m^3/s)))
# 
# barplot(width=dist2.I5$Width.Interval.m,height=dist2.I5$Q.,xlim=c(0,37), space=0,
#         bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n",col = "darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.9,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# text(x=5,y=3,labels=xp,cex=1.2)
# 
# sum(Total.Copies.Second.I5$Sums)
# xprs4=expression(bold( atop("",atop("1,948,627,474 Copies/Second","Suspended eDNA Load"))))
# 
# 
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=dist2.I5$Width.Interval.m,height = Total.Copies.Second.I5$Sums, xlim=c(0,37),
#         space=0,ylab="",xlab="",bty="n",xaxt="n",yaxt="n",col="black")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.8,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# text(x=30,y=200000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.I5$Sums)
# 
# xprs4=expression(bold( atop("",atop("295,373,718 Copies/Second","Suspended eDNA Load"))))
# 
# 
# 
# barplot(width=dist2.I5$Width.Interval.m,height=Total.Copies.Second.Cs.I5$Sums,
#         space=0, xlim=c(0,37),
#         lty="solid",ylab="",bty="n",col="white",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        col="darkgrey", horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.3, 0),xpd = TRUE)
# 
# 
# text(x=30,y=3000000,labels=xprs5, cex=1.8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##############TOH
# 
# 
# dist2.TOH=subset(dist2,Site=="TOH")
# #dist2.KLM=dist2.KLM[2:20,]
# #View(dist2.TOH)
# #dist2.KLM$Width.Interval.m[1:3]=0
# #dist2.KLM$Width.Interval.m[19]=4
# 
# KLMData.TOH=subset(KLMData,SiteN==3)
# 
# KLMData.TOH.G<-subset(KLMData.TOH,Type=="G")
# KLMData.TOH.G["Distance.Left.Bank.m"][KLMData.TOH.G["Distance.Left.Bank.m"]== 0] <-0.1
# KLMData.TOH.G=arrange(KLMData.TOH.G,Distance.Left.Bank.m)
# KLMData.TOH.P<-subset(KLMData.TOH,Type=="P")
# 
# 
# 
# 
# 
# par(mfrow=c(7,1))
# op = par(font = 2)
# opar = par(no.readonly = TRUE)
# par(mar = c(1,20,1,1))
# name.1="Depth (m)"
# name.2=c(0.8,0.6,0.2)
# name.3 = c("","At Depth", "At Surface")
# name.4 = "C.shasta"
# name.5 = "Chinook"
# #name.6=c(0,10,20,30,40)
# 
# ticks=seq(0,35,by=5)
# 
# plot(x=dist2.TOH$Distance.Left.Bank.meters,y=dist2.TOH$Depth*0.3048,xlim=c(0,37),
#      ylim = rev(range(dist2.TOH$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=dist2.TOH$Distance.Left.Bank.meters,y=dist2.TOH$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,37),ylim=c(0,1.8)
# )+
#   lines(x=dist2.TOH$Distance.Left.Bank.meters,y=dist2.TOH$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=dist2.TOH$Distance.Left.Bank.meters,y=as.numeric(dist2.TOH$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.7,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.2, 0),xpd = TRUE)
# #
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# 
# text(x=3,y=1.4,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.TOH.P$Distance.Left.Bank.m,y=KLMData.TOH.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,270),xlim=c(0,37))+
#   lines(x=KLMData.TOH.G$Distance.Left.Bank.m,y=KLMData.TOH.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# text(x=3,y=200,labels=xprs2, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.TOH.P$Distance.Left.Bank.m,y=KLMData.TOH.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,1200),xlim=c(0,37),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.TOH.G$Distance.Left.Bank.m,y=KLMData.TOH.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# #
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# 
# text(x=3,y=1100,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# cms
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("36.77" ~ (m^3/s)))
# 
# #View(dist2.TOH)
# dist2.TOH$Width.Interval.m[1:2]=2
# dist2.TOH$Width.Interval.m[18]=2
# 
# barplot(width=dist2.TOH$Width.Interval.m,height =dist2.TOH$Q.,xlim=c(0,37),
#         col="darkgrey",bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n",space=0)
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.9,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=3,y=3,labels=xp,cex=1.2)
# 
# sum(Total.Copies.Second.TOH$Sums)
# 
# 
# xprs4=expression(bold( atop("",atop("2,325,105,445 Copies/Second","Suspended eDNA Load"))))
# 
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=dist2.TOH$Width.Interval.m,height=Total.Copies.Second.TOH$Sums, xlim=c(0,37),
#         ylab="",xlab="",bty="n",xaxt="n",yaxt="n",col="black",space=0)
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.8,
#        horiz = F, bty = "n",inset = c(-.2, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=4,y=350000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.TOH$Sums)
# 
# xprs5=expression(bold( atop("",atop("16,937,433 Copies/Second","Suspended eDNA Load"))))
# 
# 
# barplot(width=dist2.TOH$Width.Interval.m,height=Total.Copies.Second.Cs.TOH$Sums,
#         xlim=c(0,37),space=0,
#         col="white",ylab="",bty="n",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        horiz = F, bty = "n",inset = c(-.2, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=4,y=1800000,labels=xprs5, cex=1.8)
# 
# 
# 
# 
# 
# 
# 
# 
# ##############BVR
# 
# 
# dist2.BVR=subset(dist2,Site=="BVR")
# #dist2.KLM=dist2.KLM[2:20,]
# #View(dist2.BVR)
# #dist2.KLM$Width.Interval.m[1:3]=0
# #dist2.KLM$Width.Interval.m[19]=4
# 
# KLMData.BVR=subset(KLMData,SiteN==4)
# 
# KLMData.BVR.G<-subset(KLMData.BVR,Type=="G")
# KLMData.BVR.G["Distance.Left.Bank.m"][KLMData.BVR.G["Distance.Left.Bank.m"]== 0] <-0.1
# KLMData.BVR.G=arrange(KLMData.BVR.G,Distance.Left.Bank.m)
# KLMData.BVR.P<-subset(KLMData.BVR,Type=="P")
# 
# #View(dist2.BVR)
# dist2.BVR$Width.Interval.m[1]=2
# dist2.BVR$Width.Interval.m[16]=0
# 
# 
# 
# 
# 
# 
# 
# 
# 
# par(mfrow=c(7,1))
# op = par(font = 2)
# opar = par(no.readonly = TRUE)
# par(mar = c(1,20,1,1))
# name.1="Depth (m)"
# name.2=c(0.8,0.6,0.2)
# name.3 = c("","At Depth", "At Surface")
# name.4 = "C.shasta"
# name.5 = "Chinook"
# #name.6=c(0,10,20,30,40)
# 
# ticks=seq(0,40,by=5)
# 
# plot(x=dist2.BVR$Distance.Left.Bank.meters,y=dist2.BVR$Depth*0.3048,xlim=c(0,40),
#      ylim = rev(range(dist2.BVR$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=dist2.BVR$Distance.Left.Bank.meters,y=dist2.BVR$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,40),ylim=c(0,1.6)
# )+
#   lines(x=dist2.BVR$Distance.Left.Bank.meters,y=dist2.BVR$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=dist2.BVR$Distance.Left.Bank.meters,y=as.numeric(dist2.BVR$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.6,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# 
# text(x=3,y=1.4,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.BVR.P$Distance.Left.Bank.m-1,y=KLMData.BVR.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,165),xlim=c(0,40))+
#   lines(x=KLMData.BVR.G$Distance.Left.Bank.m-1,y=KLMData.BVR.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# 
# text(x=3,y=150,labels=xprs2, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.BVR.P$Distance.Left.Bank.m-1,y=KLMData.BVR.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,1200),xlim=c(0,40),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.BVR.G$Distance.Left.Bank.m-1,y=KLMData.BVR.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# 
# text(x=3,y=1100,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("37.12" ~ (m^3/s)))
# 
# barplot(width=dist2.BVR$Width.Interval.m,height=dist2.BVR$Q.,xlim=c(0,40),space=0,
#         bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n",col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.9,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# text(x=33,y=3,labels=xp,cex=1.2)
# 
# sum(Total.Copies.Second.BVR$Sums)
# 
# xprs4=expression(bold( atop("",atop("2,168,143,726 Copies/Second","Suspended eDNA Load"))))
# 
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=dist2.BVR$Width.Interval.m,height=Total.Copies.Second.BVR$Sums, xlim=c(0,40),space=0,
#         col = "black",ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.6,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=33,y=250000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.BVR$Sums)
# 
# xprs5=expression(bold( atop("",atop("12,598,145 Copies/Second","Suspended eDNA Load"))))
# 
# 
# barplot(width=dist2.BVR$Width.Interval.m,height=Total.Copies.Second.Cs.BVR$Sums,space = 0,
#         xlim=c(0,40),
#         ylab="",bty="n",col="white",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.6,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=33,y=1800000,labels=xprs5, cex=1.8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##############KMN
# 
# 
# dist2.KMN=subset(dist2,Site=="KMN")
# #dist2.KLM=dist2.KLM[2:20,]
# #View(dist2.BVR)
# #dist2.KLM$Width.Interval.m[1:3]=0
# #dist2.KLM$Width.Interval.m[19]=4
# 
# KLMData.KMN=subset(KLMData,SiteN==5)
# 
# 
# 
# KLMData.KMN.G<-subset(KLMData.KMN,Type=="G")
# KLMData.KMN.G["Distance.Left.Bank.m"][KLMData.KMN.G["Distance.Left.Bank.m"]== 0] <-0.1
# KLMData.KMN.G=arrange(KLMData.KMN.G,Distance.Left.Bank.m)
# KLMData.KMN.P<-subset(KLMData.KMN,Type=="P")
# 
# #View(dist2.KMN)
# dist2.KMN$Width.Interval.m[1]=2
# dist2.KMN$Width.Interval.m[18]=2
# 
# par(mfrow=c(7,1))
# op = par(font = 2)
# opar = par(no.readonly = TRUE)
# par(mar = c(1,20,1,1))
# name.1="Depth (m)"
# name.2=c(0.8,0.6,0.2)
# name.3 = c("","At Depth", "At Surface")
# name.4 = "C.shasta"
# name.5 = "Chinook"
# #name.6=c(0,10,20,30,35)
# 
# ticks=seq(0,35,by=5)
# 
# plot(x=dist2.KMN$Distance.Left.Bank.meters,y=dist2.KMN$Depth*0.3048,xlim=c(0,35),
#      ylim = rev(range(dist2.KMN$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.37, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=dist2.KMN$Distance.Left.Bank.meters,y=dist2.KMN$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,35),ylim=c(0,1.6)
# )+
#   lines(x=dist2.KMN$Distance.Left.Bank.meters,y=dist2.KMN$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=dist2.KMN$Distance.Left.Bank.meters,y=as.numeric(dist2.KMN$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.7,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# 
# text(x=3,y=1.4,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.KMN.P$Distance.Left.Bank.m-1,y=KLMData.KMN.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,200),xlim=c(0,35))+
#   lines(x=KLMData.KMN.G$Distance.Left.Bank.m-1,y=KLMData.KMN.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# 
# text(x=3,y=190,labels=xprs2, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.KMN.P$Distance.Left.Bank.m-1,y=KLMData.KMN.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,1200),xlim=c(0,35),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.KMN.G$Distance.Left.Bank.m-1,y=KLMData.KMN.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# 
# text(x=3,y=1100,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("44.03" ~ (m^3/s)))
# 
# barplot(width=dist2.KMN$Width.Interval.m,height=dist2.KMN$Q.,xlim=c(0,35),space=0,
#         col="darkgrey",bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.9,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=3,y=4,labels=xp,cex=1.2)
# 
# sum(Total.Copies.Second.KMN$Sums)
# 
# xprs4=expression(bold( atop("",atop("3,532,269,926 Copies/Second","Suspended eDNA Load"))))
# 
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=dist2.KMN$Width.Interval.m,height=Total.Copies.Second.KMN$Sums,lwd=4,  xlim=c(0,35),
#         space=0,col="black",ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.8,
#        horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=32,y=450000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.KMN$Sums)
# 
# xprs5=expression(bold( atop("",atop("13,978,102 Copies/Second","Suspended eDNA Load"))))
# 
# 
# barplot(width=dist2.KMN$Width.Interval.m,height=Total.Copies.Second.Cs.KMN$Sums,
#         xlim=c(0,35),space=0,
#         ylab="",bty="n",col="white",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        col="darkgrey", horiz = F, bty = "n",inset = c(-.18, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=32.5,y=1550000,labels=xprs5, cex=1.8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########################SV
# 
# dist2.SV=subset(dist2,Site=="SV")
# #dist2.KLM=dist2.KLM[2:20,]
# #View(dist2.BVR)
# #dist2.KLM$Width.Interval.m[1:3]=0
# #dist2.KLM$Width.Interval.m[19]=4
# 
# KLMData.SV=subset(KLMData,SiteN==6)
# 
# 
# 
# 
# KLMData.SV.G<-subset(KLMData.SV,Type=="G")
# KLMData.SV.G["Distance.Left.Bank.m"][KLMData.SV.G["Distance.Left.Bank.m"]== 0] <-0.1
# KLMData.SV.G=arrange(KLMData.SV.G,Distance.Left.Bank.m)
# KLMData.SV.P<-subset(KLMData.SV,Type=="P")
# 
# par(mfrow=c(7,1))
# op = par(font = 2)
# opar = par(no.readonly = TRUE)
# par(mar = c(1,20,1,1))
# name.1="Depth (m)"
# name.2=c(0.8,0.6,0.2)
# name.3 = c("","At Depth", "At Surface")
# name.4 = "C.shasta"
# name.5 = "Chinook"
# #name.6=c(0,10,20,30,40)
# 
# ticks=seq(0,40,by=5)
# 
# plot(x=dist2.SV$Distance.Left.Bank.meters,y=dist2.SV$Depth*0.3048,xlim=c(0,45),
#      ylim = rev(range(dist2.SV$Depth*0.3048)),lwd=4, type="b",bty="n",
#      xaxt = "n", yaxt = "n", xlab="", ylab = "",
# )
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2)
# legend(x="left", legend = "",
#        cex =0.6,
#        lty = 1,  lwd = 4, horiz = TRUE, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# legend(x="left", legend = name.1,
#        cex =1.1, bty="n",
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xprs=expression(bold( atop("",atop("Velocity (m/s)","Proportion of Depth"))))
# 
# 
# plot(x=dist2.SV$Distance.Left.Bank.meters,y=dist2.SV$fps.8*0.3048,
#      
#      lwd=4, type="l",bty="n",ylab="",xlab="",xaxt = "n", yaxt = "n",xlim=c(0,45),ylim=c(0,2)
# )+
#   lines(x=dist2.SV$Distance.Left.Bank.meters,y=dist2.SV$fps.6*0.3048,
#         lwd=4, type="l",bty="n",col="darkgrey")+
#   lines(x=dist2.SV$Distance.Left.Bank.meters,y=as.numeric(dist2.SV$fps.2)*0.3048,
#         lwd=4, type="l",lty="dashed")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2)
# 
# legend(x="left",y=1.5, legend = c("","",""),
#        title = "",  cex=0.7,
#        lty = c(1,1,3), col=c("black","darkgrey","black"), lwd = 4, horiz = F, bty = "n",
#        inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=1.5, legend = "",
# #        title = xprs,
# #        bty = "n",cex=1.8,
# #        inset = c(-.18, 0),xpd = TRUE)
# 
# 
# text(x=3,y=1.4,labels=xprs, cex=1.8)
# 
# legend(x="left",y=1.5, legend = name.2,
#        bty = "n",cex=1.0,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs2=expression(bold(atop("",atop("C.shasta Conc.","Thousands of Copies/L"))))
# 
# plot(x=KLMData.SV.P$Distance.Left.Bank.m-2,y=KLMData.SV.P$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt = "n", yaxt = "n",ylim=c(0,200),xlim=c(0,45))+
#   lines(x=KLMData.SV.G$Distance.Left.Bank.m-2,y=KLMData.SV.G$C.shasta.vol.cor.conc/1000,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# 
# # legend(x="left",y=0.9, legend="",
# #        title = xprs2,
# #        bty = "n",cex=1.8,
# #        inset = c(-.23, 0),xpd = TRUE)
# 
# text(x=3,y=175,labels=xprs2, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# xprs3=expression(bold( atop("",atop("Chinook salmon Conc.","Copies/L"))))
# 
# plot(x=KLMData.SV.P$Distance.Left.Bank.m-2,y=KLMData.SV.P$Chinook.vol.cor.conc,lwd=4, type="b",
#      lty="solid",ylab="",xlab="",bty="n",xaxt="n",yaxt="n",ylim = c(0,1200),xlim=c(0,45),
#      cex.axis=1.7,cex.lab=2,font=2)+
#   lines(x=KLMData.SV.G$Distance.Left.Bank.m-2,y=KLMData.SV.G$Chinook.vol.cor.conc,lwd=4, type="b",
#         lty="solid",ylab="",xlab="",bty="n", col="darkgrey")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = c("","",""),
#        title = "",
#        cex = 0.6,
#        lty = c(1,1,1), col=c("white","black","darkgrey"), lwd = 4, horiz = F, bty = "n", inset = c(-.2, 0),xpd = TRUE)
# #
# # legend(x="left",y=0.2, legend="",
# #        title = xprs3,
# #        bty = "n",cex=1.8,
# #        inset = c(-.21, 0),xpd = TRUE)
# 
# text(x=10,y=1100,labels=xprs3, cex=1.8)
# 
# legend(x="left", legend = name.3,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# xl<-expression(bold(Q ~ (m^3/s)))
# xp<-expression(bold("66.77" ~ (m^3/s)))
# 
# barplot(width=dist2.SV$Width.Interval.m,height = dist2.SV$Q.,xlim=c(0,45),space=0,
#         col="darkgrey",bty="n", ylab= "" ,xlab="",xaxt="n",yaxt="n")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.9,
#        title.adj = 0.5,
#        horiz = F, bty = "n",inset = c(-.18, -.17) ,xpd = TRUE)
# 
# legend(x="left", legend = xl,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=3,y=4,labels=xp,cex=1.2)
# sum(Total.Copies.Second.SV$Sums)
# 
# xprs4=expression(bold( atop("",atop("3,696,006,554 Copies/Second","Suspended eDNA Load"))))
# 
# 
# par(mar = c(2,20,1,1))
# 
# barplot(width=dist2.SV$Width.Interval.m,height=Total.Copies.Second.SV$Sums,xlim=c(0,45),space=0,
#         col="black",ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
# axis(side=1,las=0,labels=FALSE,at=ticks)
# axis(side = 2, las = 2, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.8,
#        horiz = F, bty = "n",inset = c(-.2, 0) ,xpd = TRUE)
# legend(x="left", legend = name.4,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# text(x=5,y=400000000,labels=xprs4, cex=1.8)
# 
# 
# 
# sum(Total.Copies.Second.Cs.SV$Sums)
# 
# xprs5=expression(bold( atop("",atop("11,420,657 Copies/Second","Suspended eDNA Load"))))
# 
# 
# barplot(width=dist2.SV$Width.Interval.m,height=Total.Copies.Second.Cs.SV$Sums,space=0,
#         xlim=c(0,45),
#         ylab="",bty="n",col="white",xlab="",yaxt ="n",xaxt="n")
# 
# axis(side=1,labels=T,cex.axis = 2, las=1, font = 2,at=ticks)
# axis(side = 2, las = 1, cex.axis = 2, font = 2,srt=35)
# legend("left", legend = "",
#        cex = 0.7,
#        col="darkgrey", horiz = F, bty = "n",inset = c(-.2, 0) ,xpd = TRUE)
# legend(x="left", legend = name.5,
#        bty = "n",cex=1.1,
#        inset = c(-.35, 0),xpd = TRUE)
# 
# 
# text(x=5,y=1500000,labels=xprs5, cex=1.8)
# 
# 
# 
# 


# 
# library(akima)
# library(graphics)
# library(metR)
# #KLM Sampling Site
# 
# 
# 
# yy=data.frame(KLMData.KLM$Distance.Left.Bank.m,KLMData.KLM$D.M.,KLMData.KLM$SiteN,
#               KLMData.KLM$C.shasta.vol.cor.conc,
#               KLMData.KLM$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# dev.off()
# par(mfrow=c(1,1))
# par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("grey","grey","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("grey", "darkgreen")),
#                nlevels = 10,
#                xlim = c(6,48),
#                ylim = c(-1.4,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(KLM)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(6,48),
#                ylim = c(-1.4,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(KLM)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# 
# 
# #I5 Sampling Site
# 
# yy=data.frame(KLMData.I5$Distance.Left.Bank.m,KLMData.I5$D.M.,KLMData.I5$SiteN,
#               KLMData.I5$C.shasta.vol.cor.conc,
#               KLMData.I5$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# yy
# 
# library(ggplot2)
# 
# geom_contour(color = "black", size = 0.1)
# 
# 
# #par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# #fld=data.frame(fld)
# pal = colorRampPalette(c("white","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkgreen")),
#                nlevels = 10,
#                xlim = c(2,38),
#                ylim = c(-1.0,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(I5)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(2,38),
#                ylim = c(-1.0,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(I5)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# 
# 
# #TREE OF HeAVEN
# 
# 
# yy=data.frame(KLMData.TOH$Distance.Left.Bank.m,KLMData.TOH$D.M.,KLMData.TOH$SiteN,
#               KLMData.TOH$C.shasta.vol.cor.conc,
#               KLMData.TOH$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# 
# 
# 
# #par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkgreen")),
#                nlevels = 10,
#                xlim = c(-1,38),
#                ylim = c(-1.2,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(TOH)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(-1,38),
#                ylim = c(-1.2,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(TOH)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# 
# 
# 
# #Beaver Creek
# 
# 
# yy=data.frame(KLMData.BVR$Distance.Left.Bank.m,KLMData.BVR$D.M.,KLMData.BVR$SiteN,
#               KLMData.BVR$C.shasta.vol.cor.conc,
#               KLMData.BVR$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# 
# 
# 
# #par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkgreen")),
#                nlevels = 10,
#                xlim = c(-1,40),
#                ylim = c(-1.8,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(BVR)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(-1,40),
#                ylim = c(-1.8,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(BVR)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# 
# 
# 
# #Kinsman
# 
# 
# yy=data.frame(KLMData.KMN$Distance.Left.Bank.m,KLMData.KMN$D.M.,KLMData.KMN$SiteN,
#               KLMData.KMN$C.shasta.vol.cor.conc,
#               KLMData.KMN$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# 
# 
# 
# #par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkgreen")),
#                nlevels = 10,
#                xlim = c(-1,36),
#                ylim = c(-1.4,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(KMN)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(-1,36),
#                ylim = c(-1.4,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(KMN)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# 
# 
# 
# #Seiad Valley
# 
# 
# yy=data.frame(KLMData.SV$Distance.Left.Bank.m,KLMData.SV$D.M.,KLMData.SV$SiteN,
#               KLMData.SV$C.shasta.vol.cor.conc,
#               KLMData.SV$Chinook.vol.cor.conc)
# colnames(yy)=c("X","Y","Z","Conc.Cs","Conc.Chin")
# 
# 
# #par(mar=c(4,4,4,4))
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Cs,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkgreen"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Cs, sort(yy$Conc.Cs))
# 
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkgreen")),
#                nlevels = 10,
#                xlim = c(-1,44),
#                ylim = c(-1.2,0.1),
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected C.shasta DNA Concentraion(SV)",
#                key.title = title(main = "CONC", cex.main = 1), 
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )
# 
# fld <- with(yy, interp(x = X, y = -Y, z = Conc.Chin,duplicate="mean"))
# 
# 
# pal = colorRampPalette(c("white","darkblue"))(10)
# pal=colorRampPalette(pal)
# yy$order = findInterval(yy$Conc.Chin, sort(yy$Conc.Chin))
# 
# zlim= range(fld$z, finite = TRUE)
# nlevels=10
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette =
#                  colorRampPalette(c("white", "darkblue")),
#                
#                xlim = c(-1,44),
#                ylim = c(-1.2,0.1),
#                zlim = zlim,
#                levels = pretty(zlim, nlevels),
#                nlevels = nlevels,
#                xlab = "Distance From Left Bank",
#                ylab = "Depth",
#                main = "Volume Corrected Chinook eDNA Concentraion(SV)",
#                key.title = title(main = "CONC", cex.main = 1),
#                plot.axes = { axis(1); axis(2); 
#                  points(yy$X,yy$Y*(-1),cex = 3,pch=21, 
#                         col="black",
#                         bg=pal(nrow(yy))[yy$order])
#                }
# )

#### end xs plots ####

hist(CSD2$Chinook.c.well)

summary(CSD2$Chinook.c.well)

length(CSD2$Chinook.c.well[CSD2$Chinook.c.well >= 11.8])
length(CSD2$Chinook.c.well[CSD2$Chinook.c.well >= 11.8])/length(CSD2$Chinook.c.well)*100

length(CSD2$C.shasta.c.well[CSD2$C.shasta.c.well >= 6.3])/length(CSD2$C.shasta.c.well)*100

cs.zeros=CSD2$C.shasta.c.well[CSD2$C.shasta.c.well <=6.3]
summary(CSD2$C.shasta.c.well)

summary(CSD2$C.shasta.vol.cor.conc)
sd(CSD2$C.shasta.vol.cor.conc)
summary(CSD2$Chinook.vol.cor.conc)
sd(CSD2$Chinook.vol.cor.conc)
#LOQ

length(CSD2$Chinook.c.well[CSD2$Chinook.c.well >= 16])
length(CSD2$Chinook.c.well[CSD2$Chinook.c.well >= 16])/length(CSD2$Chinook.c.well)*100

length(CSD2$C.shasta.c.well[CSD2$C.shasta.c.well >= 18])/length(CSD2$C.shasta.c.well)*100

hist.data=
  data.frame(c(CSD2$Chinook.vol.cor.conc,CSD2$C.shasta.vol.cor.conc.lower),
             c(rep("Chinook",192),rep("C.shasta",192)))
colnames(hist.data) = c("Concentration","Species")


p1<-ggplot(subset(hist.data,Species=="C.shasta"), aes(x=Concentration/1000) ) + 
  geom_histogram(binwidth=10, color ="black",fill="darkgreen") +
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))#+
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
# facet_wrap(~ Species)

p2<-ggplot(subset(hist.data,Species=="Chinook"), aes(x=Concentration/1000) ) + 
  geom_histogram(binwidth=0.1, color = "black",fill="darkblue") + 
  theme_classic()+
  
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))#+

# geom_vline(aes(xintercept=6.3),
#            linetype="solid",colour="#BB0000",size=1.5)+
# geom_vline(aes(xintercept=16),
#            linetype="dashed",colour="#0000FF",size=1.5)+
# facet_wrap(~ Species)
library(ggpubr)
histos=ggarrange(p1+ rremove("ylab")+ rremove("xlab"),p2+ rremove("ylab")+ rremove("xlab"),nrow = 2)

annotate_figure(histos, bottom = text_grob("Thousands of Copies/L", 
                                           color = "black", face = "bold", size = 14),
                left =text_grob("Frequency", rot=90,
                                color = "black", face = "bold", size = 14))



#ggarrange(cs.klm,chin.klm,cs.i5,chin.i5,cs.toh,chin.toh,cs.bvr,chin.bvr,cs.kmn,chin.kmn,cs.sv,chin.sv,
#    ncol=2,nrow=6,legend="right")
par(mfrow=c(1,1))
#plot(mesh2d)

#summary(nospat.inla.d.m)

#summary(nospat.inla.t.m)



library(sf)
library(raster)
library(rgeos)
library(sp)
library(maptools)
library(mgcv)
library(magrittr)
require(RColorBrewer)
require(rgdal)



#### make xs with smoothers ####
inuseData=KLMData.KLM

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.klm <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")

inuseData=KLMData.I5

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.i5 <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")

inuseData=KLMData.TOH

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)



Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.TOH <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")

inuseData=KLMData.BVR

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 



Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.BVR <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.KMN

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# 

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Blues'


PP.chin.KMN <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.SV

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.SV <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=Chinook.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=Chinook.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.KLM

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Greens'

PP.CS.klm <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.I5

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# 
# head(NewData)
# 

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 



Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.i5 <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.TOH

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.TOH <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.BVR

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)


Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.BVR <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.KMN

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 



Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.KMN <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.SV

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)

Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.SV <- ggplot() +
  geom_tile(data = NewData,
            aes(x=Distance.Left.Bank.m,
                y=D.M.,
                fill=C.shasta.vol.cor.conc)) + 
  geom_contour(data = NewData,
               aes(x=Distance.Left.Bank.m,
                   y=D.M.,
                   z=C.shasta.vol.cor.conc),
               color='grey30') + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



sites13CS.fig= ggarrange( 
  PP.CS.klm+ rremove("ylab")+ rremove("xlab"),
  PP.CS.TOH+ rremove("ylab")+ rremove("xlab"),                        
  ncol=1,nrow=2, common.legend = T, legend = "top"
)

annotate_figure(sites13CS.fig, bottom = text_grob("Distance From Left Bank (m)", 
                                                  color = "black", face = "bold", size = 14),
                left =text_grob("Depth (m) ", rot=90,
                                color = "black", face = "bold", size = 14),
                right=text_grob("Sites 1 - 6",rot = -90, color = "black", face = "bold", size = 14)#,
                # top = text_grob("C. shasta   &   Chinook salmon", 
                #color = "black", face = "bold", size = 14)
)


#### end section ####

##########MAKE The XS without the isocline ####



inuseData=KLMData.KLM

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

#
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
# #


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.klm <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.I5

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.i5 <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.TOH

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
# # 

Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.TOH <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.BVR

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 



Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.BVR <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.KMN

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit


# # 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 


Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.KMN <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.SV

gamThinPlateSpline <- gam(scale(Chinook.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$Chinook.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
#

Palette <- 'Blues'
# NewData=newcoords
#head(NewData)

PP.chin.SV <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(Chinook.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.KLM

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

# 
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)



Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.klm <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")




inuseData=KLMData.I5

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

#  
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
# # 


Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.i5 <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.TOH

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

#
inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
# # 
Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.TOH <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
 
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")

inuseData=KLMData.BVR

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)


Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.BVR <- ggplot() +
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
 
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")


inuseData=KLMData.KMN

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?

newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)
# 
#


Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.KMN <- ggplot() +
  # scale_y_reverse() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
 
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



inuseData=KLMData.SV

gamThinPlateSpline <- gam(scale(C.shasta.vol.cor.conc) ~ te(D.M.,Distance.Left.Bank.m,k=3), data = inuseData) ###Log transform and GAMS?


newDepths <- seq(0,max(inuseData$D.M.)/.6,0.02) #?
newDist <- seq(0,max(inuseData$Distance.Left.Bank.m),0.01) #?




NewData <- data.frame(Distance.Left.Bank.m = rep(newDist,each=length(newDepths)),
                      D.M. = rep(newDepths,length(newDist)))
predTemp <- predict(gamThinPlateSpline,newdata = NewData, se.fit=TRUE)
NewData$C.shasta.vol.cor.conc <- predTemp$fit ###You can back transform but it will look like it dropped off the cliff 
###Log base doesn't matter
NewData$se.fit <- predTemp$se.fit

inuseData2=subset(inuseData,Type01==1)
addon=data.frame(matrix(ncol=2,nrow=2))
addon$D.M.=c(0,0)
addon$Distance.Left.Bank.m=c(0,max(inuseData$Distance.Left.Bank.m)+1)
# 
inuseData2=full_join(inuseData2,addon)



Palette <- 'Greens'
# NewData=newcoords
#head(NewData)

PP.CS.SV <- ggplot() +
  geom_line(data=inuseData2,
            aes(x=Distance.Left.Bank.m,
                y=D.M./.6),size=1.2,show.legend=F)+
  geom_point(data = inuseData,
             aes(x=Distance.Left.Bank.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=D.M.+0.1,#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
                 fill=scale(C.shasta.vol.cor.conc)),
             color='black',
             pch=21,
             size=5) + 
  # scale_y_reverse() +
  
  scale_fill_gradientn(colors=brewer.pal(9,Palette)) + 
  #scale_x_continuous(breaks=seq(3,1,-1),labels=c('L','C','R')) +
  #scale_x_continuous() +
  scale_y_reverse() +
  theme_classic()+
  xlab('Distance Left Bank (m)') +
  ylab('Depth') +
  labs(fill= "Site-Scaled Copies/L") +
  theme(legend.position="top")



#### end section ####
all.pp.fig <- ggarrange(
  PP.CS.klm+ rremove("ylab")+ rremove("xlab") ,#+ ggtitle('Site 1'),
  PP.chin.klm+ rremove("ylab")+ rremove("xlab"),
  PP.CS.i5+ rremove("ylab")+ rremove("xlab"),# + ggtitle('Site 2'),
  PP.chin.i5+ rremove("ylab")+ rremove("xlab"),
  PP.CS.TOH+ rremove("ylab")+ rremove("xlab"),# + ggtitle('Site 3'),
  PP.chin.TOH+ rremove("ylab")+ rremove("xlab"),
  PP.CS.BVR+ rremove("ylab")+ rremove("xlab"),# + ggtitle('Site 4'),
  PP.chin.BVR + rremove("ylab")+ rremove("xlab"),
  PP.CS.KMN+ rremove("ylab")+ rremove("xlab"),# + ggtitle('Site 5'),
  PP.chin.KMN+ rremove("ylab")+ rremove("xlab"),
  PP.CS.SV+ rremove("ylab")+ rremove("xlab"),# + ggtitle('Site 6'),
  PP.chin.SV+ rremove("ylab")+ rremove("xlab"),
  ncol=2,nrow=6, common.legend = T, legend = "top"
)

all.pp.fig

annotate_figure(all.pp.fig, bottom = text_grob("Distance From Left Bank (m)",
                                               color = "black", face = "bold", size = 14),
                left =text_grob("Depth (m) ", rot=90,
                                color = "black", face = "bold", size = 14),
                right=text_grob(" ",rot = -90, color = "black", face = "bold", size = 14),
                top = text_grob("C. shasta   &   Chinook salmon",
                                color = "black", face = "bold", size = 14))





###################MAKE THE DATA SUMMARY FIGURE ####

# 
# hist.data=
#   data.frame(c(CSD2$Chinook.vol.cor.conc,CSD2$C.shasta.vol.cor.conc.lower),
#              c(rep("Chinook",192),rep("C.shasta",192)))
# colnames(hist.data) = c("Concentration","Species")


pp1<-ggplot(CSD2, aes(x=D.M.) ) + 
  geom_histogram( color ="black",fill="grey", bins=10) +
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp1

pp1.1<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=D.M.,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=C.shasta.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkgreen',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp1.1

pp1.2<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=D.M.,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=Chinook.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkblue',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp1.2






pp2<-ggplot(CSD2, aes(x=DistanceFromShore.m) ) + 
  geom_histogram( color ="black",fill="grey", bins=10) +
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp2

pp2.1<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=DistanceFromShore.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=C.shasta.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkgreen',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp2.1

pp2.2<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=DistanceFromShore.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=Chinook.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkblue',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp2.2



pp3<-ggplot(CSD2, aes(x=D.thal.m) ) + 
  geom_histogram( color ="black",fill="grey", bins=10) +
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp3

pp3.1<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=D.thal.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=C.shasta.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkgreen',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp3.1

pp3.2<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=D.thal.m,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=Chinook.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkblue',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp3.2




pp4<-ggplot(CSD2, aes(x=V.M.) ) + 
  geom_histogram( color ="black",fill="grey", bins=10) +
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp4

pp4.1<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=V.M.,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=C.shasta.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkgreen',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp4.1

pp4.2<-ggplot() + 
  geom_point(data = CSD2,
             aes(x=V.M.,#+sign(replicate-2)*(replicate%%2)*0.025,  #####Whats going on here
                 y=Chinook.vol.cor.conc/1000),#+as.numeric(replicate%%2==0)*10, ###Modulus operator is basically the remainder-- take a number divided by another and wahtever is left over is the modulus. in the case of this it is either even or odd. It is good for shit like this
             #fill=scale(C.shasta.vol.cor.conc)),
             color='darkblue',
             pch=19,
             size=3) + 
  theme_classic()+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
# geom_vline(aes(xintercept=11.8),
#            linetype="solid",colour="#BB0000",size=2)+
# geom_vline(aes(xintercept=18),
#            linetype="dashed",colour="#0000FF",size=2)+
#facet_wrap(~ Species)
pp4.2

labz=c("A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","D3","D4")

library(ggpubr)
ggarrange(pp1+ rremove("ylab")+ rremove("xlab"),pp1.1+ rremove("ylab")+ rremove("xlab"),
          pp1.2+ rremove("ylab")+ rremove("xlab"),
          pp2+ rremove("ylab")+ rremove("xlab"),pp2.1+ rremove("ylab")+ rremove("xlab"),pp2.2+ rremove("ylab")+ rremove("xlab"),
          pp3+ rremove("ylab")+ rremove("xlab"),pp3.1+ rremove("ylab")+ rremove("xlab"),pp3.2+ rremove("ylab")+ rremove("xlab"),
          pp4+ rremove("ylab")+ rremove("xlab"),pp4.1+ rremove("ylab")+ rremove("xlab"),pp4.2+ rremove("ylab")+ rremove("xlab"),
          ncol=3,nrow=4,labels=labz, hjust=-9)+
  theme(plot.margin = margin(1,1,1,1, "cm")) 




#### end section ####







######### Violin Plots ####

boxplot.data$mean=boxplot.data$Q.DNA/(boxplot.data$Q*1000)

par(mar=c(4,4,2,2))
par(mfrow=c(2,2))






library(vioplot)
x1 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==1]
x2 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==2]
x3 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==3]
x4 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==4]
x5 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==5]
x6 <- CSD2$C.shasta.vol.cor.conc[CSD2$SiteN==6]
vio.cs.1=vioplot(x1, x2, x3,x4,x5,x6, names=c(1:6),
                 col="darkgreen", cex.axis=1.5, cex.lab = 10,
                 ylab="Grab-Sample C. shasta")

boxplot.data.C.shasta.Q=subset(boxplot.data.C.shasta,Method=="Q.DNA")
boxplot.data.C.shasta.Q$mean=boxplot.data.C.shasta.Q$Q.DNA/(boxplot.data.C.shasta.Q$Q*1000)
vio.cs.1



y1 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==1]
y2 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==2]
y3 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==3]
y4 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==4]
y5 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==5]
y6 <- CSD2$Chinook.vol.cor.conc[CSD2$SiteN==6]
vio.chin.1=vioplot(y1, y2, y3,y4,y5,y6, names=c(1:6),cex.axis=1.5,cex.lab = 2,
                   col="lightblue",xlab="", 
                   ylab="Grab-Sample O. tshawytscha")


boxplot.data.Chinook.Q=subset(boxplot.data.Chinook,Method=="Q.DNA")
boxplot.data.Chinook.Q$mean=boxplot.data.Chinook.Q$Q.DNA/(boxplot.data.Chinook.Q$Q*1000)
vio.chin.1#+



B=2*15000
V.target.15000=(B^2)/4
B.2=2*3000
V.target.3000=(B.2^2)/4
Sig.sq.=CSD2%>%group_by(SiteN)%>%
  summarize(N=length(C.shasta.vol.cor.conc),mean=mean(C.shasta.vol.cor.conc),
            sig.sq=(sum((C.shasta.vol.cor.conc-mean)^2))/(N-1), CV=sqrt(sig.sq)/mean*100,
            n.target.15000=ceiling((sig.sq)/((V.target.15000)+(sig.sq/N))),
            n.target.3000=ceiling((sig.sq)/((V.target.3000)+(sig.sq/N))))

B=2*107
V.target.107=(B^2)/4
B.2=2*214
V.target.214=(B.2^2)/4
Sig.sq.Chinook=CSD2%>%group_by(SiteN)%>%
  summarize(N=length(Chinook.vol.cor.conc),mean=mean(Chinook.vol.cor.conc),
            sig.sq=(sum((Chinook.vol.cor.conc-mean)^2))/(N-1), CV=sqrt(sig.sq)/mean*100,
            n.target.214=ceiling((sig.sq)/((V.target.214)+(sig.sq/N))),
            n.target.107=ceiling((sig.sq)/((V.target.107)+(sig.sq/N)))
  )



X2021FIELD_DATA_Klamath_2=read.csv("2021FIELD_DATA_Klamath.2.csv")

DI2=X2021FIELD_DATA_Klamath_2[grep("_DI_",X2021FIELD_DATA_Klamath_2$Sample.ID),]

DI2$C.shasta.vol.cor.conc=replace_na(round(DI2$Cs.Conc.Eye*20*DI2$Elution.volum..uL./15*DI2$Volume..ml./1000,0),0)
DI2$Chinook.vol.cor.conc=replace_na(round(DI2$Chin.Conc.Eye*20*DI2$Elution.volum..uL./15*DI2$Volume..ml./1000,0),0)
DI2$SiteN<-as.factor(DI2$Site)
DI2$SiteN<-recode_factor(DI2$SiteN,KLM=1,I5=2,TOH=3,BVR=4,KMN=5,SV=6)


B=2*15000
V.target.15000=(B^2)/4
B.2=2*3000
V.target.3000=(B.2^2)/4
Sig.sq.DI2=DI2%>%group_by(SiteN)%>%
  summarize(n.s=length(C.shasta.vol.cor.conc),N=n.s*3, mean=mean(C.shasta.vol.cor.conc),
            sig.sq=(sum((C.shasta.vol.cor.conc-mean)^2))/(n.s-1), CV=sqrt(sig.sq)/mean*100,
            n.target.15000=ceiling((sig.sq)/((V.target.15000)+(sig.sq/N))),
            n.target.3000=ceiling((sig.sq)/((V.target.3000)+(sig.sq/N))))


B=2*107
V.target.107=(B^2)/4
B.2=2*214
V.target.214=(B.2^2)/4

Sig.sq.Chinook.DI2=DI2%>%group_by(SiteN)%>%
  summarize(n.s=length(Chinook.vol.cor.conc),N=n.s*3,mean=mean(Chinook.vol.cor.conc),
            sig.sq=(sum((Chinook.vol.cor.conc-mean)^2))/(n.s-1), CV=sqrt(sig.sq)/mean*100,
            n.target.214=ceiling((sig.sq)/((V.target.214)+(sig.sq/N))),
            n.target.107=ceiling((sig.sq)/((V.target.107)+(sig.sq/N))))



x1 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==1]
x2 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==2]
x3 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==3]
x4 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==4]
x5 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==5]
x6 <- DI2$C.shasta.vol.cor.conc[DI2$SiteN==6]
vio.cs.2=vioplot(x1, x2, x3,x4,x5,x6, 
                 names=c(1:6), ylim= c(0,250000),cex.axis=1.5,
                 col="darkgreen", ylab="Depth-Integrated C. shasta")
vio.cs.2



y1 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==1]
y2 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==2]
y3 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==3]
y4 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==4]
y5 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==5]
y6 <- DI2$Chinook.vol.cor.conc[DI2$SiteN==6]
vio.chin.2=vioplot(y1, y2, y3,y4,y5,y6, names=c(1:6),ylim=c(0,3000),
                   col="lightblue",ylab="Depth-Integrated O. tshawytscha", cex.axis=1.5,
                   xlab="")
vio.chin.2

#### end section ####

#### ANCOVA for treatment type ####

#make the dataset
species=c(rep("C.shasta",6),rep("O.tshawytcha",6))
method=c(rep("Grab Samples",12),rep("Depth Integrated Samples",12))
one=dplyr::full_join(Sig.sq.,Sig.sq.Chinook)
two=dplyr::full_join(Sig.sq.DI2,Sig.sq.Chinook.DI2)
#one=one[,1:5]
one$target=species
#two=two[,1:6]
#two=two[,-2]
two$target=species
two$SiteN=as.integer(two$SiteN)
three=dplyr::full_join(one,two)

three$method=method
three$n=c(three$N[1:12],6,5,5,4,5,6,6,5,5,4,5,6)

three=three%>%relocate(n, .after = N)

summary(three)

colnames(three)=c("Site","N","n","mu","sig.sq","CV","n.target.15000",
                  "n.target.3000","n.target.214","n.target.107","target",
                  "n.s","method")

library(dplyr)
three%>%
  group_by(target) %>%  
  summarise(mean = mean(CV),
            sd = sd(CV),
            mean_mean = mean(mu),
            sd_mean = sd(mu))
par(mfrow=c(1,1))
boxplot(CV~method,
        data=three)


## Check assumptions

library(car)


three$Site=as.factor(three$Site)
three$target=as.factor(three$target)
three$method=as.factor(three$method)



leveneTest(CV~target,data=three)
leveneTest(CV~method,data=three)
leveneTest(CV~Site, data=three)

#pass the homogeneity of variance assumption

#Fit the model
ancova_model<-aov(CV~method+target+Site,data=three)
car::Anova(ancova_model,type="III",test.statistic="Chisq")


quantile(three$CV, prob=0.0001)

colors=c("darkgreen","darkblue")
library(ggplot2)

ppp= ggplot(three,aes(x=as.factor(method),y=CV,fill=target))+
  theme_classic()+
  geom_boxplot()+scale_fill_manual(values=colors, name="")+labs(title="",
                                                                
                                                                x="Method", 
                                                                y = "Coefficient of Variation % (CV) ")



ppp+ theme(legend.position="top",axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),legend.text = element_text(size=16),
           axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))

library(multcomp)

postHocs2 <- multcomp::glht(ancova_model, linfct = mcp(method = "Tukey"))

summary(postHocs)
confint(postHocs)

postHocs2 <- multcomp::glht(ancova_model, linfct = mcp(Site = "Tukey"))

summary(postHocs2)
confint(postHocs2)



head(three)

three.cs=subset(three, target=="C.shasta")
head(three.cs)
three.cs=three.cs[,-9]
three.cs=three.cs[,-9]
three.cs=three.cs[,-9]
three.cs=three.cs[,-9]

library(kableExtra)
kbl(three.cs, caption = "<b>Table 3: <i>C. shasta</i>
    optimal sample size calculation for 6 sites, where 
    <i>N</i> is the total number of eDNA samples possible per sampling method, <i>n</i> is the sample size, 
    $\\mu$ is the mean or estimate of the mean, $\\sigma^{2}$ is the 
    finite population variance or estimate, CV is the coefficient of variation as a percent, $n_o$ is
    the optimal sample size for the specified standard error (SE), and method is the sampling strategy.<b/>",
    col.names = c("Site", "N","n", "$\\mu$ (copies/L)","$\\sigma^{2}$","CV%",
                  "$n_o$ (SE = 15,000)","$n_o$ (SE = 3,000)","method"))%>%  
  kable_classic(full_width = F, html_font = "Cambria")



three.ot=subset(three, target=="O.tshawytcha")
three.ot=three.ot[,-7]
three.ot=three.ot[,-7]
three.ot=three.ot[,-9]

three.ot=three.ot[,-9]

library(kableExtra)
kbl(three.ot, caption = "<b>Table 4: <i>O. tshawytcha</i>
    optimal sample size calculation for 6 sites, where 
    <i>N</i> is the total number of eDNA samples possible per sampling method, <i>n</i> is the sample size, 
    $\\mu$ is the mean or estimate of the mean, $\\sigma^{2}$ is the 
    finite population variance or estimate, CV is the coefficient of variation as a percent, $n_o$ is
    the optimal sample size for the specified standard error (SE), and method is the sampling strategy.<b/>",
    col.names = c("Site", "N","n", "$\\mu$ (copies/L)","$\\sigma^{2}$","CV%",
                  "$n_o$ (SE = 214)","$n_o$ (SE = 107)","method"))%>%  
  kable_classic(full_width = F, html_font = "Cambria")



#### end section ####


## In order to show the paired differences between surface and depth samples


bardata=CSD2%>%group_by(Distance.Left.Bank.m,SiteN,Type01)%>%summarize(Chin=mean(Chinook.vol.cor.conc))
bardata2=bardata%>%group_by(Distance.Left.Bank.m,SiteN)%>%
  summarize(Cs.diff=Chin[1]-Chin[2])



hist(as.numeric(bardata2$Cs.diff),na.rm=T, breaks=100)


bardata2=drop_na(bardata2)
bin=ifelse(bardata2$Chin>0,1,0)


zee_ave=bardata2%>%group_by(SiteN)%>%
  summarize(bin=mean(ifelse(Cs.diff>0,1,0)),percent=bin*100,label=paste0(round(percent,1),"%"))

ggplot(bardata2, aes(x=Distance.Left.Bank.m, y=Cs.diff))+
  geom_bar(stat='identity', fill="darkblue",width=2)+theme_classic()+
  facet_wrap(~SiteN,nrow = 3,ncol = 2,scales="free")+geom_abline(slope=0,intercept=0)+
  #geom_smooth(aes(x=Distance.Left.Bank.m, y=Cs.diff),linetype="dashed", color="red",se=F)+
  labs(title="",
       x="Distance From Left Bank (m)", 
       y = "Surface - Depth")+
  theme(legend.position="top",axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),legend.text = element_text(size=16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))+theme(strip.background = element_blank())+
  # Add individual text to plot
  geom_text(data = zee_ave,size=5,
            mapping = aes(x = c(15,10,7,10,10,10),
                          y = c(3000,1000,1000,800,750,200),
                          label = label))



#### end ####


#### marginal calcs ####
## chose a marginal and compare the with the results computed by the
## inla-program
r = nospat.inla.d.m$summary.fixed["x",]
m = nospat.inla.d.m$marginals.fixed$DistanceFromShore.m.z.z

## compute the 95% HPD interval
inla.hpdmarginal(0.95, m)

x = seq(-6, 6, len = 1000)
y = dnorm(x)
inla.hpdmarginal(0.95, list(x=x, y=y))

## compute the the density for exp(r), version 1
r.exp = inla.tmarginal(exp, m)
## or version 2
r.exp = inla.tmarginal(function(x) exp(x), m)

## to plot the marginal, we use the inla.smarginal, which interpolates (in
## log-scale). Compare with some samples.
plot(inla.smarginal(m), type="l",xlab="marginal",ylab="probability")
abline(v=0,col="red")

s = inla.rmarginal(1000, m)
hist(inla.rmarginal(1000, m), add=TRUE, prob=TRUE)
lines(density(s), lty=2)

m1 = inla.emarginal(function(x) x^1, m)
m2 = inla.emarginal(function(x) x^2, m)
stdev = sqrt(m2 - m1^2)
q = inla.qmarginal(c(0,0.999), m)

## inla-program results
print(r)

## inla.marginal-results (they shouldn't be perfect!)
print(c(mean=m1, sd=stdev, "0.0quant" = q[1], "0.999quant" = q[2]))
## using the buildt-in function
inla.zmarginal(m)

#### end ####