#1. Packages ----

library(MASS); library(scales); library(gridExtra)
library(methods); library(devtools); library(pbkrtest)
library(car); library(colorspace); library(fields)
library(tidyverse); library(broom); library(mgcv)
library(dplyr); library(mgcv); library(lsmeans)
library(ggplot2); library(performance)

#2. Load data ----
#file = "Stroudetal_GitHub_data.csv"
res=read.csv(file.choose(""))

#6. Fitness surface ----

#load the magic that is fields
library(fields)

##6a. Prepare data----
#estimate fitness (survival probability) across all morphospace
#select axes
x = cbind(res$LD1, res$LD3)
#turn it into a matrix
x2 = as.matrix(x, ncol=2)
#turn your Z axis variable into a vector
#this is your metric of fitness, here binomial survival data
survival = c(res$survival)

##6b. Estimate surface----
#fit survival data to the morph space
#using a thin plate spline (Tps)
#this can take some time depending on your computer
#mine is terrible, so this is always painful
comm.fit = Tps(x2, survival)
#check model details
summary(comm.fit)
#what's the edf of the surface?
comm.fit$eff.df
#ooh that's high!
#quick look at the surface
surface(comm.fit)
#surface looks cool!

#predict true surface
comm.out = predictSurface(comm.fit)
comm.out

##6c. Plot surface----
###6c(i) 2D surface----
#take a look at the selection surface
par(mfrow = c(1,1))
#2D surface plot
image.plot(comm.out, xlab = "", 
           ylab = "", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1), 
           ylim = c(min(res$LD3)-1, max(res$LD3)+1))

#add contours so the topography is a little bit clearer
contour(comm.out, add=T, col = "gray50", nlevels=12,
        drawlabels = FALSE
        #,drawlabels=TRUE #you can add labels to the contour lines if you like
)
#convex hulls of distributions
#sag
sag.dist = data.frame(subset(res, species == "sagrei")$LD1, 
                      subset(res, species == "sagrei")$LD3)
sag.chx <- chull(sag.dist)
sag.hull <- rbind(sag.dist = sag.dist[sag.chx, ], sag.dist[sag.chx[1], ])
lines(sag.hull, lwd = 4, col = "brown")

#caro
caro.dist = data.frame(subset(res, species == "carolinensis")$LD1, 
                       subset(res, species == "carolinensis")$LD3)
caro.chx <- chull(caro.dist)
caro.hull <- rbind(caro.dist = caro.dist[caro.chx, ], caro.dist[caro.chx[1], ])
lines(caro.hull, lwd = 4, col = "green")

#dist
dist.dist = data.frame(subset(res, species == "distichus")$LD1, 
                       subset(res, species == "distichus")$LD3)
dist.chx <- chull(dist.dist)
dist.hull <- rbind(dist.dist = dist.dist[dist.chx, ], dist.dist[dist.chx[1], ])
lines(dist.hull, lwd = 4, col = "gray50")

#equ
eque.dist = data.frame(subset(res, species == "equestris")$LD1, 
                       subset(res, species == "equestris")$LD3)
eque.chx <- chull(eque.dist)
eque.hull <- rbind(eque.dist = eque.dist[eque.chx, ], eque.dist[eque.chx[1], ])
lines(eque.hull, lwd = 4, col = "black")

#add species centroids
head(sag.dist)
points(mean(sag.dist[,1]), mean(sag.dist[,2]), pch=21, bg="chocolate", cex=3)
points(mean(dist.dist[,1]), mean(dist.dist[,2]), pch=21, bg="gray75", cex=3)
points(mean(eque.dist[,1]), mean(eque.dist[,2]), pch=21, bg="black", cex=3)
points(mean(caro.dist[,1]), mean(caro.dist[,2]), pch=21, bg="green", cex=3)

#find local peaks and plot those points
res$fit = comm.fit$fitted.values
sag.fit = droplevels(subset(res, species == "sagrei"))
sag.peak = sag.fit[which.max(sag.fit$fit),]
points(sag.peak$LD1, sag.peak$LD3, pch=23, bg="white", cex=3)

car.fit = droplevels(subset(res, species == "carolinensis"))
car.peak = car.fit[which.max(car.fit$fit),]
points(car.peak$LD1, car.peak$LD3, pch=23, bg="white", cex=3)

dist.fit = droplevels(subset(res, species == "distichus"))
dist.peak = dist.fit[which.max(dist.fit$fit),]
points(dist.peak$LD1, dist.peak$LD3, pch=23, bg="white", cex=3)

eque.fit = droplevels(subset(res, species == "equestris"))
eque.peak = eque.fit[which.max(eque.fit$fit),]
points(eque.peak$LD1, eque.peak$LD3, pch=23, bg="white", cex=3)

###6c(ii) 3D surface----

#3D plot
par(mfrow=c(1,1))
drape.plot = drape.plot(comm.out, zlab = "", xlab = "", ylab = "",
                        #zlim = c(0.025,0.325), 
                        theta=255, phi=20, 
                        add.legend=F,
                        xlim = c(min(res$LD1)-1, max(res$LD1)+2), 
                        ylim = c(min(res$LD3)-1, max(res$LD3)+2))

#use the pushpin function to locate species mean locations
pushpin(mean(sag.dist[,1]), mean(sag.dist[,2]), 0.25,
        p.out=drape.plot, height=0.15,
        pch=21, col = "brown", cex=4
        #,text = "TrGr"
)
pushpin(mean(caro.dist[,1]), mean(caro.dist[,2]), 0.24,
        p.out=drape.plot, height=0.15,
        col = "darkgreen", cex=4
        #,text = "TrCr", pos=2
)
pushpin(mean(dist.dist[,1]), mean(dist.dist[,2]), 0.26,
        p.out=drape.plot, height=0.15,
        col = "darkgray", cex=4
        #,text = "Tr", pos=2
)
pushpin(mean(eque.dist[,1]), mean(eque.dist[,2]), 0.22,
        p.out=drape.plot, height=0.15,
        col = "black", cex=4
        #,text = "CrGi", pos=2
)

#plot points
res$z1 = 0.025
res$z2 = 0.025
sag = droplevels(subset(res, species == "sagrei"))
caro = droplevels(subset(res, species == "carolinensis"))
dist = droplevels(subset(res, species == "distichus"))
eque = droplevels(subset(res, species == "equestris"))

mypoints <- trans3d(res$LD1, res$LD3, res$z1, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)

#sagrei
mypoints <- trans3d(sag$LD1, sag$LD3, sag$z1, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(sag, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z2, pmat=drape.plot)
points(mypoints, pch=21, bg="brown", cex=2.5)

#distichus
mypoints <- trans3d(dist$LD1, dist$LD3, dist$z1, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(dist, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z2, pmat=drape.plot)
points(mypoints, pch=21, bg="darkgray", cex=2.5)

#equestris
mypoints <- trans3d(eque$LD1, eque$LD3, eque$z1, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(eque, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z2, pmat=drape.plot)
points(mypoints, pch=21, bg="black", cex=2.5)

#carolinensis
mypoints <- trans3d(caro$LD1, caro$LD3, caro$z1, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(caro, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z2, pmat=drape.plot)
points(mypoints, pch=21, bg="darkgreen", cex=2.5)

#3D surface plot

#cool!

##6d. LDA figure----
#### . Fig 1: LDA ----
par(mfrow=c(1,1))
plot(res$LD1, res$LD3,col = 'white', xlab = "Linear Discriminant Axis 1 (79.9%)",
     ylab = "Linear Discriminant Axis 3 (7.7%)", 
     xlim = c(min(res$LD1)-1, max(res$LD1)+1), 
     ylim = c(min(res$LD3)-1, max(res$LD3)+2))
#draw convex hull around the whole community
lines(sag.hull, lwd = 4, col = "brown")
lines(dist.hull, lwd = 4, col = "gray50")
lines(caro.hull, lwd = 4, col = "forestgreen")
lines(eque.hull, lwd = 4, col = "black")
points(sag$LD1, sag$LD3, pch=21, bg = "chocolate")
points(caro$LD1, caro$LD3, pch=21, bg = "green")
points(dist$LD1, dist$LD3, pch=21, bg = "gray75")
points(eque$LD1, eque$LD3, pch=21, bg = "black")

#### . LDA w/ arrows ----
plot(lda.values$x[,1], lda.values$x[,3], 
     col="gray90", 
     cex = 1.5,
     xlim = c(-17.5,10),
     ylim = c(-15,7.5),
     xlab = "Linear Discriminant Axis 1",
     ylab = "Linear Discriminant Axis 3")
lines(sag.hull, lwd = 1, col = "brown", lty=2)
lines(dist.hull, lwd = 1, col = "gray50", lty=2)
lines(caro.hull, lwd = 1, col = "forestgreen", lty=2)
lines(eque.hull, lwd = 1, col = "black", lty=2)
lda.arrows(ld, col = "gray50", 
           myscale = 6)

#### . Other LD axes ----
library(plyr)
find_hull <- function(res) res[chull(res$LD1, res$LD2), ]
hulls <- ddply(res, "species", find_hull)
ld1.ld2 <- ggplot(data = res, aes(x = LD1, y = LD2, colour=species, fill = species)) +
  geom_point() + 
  geom_polygon(data = hulls, alpha = 0.5) +
  labs(x = "Linear discriminant axis 1 (79.9%)", 
       y = "Linear discriminant axis 2 (12.5%)")+
  scale_fill_manual(values=c("forestgreen", "gray75", "black", "chocolate"))+
  scale_color_manual(values=c("forestgreen", "gray75", "black", "brown"))+
  theme_bw()
ld1.ld2

find_hull <- function(res) res[chull(res$LD2, res$LD3), ]
hulls <- ddply(res, "species", find_hull)
ld2.ld3 <- ggplot(data = res, aes(x = LD2, y = LD3, colour=species, fill = species)) +
  geom_point() + 
  geom_polygon(data = hulls, alpha = 0.5) +
  labs(x = "Linear discriminant axis 2 (12.5%)", 
       y = "Linear discriminant axis 3 (7.7%)")+
  scale_fill_manual(values=c("forestgreen", "gray75", "black", "chocolate"))+
  scale_color_manual(values=c("forestgreen", "gray75", "black", "brown"))+
  theme_bw()
ld2.ld3

find_hull <- function(res) res[chull(res$LD1, res$LD3), ]
hulls <- ddply(res, "species", find_hull)
ld1.ld3 <- ggplot(data = res, aes(x = LD1, y = LD3, colour=species, fill = species)) +
  geom_point() + 
  geom_polygon(data = hulls, alpha = 0.5) +
  labs(x = "Linear discriminant axis 2 (12.5%)", 
       y = "Linear discriminant axis 3 (7.7%)")+
  scale_fill_manual(values=c("forestgreen", "gray75", "black", "chocolate"))+
  scale_color_manual(values=c("forestgreen", "gray75", "black", "brown"))+
  theme_bw()
ld1.ld3


library(ggpubr)
ggarrange(ld1.ld2, ld2.ld3, ld1.ld3, ncol =3, nrow=1,
          labels = c("(a)","(b)", "(c)"))

#7. Null model ----
#comparing empirical surface with null model simulation data 

#see separate code for how to run the permutation null model
#here we will compare the curvature of the empirical selection surface
#relative to the curvature of surfaces calculated from randomly assigned survival data
#survival data always match the species-level observed survival rates

#file = "communitysurface_nullmodel_permutations.csv"

sims = read.csv(file.choose(""))
str(sims)

#get the effective degrees of freedom from the empirical model
empirical.edfs = comm.fit$eff.df
empirical.edfs

#see how many simulations produced surfaces less complex
#than the empirical estimate
similar.sims = length(subset(sims, edf.ests < empirical.edfs)$edf.ests)
#calculate simple p-value
p = 1-(similar.sims/length(sims$edf.ests)) 
p
#p=0.0344

#plot simulated data vs. (red line) empirical surface complexity
#simple plot
ggplot(sims, aes(x=edf.ests))+
  geom_density(color="darkblue", fill="lightblue")+
  xlim(0,max(sims$edf.ests))+
  geom_vline(xintercept = empirical.edfs, lty=2, col = "red")+
  theme_bw()+
  ylab("")+
  xlab("")


#fancy plot if needed
ggplot(sims, aes(x=edf.ests))+
  geom_density(color="darkblue", fill="lightblue")+
  xlim(0,max(sims$edf.ests))+
  geom_vline(xintercept = 11.0623, lty=2, col = "red")+
  theme_bw()+
  ylab("Frequency 
of simulations")+
  theme(axis.title.y = element_text(angle = 0))+
  xlab("Surface curvature
(Effective Degrees of Freedom [EDF] of Thin Plate Spline surfaces)")+
  annotate(
    "text", label = "Empirical surface",
    x = empirical.edfs+6.5, y = 0.68, size = 4.5, colour = "red")

#ok, so the empirical surface is incredibly unlikely to have arisen by random survival
#i.e. by chance


#8. Species fitness surfaces---- 

#to check the community wide pattern wasn't just a weird phenomenon
#driven by the unique situation of all 4 species data being combined
#we tested the hypothesis that the fitness surfaces for each species (estimated independently)
#should each individually represent what is estimate on the community-wide surface

#fitness surfaces for each species independently 
sag = droplevels(subset(res, species == "sagrei"))
caro = droplevels(subset(res, species == "carolinensis"))
dist = droplevels(subset(res, species == "distichus"))
eque = droplevels(subset(res, species == "equestris"))

#community wide minimum convex polygon
comm.dist = data.frame(res$LD1, res$LD3)
comm.chx <- chull(comm.dist)
comm.hull <- rbind(comm.dist = comm.dist[comm.chx, ], comm.dist[comm.chx[1], ])

par(mfrow = c(1,2))

##8a. Anolis sagrei----
# sagrei #
x = cbind(sag$LD1, sag$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag$survival)
sag.fit = Tps(x2, survival)

sag.out = predictSurface(sag.fit)
#2D surface plot
image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3",
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1)
)
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3",
           xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1)
)
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(sag.hull, lwd = 2, col = "black")


drape.plot(sag.out, zlab = "Fitness
(Survival prob.)", xlab = "LD1", ylab = "LD3",
           zlim = c(-0.05,0.4), 
           theta=10, phi=25, 
           xlim = c(min(res$LD1)+4, max(res$LD1)+2),
           ylim = c(min(res$LD3)+4, max(res$LD3)),
           horiz = F)

#add survival points
#sagrei
sag$z11 <- 0
mypoints <- trans3d(sag$LD1, sag$LD3, sag$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(sag, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="brown", cex=2.5)


##8b. Anolis carolinensis----
# carolinensis #
x = cbind(caro$LD1, caro$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro$survival)
caro.fit = Tps(x2, survival)
caro.out = predictSurface(caro.fit)
#2D surface plot
image.plot(caro.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1))
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

image.plot(caro.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3",
           xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1)
)
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(caro.hull, lwd = 2, col = "black")


drape.plot(caro.out, zlab = "Fitness
           (Survival prob.)", xlab = "LD1", ylab = "LD3",
           zlim = c(0,0.4), 
           theta=10, phi=25, 
           xlim = c(min(res$LD1)+2, max(res$LD1)-2),
           ylim = c(min(res$LD3)+2, max(res$LD3)-2), 
           lwd=0.0001, horiz = T, lwd=0.3)
#carolinensis
caro$z11 <- 0
mypoints <- trans3d(caro$LD1, caro$LD3, caro$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(caro, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="darkgreen", cex=2.5)

##8c. Anolis distichus----
# distichus #
x = cbind(dist$LD1, dist$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(dist$survival)
dist.fit = Tps(x2, survival)
dist.out = predictSurface(dist.fit)
#2D surface plot
image.plot(dist.out, xlab = "LD1: Body size and foot length", 
           ylab = "LD3: Limb dimensions", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1))
contour(dist.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")
points(mean(dist$LD1), mean(dist$LD3), pch=21, bg="lightgray", cex=1.5)

drape.plot(dist.out, zlab = "Fitness
(Survival prob.)", xlab = "LD1", ylab = "LD3",
           zlim = c(0,0.5), 
           theta=10, phi=25, 
           xlim = c(min(res$LD1), max(res$LD1)),
           ylim = c(min(res$LD3), max(res$LD3)), horiz = T)
#distichus
dist$z11 <- 0
mypoints <- trans3d(dist$LD1, dist$LD3, caro$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(dist, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray40", cex=2.5)


##8d. Anolis equestris----
# equestris #
x = cbind(eque$LD1, eque$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(eque$survival)
eque.fit = Tps(x2, survival)
eque.out = predictSurface(eque.fit)
#2D surface plot
image.plot(eque.out, xlab = "LD1: Body size and foot length", 
           ylab = "LD3: Limb dimensions", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           zlim = c(-0.05,0.625))
contour(eque.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")
points(mean(eque$LD1), mean(eque$LD3), pch=21, bg="black", cex=1.5)

drape.plot(eque.out, zlab = "Fitness
(Survival prob.)", xlab = "LD1", ylab = "LD3",
           #zlim = c(0,0.4), 
           theta=10, phi=25, 
           xlim = c(min(res$LD1), max(res$LD1)),
           ylim = c(min(res$LD3), max(res$LD3)), horiz = T)
#equestris
eque$z11 <- 0
mypoints <- trans3d(eque$LD1, eque$LD3, caro$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray95", cex=2.5)
ss = droplevels(subset(eque, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z11, pmat=drape.plot)
points(mypoints, pch=21, bg="gray40", cex=2.5)

##8e. Multi-plot comparison----
#now compare community-wide surface vs. combination of individual surfaces
par(mfrow = c(1,2))

#community-wide
image(comm.out, xlab = "LD1: Body size and foot length", 
      ylab = "LD3: Limb dimensions",
      xlim = c(min(res$LD1)-1, max(res$LD1)+1),
      ylim = c(min(res$LD3)-1, max(res$LD3)+1),
      main = "Community-wide surface",
      col = tim.colors())
#add contours so the topography is a little bit clearer
contour(comm.out, add=T, col = "gray50", nlevels=10,
        drawlabels = FALSE
        #,drawlabels=TRUE #you can add labels to the contour lines if you like
)
lines(sag.hull, col = "brown", lwd=3)
lines(caro.hull, col = "darkgreen", lwd=3)
lines(dist.hull, col = "gray50", lwd=3)
lines(eque.hull, col = "black", lwd=3)

#individuals
#all of the individual surfaces are
#set to the same colour scale for fitness
#that's why some are not as red as when plotted individually above
plot(res$LD1, res$LD3,
     xlab = "", 
     ylab = "", 
     xlim = c(min(res$LD1), max(res$LD1)), 
     ylim = c(min(res$LD3), max(res$LD3)),
     main = "", type="n")
lines(comm.hull, lty=2, col = "black")
image(sag.out, add=T, col = tim.colors(), zlim = c(-0.05,0.625))
contour(sag.out, add=T, col = "gray50", nlevels=10,
        drawlabels = FALSE)
lines(sag.hull, col = "brown", lwd=3)
image(eque.out, add=T, col = tim.colors(), zlim = c(-0.05,0.625))
contour(eque.out, add=T, col = "gray50", nlevels=10,
        drawlabels = FALSE)
lines(eque.hull, col = "black", lwd=3)
image(caro.out, add=T, col = tim.colors(), zlim = c(-0.05,0.625))
contour(caro.out, add=T, col = "gray50", nlevels=10,
        drawlabels = FALSE)
lines(caro.hull, col = "green", lwd=3)
image(dist.out, add=T, col = tim.colors(), zlim = c(-0.05,0.625))
contour(dist.out, add=T, col = "gray50", nlevels=10,
        drawlabels = FALSE)
lines(dist.hull, col = "gray50", lwd=3)

#add species centroids
points(mean(sag.dist[,1]), mean(sag.dist[,2]), pch=21, bg="chocolate", cex=3)
points(mean(dist.dist[,1]), mean(dist.dist[,2]), pch=21, bg="gray75", cex=3)
points(mean(eque.dist[,1]), mean(eque.dist[,2]), pch=21, bg="black", cex=3)
points(mean(caro.dist[,1]), mean(caro.dist[,2]), pch=21, bg="green", cex=3)

#find local peaks and plot those points
points(sag.peak$LD1, sag.peak$LD3, pch=23, bg="white", cex=3)
points(car.peak$LD1, car.peak$LD3, pch=23, bg="white", cex=3)
points(dist.peak$LD1, dist.peak$LD3, pch=23, bg="white", cex=3)
points(eque.peak$LD1, eque.peak$LD3, pch=23, bg="white", cex=3)

#new individual selection surface peaks
#sagrei
sag$fit = sag.fit$fitted.values
sag.peak2 = sag[which.max(sag$fit),]
points(sag.peak2$LD1, sag.peak2$LD3, pch=24, bg="red", cex=2)

#carolinensis
caro$fit = caro.fit$fitted.values
caro.peak2 = caro[which.max(caro$fit),]
points(caro.peak2$LD1, caro.peak2$LD3, pch=24, bg="red", cex=2)

#distichus
dist$fit = dist.fit$fitted.values
dist.peak2 = dist[which.max(dist$fit),]
points(dist.peak2$LD1, dist.peak2$LD3, pch=24, bg="red", cex=2)

#equestris
eque$fit = eque.fit$fitted.values
eque.peak2 = eque[which.max(eque$fit),]
points(eque.peak2$LD1, eque.peak2$LD3, pch=24, bg="red", cex=2)

#highly congruent!
#nice.

#9. Temporal variation----
library(mgcv); library(lme4)

#first, let's see the sagrei peak

## 9a. Anolis sagrei ----
par(mfrow = c(1,1))
sag$z1 = -0.01
pmat = drape.plot(sag.out, 
                  zlim2 = c(0,0.65),
                  add.legend=F,
                  xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
                  ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
                  zlim=c(-0.01,0.3),
                  xlab="",ylab="",zlab="",
                  theta = 25, phi=25)
mypoints <- trans3d(sag$LD1, sag$LD3, sag$z1, pmat=pmat)
points(mypoints, pch=21, bg="gray95", cex=2.5)

ss = droplevels(subset(sag, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z1, pmat=pmat)
points(mypoints, pch=21, bg="gray40", cex=2.5)

#ok, now we are going to visualise the distribution of fitness through time
#and to do this, we will estimate fitness surfaces
#for each sampling period independently

par(mfrow=c(2,3))
image.plot(sag.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Cumulative")
contour(sag.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)


###9a(i). Surfaces through time ----
sag.f15 = droplevels(subset(sag, sampling.session == "fall.2015"))
x = cbind(sag.f15$LD1, sag.f15$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f15$survival)
sag.15.fit = Tps(x2, survival)
sag.f15.out = predictSurface(sag.15.fit)
image.plot(sag.f15.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2015")
contour(sag.f15.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)

#0.35
sag.s16 = droplevels(subset(sag, sampling.session == "spring.2016"))
x = cbind(sag.s16$LD1, sag.s16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.s16$survival)
fit = Tps(x2, survival)
sag.s16.out = predictSurface(fit)
image.plot(sag.s16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Spring 2016")
contour(sag.s16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)

#0.35
sag.f16 = droplevels(subset(sag, sampling.session == "fall.2016"))
x = cbind(sag.f16$LD1, sag.f16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f16$survival)
fit = Tps(x2, survival)
sag.f16.out = predictSurface(fit)
image.plot(sag.f16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2016")
contour(sag.f16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
#0.6

sag.s17 = droplevels(subset(sag, sampling.session == "spring.2017"))
x = cbind(sag.s17$LD1, sag.s17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.s17$survival)
fit = Tps(x2, survival)
sag.s17.out = predictSurface(fit)
image.plot(sag.s17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Spring 2017")
contour(sag.s17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
#0.25

sag.f17 = droplevels(subset(sag, sampling.session == "fall.2017"))
x = cbind(sag.f17$LD1, sag.f17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f17$survival)
fit = Tps(x2, survival)
sag.f17.out = predictSurface(fit)
image.plot(sag.f17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2017")
contour(sag.f17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
#0.35

#cool.

##9b. Anolis carolinensis ----
par(mfrow = c(1,1))
caro$z1 = -0.01
pmat = drape.plot(caro.out, 
                  zlim2 = c(0,0.65),
                  add.legend=F,
                  xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
                  ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
                  zlim=c(-0.01,0.3),
                  xlab="",ylab="",zlab="",
                  theta = 25, phi=25)
mypoints <- trans3d(caro$LD1, caro$LD3, caro$z1, pmat=pmat)
points(mypoints, pch=21, bg="gray95", cex=2.5)

ss = droplevels(subset(caro, survival == "1"))
mypoints <- trans3d(ss$LD1, ss$LD3, ss$z1, pmat=pmat)
points(mypoints, pch=21, bg="gray40", cex=2.5)

#ok, now we are going to visualise the distribution of fitness through time
#and to do this, we will estimate fitness surfaces
#for each sampling period independently

par(mfrow=c(2,3))
image.plot(caro.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Cumulative")
contour(caro.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)


###9b(i). Surfaces through time ----
caro.f15 = droplevels(subset(caro, sampling.session == "fall.2015"))
x = cbind(caro.f15$LD1, caro.f15$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f15$survival)
caro.15.fit = Tps(x2, survival)
caro.f15.out = predictSurface(caro.15.fit)
image.plot(caro.f15.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2015")
contour(caro.f15.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)

#0.35
caro.s16 = droplevels(subset(caro, sampling.session == "spring.2016"))
x = cbind(caro.s16$LD1, caro.s16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.s16$survival)
fit = Tps(x2, survival)
caro.s16.out = predictSurface(fit)
image.plot(caro.s16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Spring 2016")
contour(caro.s16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)

#0.35
caro.f16 = droplevels(subset(caro, sampling.session == "fall.2016"))
x = cbind(caro.f16$LD1, caro.f16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f16$survival)
fit = Tps(x2, survival)
caro.f16.out = predictSurface(fit)
image.plot(caro.f16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2016")
contour(caro.f16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
#0.6

caro.s17 = droplevels(subset(caro, sampling.session == "spring.2017"))
x = cbind(caro.s17$LD1, caro.s17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.s17$survival)
fit = Tps(x2, survival)
caro.s17.out = predictSurface(fit)
image.plot(caro.s17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Spring 2017")
contour(caro.s17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
#0.25

caro.f17 = droplevels(subset(caro, sampling.session == "fall.2017"))
x = cbind(caro.f17$LD1, caro.f17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f17$survival)
fit = Tps(x2, survival)
caro.f17.out = predictSurface(fit)
image.plot(caro.f17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2017")
contour(caro.f17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
#0.35

#cool.

#10. Survival rates----

##10a. Anolis sagrei----
#survival percentage through time
sag$sampling.session <- factor(sag$sampling.session, levels = 
                                 c("fall.2015", 
                                   "spring.2016", "fall.2016", 
                                   "spring.2017", "fall.2017"))
tab.sag = as.data.frame.matrix(table(sag$sampling.session, sag$survival))
colnames(tab.sag) = c('dead','alive')
tab.sag$prop.survival = tab.sag$alive/(tab.sag$dead+tab.sag$alive)
tab.sag$prop.survival = round(tab.sag$prop.survival,3)
tab.sag
#overall survival
l = length(subset(sag, survival == "1")$svl)
t = length(sag$svl)
g = l/t
sag.overall.survival = round(g,3)
#percentage survive (%)
sag.overall.survival

tab.sag <- subset(tab.sag, select = -c(dead, alive))
tab.sag = as.data.frame(tab.sag)
tab.sag[nrow(tab.sag) + c(1:2),] = c(sag.overall.survival)
rownames(tab.sag)[rownames(tab.sag) == "6"] <- "cumulative"
rownames(tab.sag)[rownames(tab.sag) == "7"] <- "cumulative.ran"
tab.sag

##10b. Anolis carolinensis----
#survival percentage through time
caro$sampling.session <- factor(caro$sampling.session, levels = 
                                  c("fall.2015", 
                                    "spring.2016", "fall.2016", 
                                    "spring.2017", "fall.2017"))
tab.caro = as.data.frame.matrix(table(caro$sampling.session, caro$survival))
colnames(tab.caro) = c('dead','alive')
tab.caro$prop.survival = tab.caro$alive/(tab.caro$dead+tab.caro$alive)
tab.caro$prop.survival = round(tab.caro$prop.survival,3)
tab.caro
#overall survival
l = length(subset(caro, survival == "1")$svl)
t = length(caro$svl)
g = l/t
caro.overall.survival = round(g,3)
#percentage survive (%)
caro.overall.survival

tab.caro <- subset(tab.caro, select = -c(dead, alive))
tab.caro = as.data.frame(tab.caro)
tab.caro[nrow(tab.caro) + c(1:2),] = c(caro.overall.survival)
rownames(tab.caro)[rownames(tab.caro) == "6"] <- "cumulative"
rownames(tab.caro)[rownames(tab.caro) == "7"] <- "cumulative.ran"
tab.caro

#11. LD1: Selection ----

##11a. Anolis sagrei----
library(lme4); library(nlme); library(mgcv)


###11a(i) Visualize surfaces ----

par(mfrow=c(1,2))
#LD1
#visualize selection surface using a cubic spline
sag.LD1.gam <- gam(survival ~ s(LD1), data = sag, 
                   family = binomial(link = "logit"), 
                   method = "GCV.Cp")

plot(sag.LD1.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.LD1.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "LD Axis 1", #main = "Cumulative: Anolis sagrei",
     col = "brown",las = 1,#xaxt="n", 
     ylim = c(0,0.4)
)
rug(sag$LD1, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$LD1, side=3, col = "brown", lwd=0.01)

#suggesting stabilizing selection

###11a(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = sag, 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#time as random effect (linear)
sag.glm.lin.ran <- glmer(survival ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                         data = sag, 
                         family = binomial(link = "logit"))

sag.glm.lin.ran.drop <- glmer(survival ~ #scale(LD1)+
                                scale(LD3) + 
                                (1|sampling.session), 
                              data = sag, 
                              family = binomial(link = "logit"))

s.sag.glm.lin.ran = anova(sag.glm.lin.ran, sag.glm.lin.ran.drop)
s.sag.glm.lin.ran

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = sag, 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#time as random effect (including nonlinear and correlational terms)
sag.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = sag, 
                          family = binomial(link = "logit"))
#model with term of interest dropped
sag.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                 #I(scale(LD1)^2)+
                                 I(scale(LD3)^2)+
                                 scale(LD1):scale(LD3)+
                                 (1|sampling.session), 
                               data = sag, 
                               family = binomial(link = "logit"))
#LLR of differences
s.sag.glm.full.ran = anova(sag.glm.full.ran, sag.glm.full.ran.drop)
s.sag.glm.full.ran

###11a(iii) LMs ----

#standard model (linear)
str(data)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = sag)
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#time as random effect (linear)
sag.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                        data = sag)
s.sag.ols.lin.ran = summary(sag.ols.lin.ran)
s.sag.ols.lin.ran

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = sag)
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full
#time as random effect (including nonlinear and correlational terms)
sag.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                           I(scale(LD1)^2)+I(scale(LD3)^2)+
                           scale(LD1):scale(LD3)+
                           (1|sampling.session), 
                         data = sag)

s.sag.ols.full.ran = summary(sag.ols.full.ran)
s.sag.ols.full.ran$coefficients


###11a(iv) Results table----
#let's combine all of this together
tab.sag
sagrei.table = tab.sag
sagrei.table$LD1.lin <- NA
sagrei.table$LD1.lin.SE <- NA
sagrei.table$LD1.lin.p <- NA
sagrei.table$LD1.quad <- NA
sagrei.table$LD1.quad.SE <- NA
sagrei.table$LD1.quad.p <- NA
sagrei.table

s.sag.ols.lin$coefficients
#add in the standard model results first
sagrei.table[6,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[6,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[6,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[6,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[6,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[6,7] <- round(s.sag.glm.full[3,3],3)
sagrei.table

#and now the random effects models
sagrei.table[7,2] <- round(s.sag.ols.lin.ran$coefficients[2,1],3)
sagrei.table[7,3] <- round(s.sag.ols.lin.ran$coefficients[2,2],3)
sagrei.table[7,4] <- round(s.sag.glm.lin.ran[2,8],3)
sagrei.table[7,5] <- round(s.sag.ols.full.ran$coefficients[4,1],3)
sagrei.table[7,6] <- round(s.sag.ols.full.ran$coefficients[4,2],3)
sagrei.table[7,7] <- round(s.sag.glm.full.ran[2,8],3)
sagrei.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2015"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[1,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[1,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[1,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[1,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[1,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[1,7] <- round(s.sag.glm.full[3,3],3)
sagrei.table

#spring 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[2,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[2,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[2,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[2,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[2,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[2,7] <- round(s.sag.glm.full[3,3],3)
sagrei.table

#fall 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[3,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[3,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[3,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[3,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[3,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[3,7] <- round(s.sag.glm.full[3,3],3)
sagrei.table

#spring 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[4,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[4,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[4,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[4,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[4,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[4,7] <- round(s.sag.glm.full[3,3],3)
sagrei.table

#fall 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[5,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[5,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[5,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[5,5] <- round(s.sag.ols.full$coefficients[4,1],3)
sagrei.table[5,6] <- round(s.sag.ols.full$coefficients[4,2],3)
sagrei.table[5,7] <- round(s.sag.glm.full[3,3],3)

###11a(v) LD1 final table----
sagrei.table$LD1.quad = sagrei.table$LD1.quad*2
sagrei.table$LD1.quad.SE = sagrei.table$LD1.quad.SE*2
sagrei.table
LD1.sag.table = sagrei.table
LD1.sag.table

##11b. Anolis carolinensis----


###11b(i) Visualize surfaces ----

par(mfrow=c(1,2))
#LD1
#visualize selection surface using a cubic spline
caro.LD1.gam <- gam(survival ~ s(LD1), data = caro, 
                    family = binomial(link = "logit"), 
                    method = "GCV.Cp")

plot(caro.LD1.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.LD1.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "LD Axis 1", #main = "Cumulative: Anolis caro",
     col = "forestgreen",las = 1,#xaxt="n", 
     ylim = c(0,0.6)
)
rug(caro$LD1, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, survival == "1")$LD1, side=3, col = "forestgreen", lwd=0.01)

#LD1 is funky

###11b(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = caro, 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#time as random effect (linear)
caro.glm.lin.ran <- glmer(survival ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                          data = caro, 
                          family = binomial(link = "logit"))
caro.glm.lin.ran.drop <- glmer(survival ~ #scale(LD1)+
                                 scale(LD3) + (1|sampling.session), 
                               data = caro, 
                               family = binomial(link = "logit"))
s.caro.glm.lin.ran = anova(caro.glm.lin.ran, caro.glm.lin.ran.drop)
s.caro.glm.lin.ran

#full model
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = caro, 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type ="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
#time as random effect (including nonlinear and correlational terms)
caro.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                             I(scale(LD1)^2)+I(scale(LD3)^2)+
                             scale(LD1):scale(LD3)+
                             (1|sampling.session), 
                           data = caro, 
                           family = binomial(link = "logit"))
#drop term of interest
caro.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                  #I(scale(LD1)^2)+
                                  I(scale(LD3)^2)+
                                  scale(LD1):scale(LD3)+
                                  (1|sampling.session), 
                                data = caro, 
                                family = binomial(link = "logit"))
s.caro.glm.full.ran = anova(caro.glm.full.ran, caro.glm.full.ran.drop)


###11b(iii) LMs ----

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = caro)
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#time as random effect (linear)
caro.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                         data = caro)
s.caro.ols.lin.ran = summary(caro.ols.lin.ran)
s.caro.ols.lin.ran

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = caro)
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full
#time as random effect (including nonlinear and correlational terms)
caro.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = caro)
s.caro.ols.full.ran = summary(caro.ols.full.ran)
s.caro.ols.full.ran$coefficients


###11b(iv) Results table----
#let's combine all of this together
tab.caro
caro.table = tab.caro
caro.table$LD1.lin <- NA
caro.table$LD1.lin.SE <- NA
caro.table$LD1.lin.p <- NA
caro.table$LD1.quad <- NA
caro.table$LD1.quad.SE <- NA
caro.table$LD1.quad.p <- NA
caro.table

s.caro.ols.lin$coefficients
#add in the standard model results first
caro.table[6,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[6,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[6,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[6,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[6,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[6,7] <- round(s.caro.glm.full[3,3],3)
caro.table

#and now the random effects models
caro.table[7,2] <- round(s.caro.ols.lin.ran$coefficients[2,1],3)
caro.table[7,3] <- round(s.caro.ols.lin.ran$coefficients[2,2],3)
caro.table[7,4] <- round(s.caro.glm.lin.ran[2,8],3)
caro.table[7,5] <- round(s.caro.ols.full.ran$coefficients[4,1],3)
caro.table[7,6] <- round(s.caro.ols.full.ran$coefficients[4,2],3)
caro.table[7,7] <- round(s.caro.glm.full.ran[2,8],3)
caro.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2015"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[1,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[1,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[1,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[1,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[1,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[1,7] <- round(s.caro.glm.full[3,3],3)
caro.table

#spring 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[2,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[2,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[2,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[2,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[2,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[2,7] <- round(s.caro.glm.full[3,3],3)
caro.table

#fall 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[3,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[3,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[3,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[3,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[3,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[3,7] <- round(s.caro.glm.full[3,3],3)
caro.table

#spring 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[4,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[4,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[4,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[4,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[4,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[4,7] <- round(s.caro.glm.full[3,3],3)
caro.table

#fall 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[5,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
caro.table[5,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
caro.table[5,4] <- round(s.caro.glm.lin[1,3],3)
caro.table[5,5] <- round(s.caro.ols.full$coefficients[4,1],3)
caro.table[5,6] <- round(s.caro.ols.full$coefficients[4,2],3)
caro.table[5,7] <- round(s.caro.glm.full[3,3],3)

###11b(v) LD1 final table----
caro.table$LD1.quad = caro.table$LD1.quad*2
caro.table$LD1.quad.SE = caro.table$LD1.quad.SE*2
caro.table
LD1.caro.table = caro.table


#12. LD3: Selection----

##12a. Anolis sagrei----
library(lme4); library(nlme); library(mgcv)


###12a(i) Visualize surfaces ----

par(mfrow=c(1,1))
#LD3
#visualize selection surface using a cubic spline
sag.LD3.gam <- gam(survival ~ s(LD3), data = sag, 
                   family = binomial(link = "logit"), 
                   method = "GCV.Cp")

plot(sag.LD3.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.LD3.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "LD Axis 3", #main = "Cumulative: Anolis sagrei",
     col = "brown",las = 1,#xaxt="n", 
     ylim = c(0,0.4)
)
rug(sag$LD3, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$LD3, side=3, col = "brown", lwd=0.01)

#ok
#suggesting stabilizing selection

###12a(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = sag, 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#time as random effect (linear)
sag.glm.lin.ran <- glmer(survival ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                         data = sag, 
                         family = binomial(link = "logit"))
sag.glm.lin.ran.drop <- glmer(survival ~ scale(LD1)+#scale(LD3) + 
                                (1|sampling.session), 
                              data = sag, 
                              family = binomial(link = "logit"))
s.sag.glm.lin.ran = anova(sag.glm.lin.ran, sag.glm.lin.ran.drop)
s.sag.glm.lin.ran



#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = sag, 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#time as random effect (including nonlinear and correlational terms)
sag.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = sag, 
                          family = binomial(link = "logit"))
#drop term of interest
sag.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                 I(scale(LD1)^2)+
                                 #I(scale(LD3)^2)+
                                 scale(LD1):scale(LD3)+
                                 (1|sampling.session), 
                               data = sag, 
                               family = binomial(link = "logit"))
s.sag.glm.full.ran= anova(sag.glm.full.ran, sag.glm.full.ran.drop)

###12a(iii) LMs ----

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = sag)
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#time as random effect (linear)
sag.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                        data = sag)
s.sag.ols.lin.ran = summary(sag.ols.lin.ran)
s.sag.ols.lin.ran

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = sag)
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full
#time as random effect (including nonlinear and correlational terms)
sag.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                           I(scale(LD1)^2)+I(scale(LD3)^2)+
                           scale(LD1):scale(LD3)+
                           (1|sampling.session), 
                         data = sag)
s.sag.ols.full.ran = summary(sag.ols.full.ran)
s.sag.ols.full.ran$coefficients

###12a(iv) Results table----
#let's combine all of this together
tab.sag
sagrei.table = tab.sag
sagrei.table$LD3.lin <- NA
sagrei.table$LD3.lin.SE <- NA
sagrei.table$LD3.lin.p <- NA
sagrei.table$LD3.quad <- NA
sagrei.table$LD3.quad.SE <- NA
sagrei.table$LD3.quad.p <- NA
sagrei.table

#add in the standard model results first
sagrei.table[6,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[6,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[6,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[6,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[6,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[6,7] <- round(s.sag.glm.full[4,3],3)
sagrei.table

#and now the random effects models
sagrei.table[7,2] <- round(s.sag.ols.lin.ran$coefficients[3,1],3)
sagrei.table[7,3] <- round(s.sag.ols.lin.ran$coefficients[3,2],3)
sagrei.table[7,4] <- round(s.sag.glm.lin.ran[2,8],3)
sagrei.table[7,5] <- round(s.sag.ols.full.ran$coefficients[5,1],3)
sagrei.table[7,6] <- round(s.sag.ols.full.ran$coefficients[5,2],3)
sagrei.table[7,7] <- round(s.sag.glm.full.ran[2,8],3)
sagrei.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2015"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[1,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[1,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[1,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[1,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[1,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[1,7] <- round(s.sag.glm.full[4,3],3)
sagrei.table

#spring 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[2,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[2,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[2,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[2,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[2,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[2,7] <- round(s.sag.glm.full[4,3],3)
sagrei.table

#fall 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[3,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[3,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[3,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[3,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[3,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[3,7] <- round(s.sag.glm.full[4,3],3)
sagrei.table

#spring 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[4,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[4,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[4,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[4,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[4,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[4,7] <- round(s.sag.glm.full[4,3],3)
sagrei.table

#fall 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                  data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[5,2] <- round(s.sag.ols.lin$coefficients[3,1],3)
sagrei.table[5,3] <- round(s.sag.ols.lin$coefficients[3,2],3)
sagrei.table[5,4] <- round(s.sag.glm.lin[2,3],3)
sagrei.table[5,5] <- round(s.sag.ols.full$coefficients[5,1],3)
sagrei.table[5,6] <- round(s.sag.ols.full$coefficients[5,2],3)
sagrei.table[5,7] <- round(s.sag.glm.full[4,3],3)

###12a(v) LD3 final table----
sagrei.table
#we need to double the quadratic coefficients (and associated SE's)
sagrei.table$LD3.quad = sagrei.table$LD3.quad*2
sagrei.table$LD3.quad.SE = sagrei.table$LD3.quad.SE*2
sagrei.table
LD3.sag.table = sagrei.table

##12b. Anolis carolinensis----
library(lme4); library(nlme); library(mgcv)


###12b(i) Visualize surfaces ----

par(mfrow=c(1,1))
#LD3
#visualize selection surface using a cubic spline
caro.LD3.gam <- gam(survival ~ s(LD3), data = caro, 
                    family = binomial(link = "logit"), 
                    method = "GCV.Cp")

plot(caro.LD3.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.LD3.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "LD Axis 3", #main = "Cumulative: Anolis caro",
     col = "forestgreen",las = 1,#xaxt="n", 
     ylim = c(0,0.4)
)
rug(caro$LD3, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, survival == "1")$LD3, side=3, col = "forestgreen", lwd=0.01)

#ok
#suggesting stabilizing selection

###12b(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = caro, 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#time as random effect (linear)
caro.glm.lin.ran <- glmer(survival ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                          data = caro, 
                          family = binomial(link = "logit"))
caro.glm.lin.ran.drop <- glmer(survival ~ scale(LD1)+#scale(LD3) + 
                                 (1|sampling.session), 
                               data = caro, 
                               family = binomial(link = "logit"))
s.caro.glm.lin.ran = anova(caro.glm.lin.ran, caro.glm.lin.ran.drop)
s.caro.glm.lin.ran

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = caro, 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#time as random effect (including nonlinear and correlational terms)
caro.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                             I(scale(LD1)^2)+I(scale(LD3)^2)+
                             scale(LD1):scale(LD3)+
                             (1|sampling.session), 
                           data = caro, 
                           family = binomial(link = "logit"))
#drop term of interest
caro.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                  I(scale(LD1)^2)+
                                  #I(scale(LD3)^2)+
                                  scale(LD1):scale(LD3)+
                                  (1|sampling.session), 
                                data = caro, 
                                family = binomial(link = "logit"))
#anova
s.caro.glm.full.ran = anova(caro.glm.full.ran, caro.glm.full.ran.drop)


###12b(iii) LMs ----

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = caro)
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#time as random effect (linear)
caro.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(LD1)+scale(LD3) + (1|sampling.session), 
                         data = caro)
s.caro.ols.lin.ran = summary(caro.ols.lin.ran)
s.caro.ols.lin.ran

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = caro)
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full
#time as random effect (including nonlinear and correlational terms)
caro.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = caro)
s.caro.ols.full.ran = summary(caro.ols.full.ran)
s.caro.ols.full.ran$coefficients


###12b(iv) Results table----
#let's combine all of this together
tab.caro
caro.table = tab.caro
caro.table$LD3.lin <- NA
caro.table$LD3.lin.SE <- NA
caro.table$LD3.lin.p <- NA
caro.table$LD3.quad <- NA
caro.table$LD3.quad.SE <- NA
caro.table$LD3.quad.p <- NA
caro.table

#add in the standard model results first
caro.table[6,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[6,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[6,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[6,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[6,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[6,7] <- round(s.caro.glm.full[4,3],3)
caro.table

#and now the random effects models
caro.table[7,2] <- round(s.caro.ols.lin.ran$coefficients[3,1],3)
caro.table[7,3] <- round(s.caro.ols.lin.ran$coefficients[3,2],3)
caro.table[7,4] <- round(s.caro.glm.lin.ran[2,8],3)
caro.table[7,5] <- round(s.caro.ols.full.ran$coefficients[5,1],3)
caro.table[7,6] <- round(s.caro.ols.full.ran$coefficients[5,2],3)
caro.table[7,7] <- round(s.caro.glm.full.ran[2,8],3)
caro.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2015"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[1,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[1,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[1,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[1,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[1,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[1,7] <- round(s.caro.glm.full[4,3],3)
caro.table

#spring 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[2,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[2,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[2,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[2,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[2,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[2,7] <- round(s.caro.glm.full[4,3],3)
caro.table

#fall 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[3,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[3,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[3,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[3,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[3,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[3,7] <- round(s.caro.glm.full[4,3],3)
caro.table

#spring 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[4,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[4,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[4,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[4,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[4,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[4,7] <- round(s.caro.glm.full[4,3],3)
caro.table

#fall 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(LD1)+scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3), 
                   data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[5,2] <- round(s.caro.ols.lin$coefficients[3,1],3)
caro.table[5,3] <- round(s.caro.ols.lin$coefficients[3,2],3)
caro.table[5,4] <- round(s.caro.glm.lin[2,3],3)
caro.table[5,5] <- round(s.caro.ols.full$coefficients[5,1],3)
caro.table[5,6] <- round(s.caro.ols.full$coefficients[5,2],3)
caro.table[5,7] <- round(s.caro.glm.full[4,3],3)

###12b(v) LD3 final table----

caro.table
#we need to double the quadratic coefficients (and associated SE's)
caro.table$LD3.quad = caro.table$LD3.quad*2
caro.table$LD3.quad.SE = caro.table$LD3.quad.SE*2
caro.table
LD3.caro.table = caro.table

#13. Correlational selection ----

##13a. Anolis sagrei----
###13a(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = sag, 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#time as random effect (including nonlinear and correlational terms)
sag.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = sag, 
                          family = binomial(link = "logit"))

sag.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                 I(scale(LD1)^2)+I(scale(LD3)^2)+
                                 #scale(LD1):scale(LD3)+
                                 (1|sampling.session), 
                               data = sag, 
                               family = binomial(link = "logit"))
su.sag.ss.ran.glm = anova(sag.glm.full.ran,sag.glm.full.ran.drop)
su.sag.ss.ran.glm


###13a(iii) LMs ----

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = sag)
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full
#time as random effect (including nonlinear and correlational terms)
sag.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                           I(scale(LD1)^2)+I(scale(LD3)^2)+
                           scale(LD1):scale(LD3)+
                           (1|sampling.session), 
                         data = sag)
s.sag.ols.full.ran = summary(sag.ols.full.ran)
s.sag.ols.full.ran$coefficients

###13a(iv) Results table----
#let's combine all of this together
tab.sag
sagrei.table = tab.sag
sagrei.table$LD1xLD3 <- NA
sagrei.table$LD1xLD3.SE <- NA
sagrei.table$LD1xLD3.p <- NA
sagrei.table

#add in the standard model results first
sagrei.table[6,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[6,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[6,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

#and now the random effects models
sagrei.table[7,2] <- round(s.sag.ols.full.ran$coefficients[6,1],3)
sagrei.table[7,3] <- round(s.sag.ols.full.ran$coefficients[6,2],3)
sagrei.table[7,4] <- round(s.sag.glm.full.ran[2,8],3)
sagrei.table

#and now let's start the individual sampling sessions

#fall 2015
#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[1,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[1,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[1,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

#spring.2016
#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[2,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[2,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[2,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

#fall.2016
#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[3,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[3,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[3,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

#spring.2017
#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[4,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[4,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[4,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

#fall.2017
#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(sag, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                     I(scale(LD1)^2)+I(scale(LD3)^2)+
                     scale(LD1):scale(LD3), 
                   data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[5,2] <- round(s.sag.ols.full$coefficients[6,1],3)
sagrei.table[5,3] <- round(s.sag.ols.full$coefficients[6,2],3)
sagrei.table[5,4] <- round(s.sag.glm.full[5,3],3)
sagrei.table

###13a(v) LD3 final table----
sagrei.table
LDcorr.sag.table = sagrei.table

##13b. Anolis carolinensis----
###13b(ii) GLMs ----

# for selection analyses:
# we will get the p-values from binomial GLMs
# and then estimate the strength of the selection surfaces
# from OLS models of relative survival

#first, let's get the p-values
#linear terms must be estimated independently from
#models that include nonlinear and correlational terms

#and i want to run 2x models
#1. stanard mixed models
#2. the same models with time as a random effect

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = caro, 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#time as random effect (including nonlinear and correlational terms)
caro.glm.full.ran <- glmer(survival ~ scale(LD1)+scale(LD3)+
                             I(scale(LD1)^2)+I(scale(LD3)^2)+
                             scale(LD1):scale(LD3)+
                             (1|sampling.session), 
                           data = caro, 
                           family = binomial(link = "logit"))
caro.glm.full.ran.drop <- glmer(survival ~ scale(LD1)+scale(LD3)+
                                  I(scale(LD1)^2)+I(scale(LD3)^2)+
                                  #scale(LD1):scale(LD3)+
                                  (1|sampling.session), 
                                data = caro, 
                                family = binomial(link = "logit"))
s.caro.glm.full.ran = anova(caro.glm.full.ran, caro.glm.full.ran.drop)
s.caro.glm.full.ran


###13b(iii) LMs ----

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = caro)
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full
#time as random effect (including nonlinear and correlational terms)
caro.ols.full.ran <- lmer(survival ~ scale(LD1)+scale(LD3)+
                            I(scale(LD1)^2)+I(scale(LD3)^2)+
                            scale(LD1):scale(LD3)+
                            (1|sampling.session), 
                          data = caro)
s.caro.ols.full.ran = summary(caro.ols.full.ran)
s.caro.ols.full.ran$coefficients

###13b(iv) Results table----
#let's combine all of this together
tab.caro
caro.table = tab.caro
caro.table$LD1xLD3 <- NA
caro.table$LD1xLD3.SE <- NA
caro.table$LD1xLD3.p <- NA
caro.table

#add in the standard model results first
caro.table[6,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[6,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[6,4] <- round(s.caro.glm.full[5,3],3)
caro.table

#and now the random effects models
caro.table[7,2] <- round(s.caro.ols.full.ran$coefficients[6,1],3)
caro.table[7,3] <- round(s.caro.ols.full.ran$coefficients[6,2],3)
caro.table[7,4] <- round(s.caro.glm.full.ran[2,8],3)
caro.table

#and now let's start the individual sampling sessions

#fall 2015
#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2015"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[1,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[1,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[1,4] <- round(s.caro.glm.full[5,3],3)
caro.table

#spring.2016
#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[2,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[2,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[2,4] <- round(s.caro.glm.full[5,3],3)
caro.table

#fall.2016
#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[3,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[3,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[3,4] <- round(s.caro.glm.full[5,3],3)
caro.table

#spring.2017
#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "spring.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[4,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[4,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[4,4] <- round(s.caro.glm.full[5,3],3)
caro.table

#fall.2017
#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(LD1)+scale(LD3)+
                       I(scale(LD1)^2)+I(scale(LD3)^2)+
                       scale(LD1):scale(LD3), 
                     data = subset(caro, sampling.session == "fall.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(LD1)+scale(LD3)+
                      I(scale(LD1)^2)+I(scale(LD3)^2)+
                      scale(LD1):scale(LD3), 
                    data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
caro.table[5,2] <- round(s.caro.ols.full$coefficients[6,1],3)
caro.table[5,3] <- round(s.caro.ols.full$coefficients[6,2],3)
caro.table[5,4] <- round(s.caro.glm.full[5,3],3)
caro.table

###13b(v) LD3 final table----
caro.table
LDcorr.caro.table = caro.table

#14. LD selection tables----

##14a. LD1 sagrei----
LD1.sag.table
write.table(LD1.sag.table, 
            "clipboard", sep="\t", row.names=TRUE)
##14b. LD3 sagrei----
LD3.sag.table
write.table(LD3.sag.table, 
            "clipboard", sep="\t", row.names=TRUE)
##14c. LD1 carolinensis----
LD1.caro.table
write.table(LD1.caro.table, 
            "clipboard", sep="\t", row.names=TRUE)
##14d. LD3 carolinensis----
LD3.caro.table
write.table(LD3.caro.table, 
            "clipboard", sep="\t", row.names=TRUE)
##14e. Corr sagrei----
LDcorr.sag.table
write.table(LDcorr.sag.table, 
            "clipboard", sep="\t", row.names=TRUE)
##14f. Corr carolinensis----
LDcorr.caro.table
write.table(LDcorr.caro.table, 
            "clipboard", sep="\t", row.names=TRUE)


#15. Euclidean selection----

#now we will test for stabilizing selection in a different way
#we will identify the mean phenotype
#and then calculate the euclidean distance of each individual to that phenotype
#and estiamte survival as a (linear) function of euclidean distance~mean

##15a. Anolis sagrei----

###15a(i). Calculate Euclidean distances----
#stabilizing selection to a centroid
sag.points <- data.frame(x=sag$LD1, y = sag$LD3)
sag.cnt = c(mean(sag.points[,1]),mean(sag.points[,2]))
sag$euclid = apply(sag.points,1,function(x,sag.cnt) {
  (sqrt((x[1] - sag.cnt[1])^2+(x[2]-sag.cnt[2])^2))
},sag.cnt)
sag$scale.euclid = scale(sag$euclid)

#plot to make sure it looks like the code is working
#mean phenotype = red dot
#most peripheral phenotypes should be blue
#phenotypes near the mean should be yellow
ggplot(sag, aes(x=LD1, y=LD3)) + 
  geom_point(aes(fill=euclid), pch=21, size=5) +
  scale_fill_gradient(low = "yellow", high = "blue", na.value = NA)+
  geom_point(aes(mean(sag.points[,1]),mean(sag.points[,2])), color = "red", size = 8)+
  theme_bw()

###15a(ii). Visualize selection surface----

#visualize selection surface using a cubic spline
par(mfrow=c(1,1))
sag.ss.gam <- gam(survival ~ s(euclid), data = sag, 
                  family = binomial(link = "logit"), 
                  method = "GCV.Cp")

plot(sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance from mean phenotype", #main = "Cumulative: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.05,0.35)
)
rug(sag$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#this is cool and really strong, so let's see how it compares to all other seasons
par(mfrow=c(2,3))
plot(sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Cumulative: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6)
)
rug(sag$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2015
f15.sag.ss.gam <- gam(survival ~ s(euclid), data = subset(sag, sampling.session == "fall.2015"), 
                      family = binomial(link = "logit"), 
                      method = "GCV.Cp")
plot(f15.sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f15.sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2015: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(sag$euclid), max(sag$euclid))
)
rug(subset(sag, sampling.session == "fall.2015")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session == "fall.2015" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#spring.2016
s16.sag.ss.gam <- gam(survival ~ s(euclid), data = subset(sag, sampling.session == "spring.2016"), 
                      family = binomial(link = "logit"), 
                      method = "GCV.Cp")
plot(s16.sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(s16.sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Spring 2016: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(sag$euclid), max(sag$euclid))
)
rug(subset(sag, sampling.session == "spring.2016")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session == "spring.2016" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2016
f16.sag.ss.gam <- gam(survival ~ s(euclid), data = subset(sag, sampling.session == "fall.2016"), 
                      family = binomial(link = "logit"), 
                      method = "GCV.Cp")
plot(f16.sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f16.sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2016: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(sag$euclid), max(sag$euclid))
)
rug(subset(sag, sampling.session == "fall.2016")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session == "fall.2016" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#spring.2017
s17.sag.ss.gam <- gam(survival ~ s(euclid), data = subset(sag, sampling.session == "spring.2017"), 
                      family = binomial(link = "logit"), 
                      method = "GCV.Cp")
plot(s17.sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(s17.sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Spring 2017: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(sag$euclid), max(sag$euclid))
)
rug(subset(sag, sampling.session == "spring.2017")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session == "spring.2017" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2017
f17.sag.ss.gam <- gam(survival ~ s(euclid), data = subset(sag, sampling.session == "fall.2017"), 
                      family = binomial(link = "logit"), 
                      method = "GCV.Cp")
plot(f17.sag.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f17.sag.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2017: Anolis sagrei",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(sag$euclid), max(sag$euclid))
)
rug(subset(sag, sampling.session == "fall.2017")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session == "fall.2017" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#this is pretty cool
#definitely variation in the strength of stabilizing selection

###15a(iii). Estimate selection----

#make blank table
sag.table=data.frame(time=c("fall.2015",
                            "spring.2016", "fall.2016",
                            "spring.2017", "fall.2017",
                            "cumulative",
                            "cumulative*"))
sag.table
sag.table$n <- NA
sag.table$Beta <- NA
sag.table$S.E. <- NA
sag.table$z <- NA
sag.table$p <- NA
sag.table

#sample sizes
sag.table[1,2] <- length(subset(sag, sampling.session == "fall.2015")$sampling.session)
sag.table[2,2] <- length(subset(sag, sampling.session == "spring.2016")$sampling.session)
sag.table[3,2] <- length(subset(sag, sampling.session == "fall.2016")$sampling.session)
sag.table[4,2] <- length(subset(sag, sampling.session == "spring.2017")$sampling.session)
sag.table[5,2] <- length(subset(sag, sampling.session == "fall.2017")$sampling.session)
sag.table[6:7,2] <- length(sag$sampling.session)
sag.table

# test for stabilizing selection
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = sag, 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = sag)
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

# include sampling period as a random effect
sag.ss.ran.glm <- glmer(survival ~ scale(euclid) + 
                          (1|sampling.session), 
                        data = sag, 
                        family = binomial(link = "logit"))
sag.ss.ran.glm.drop <- glmer(survival ~ 1 + 
                               (1|sampling.session), 
                             data = sag, 
                             family = binomial(link = "logit"))
su.sag.ss.ran.glm = anova(sag.ss.ran.glm,sag.ss.ran.glm.drop)
su.sag.ss.ran.glm

#OLS with random effect
sag.ss.ran.ols <- lmer(survival/mean(survival) ~ scale(euclid) + (1|sampling.session), 
                       data = sag)
s.sag.ss.ran.ols <- summary(sag.ss.ran.ols)

#insert simple model results first (no random effects)
sag.table[6,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[6,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[6,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[6,6] <- round(s.sag.ss.simple.glm[1,3],3)

#now random effects models
sag.table[7,3] <- round(s.sag.ss.ran.ols$coefficients[2,1],3)
sag.table[7,4] <- round(s.sag.ss.ran.ols$coefficients[2,2],3)
sag.table[7,5] <- round(su.sag.ss.ran.glm[2,6],3)
sag.table[7,6] <- round(su.sag.ss.ran.glm[2,8],3)
su.sag.ss.ran.glm
sag.table

#fall.2015
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2015"), 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2015"))
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

#insert simple model results first (no random effects)
sag.table[1,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[1,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[1,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[1,6] <- round(s.sag.ss.simple.glm[1,3],3)
sag.table

#spring.2016
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = subset(sag, sampling.session == "spring.2016"), 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = subset(sag, sampling.session == "spring.2016"))
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

#insert simple model results first (no random effects)
sag.table[2,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[2,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[2,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[2,6] <- round(s.sag.ss.simple.glm[1,3],3)
sag.table

#fall.2016
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2016"), 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2016"))
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

#insert simple model results first (no random effects)
sag.table[3,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[3,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[3,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[3,6] <- round(s.sag.ss.simple.glm[1,3],3)
sag.table

#spring.2017
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = subset(sag, sampling.session == "spring.2017"), 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = subset(sag, sampling.session == "spring.2017"))
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

#insert simple model results first (no random effects)
sag.table[4,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[4,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[4,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[4,6] <- round(s.sag.ss.simple.glm[1,3],3)
sag.table

#fall.2017
#simple GLM
sag.ss.simple.glm <- glm(survival ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2017"), 
                         family = binomial(link = "logit"))
s.sag.ss.simple.glm = Anova(sag.ss.simple.glm, type="III")
su.sag.ss.simple.glm = summary(sag.ss.simple.glm)

#OLS
sag.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                         data = subset(sag, sampling.session == "fall.2017"))
s.sag.ss.simple.ols <- summary(sag.ss.simple.ols)

#insert simple model results first (no random effects)
sag.table[5,3] <- round(s.sag.ss.simple.ols$coefficients[2,1],3)
sag.table[5,4] <- round(s.sag.ss.simple.ols$coefficients[2,2],3)
sag.table[5,5] <- round(su.sag.ss.simple.glm$coefficients[2,3],3)
sag.table[5,6] <- round(s.sag.ss.simple.glm[1,3],3)
sag.table

###15a(iv). Sagrei stabilizing selection----
sag.table
sagrei.euclidean.selection <- sag.table
sagrei.euclidean.selection

##15b. Anolis carolinensis----

###15b(i). Calculate Euclidean distances----
#stabilizing selection to a centroid
caro.points <- data.frame(x=caro$LD1, y = caro$LD3)
caro.cnt = c(mean(caro.points[,1]),mean(caro.points[,2]))
caro$euclid = apply(caro.points,1,function(x,caro.cnt) {
  (sqrt((x[1] - caro.cnt[1])^2+(x[2]-caro.cnt[2])^2))
},caro.cnt)
caro$scale.euclid = scale(caro$euclid)

#plot to make sure it looks like the code is working
#mean phenotype = red dot
#most peripheral phenotypes should be blue
#phenotypes near the mean should be yellow
ggplot(caro, aes(x=LD1, y=LD3)) + 
  geom_point(aes(fill=euclid), pch=21, size=5) +
  scale_fill_gradient(low = "yellow", high = "blue", na.value = NA)+
  geom_point(aes(mean(caro.points[,1]),mean(caro.points[,2])), color = "red", size = 8)+
  theme_bw()

###15b(ii). Visualize selection surface----

#visualize selection surface using a cubic spline
par(mfrow=c(1,1))
caro.ss.gam <- gam(survival ~ s(euclid), data = caro, 
                   family = binomial(link = "logit"), 
                   method = "GCV.Cp")

plot(caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance from mean phenotype",# main = "Cumulative: Anolis carolinensis",
     col = "forestgreen",las = 1,xaxt="n"#, ylim = c(0.0,0.5)
)
rug(caro$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#this is cool.
#pretty weird.
#so let's see how it compares to all other seasons

par(mfrow=c(2,3))
plot(caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Cumulative: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6)
)
rug(caro$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2015
f15.caro.ss.gam <- gam(survival ~ s(euclid), data = subset(caro, sampling.session == "fall.2015"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp")
plot(f15.caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f15.caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2015: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(caro$euclid), max(caro$euclid))
)
rug(subset(caro, sampling.session == "fall.2015")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session == "fall.2015" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#spring.2016
s16.caro.ss.gam <- gam(survival ~ s(euclid), data = subset(caro, sampling.session == "spring.2016"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp")
plot(s16.caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(s16.caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Spring 2016: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(caro$euclid), max(caro$euclid))
)
rug(subset(caro, sampling.session == "spring.2016")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session == "spring.2016" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2016
f16.caro.ss.gam <- gam(survival ~ s(euclid), data = subset(caro, sampling.session == "fall.2016"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp")
plot(f16.caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f16.caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2016: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(caro$euclid), max(caro$euclid))
)
rug(subset(caro, sampling.session == "fall.2016")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session == "fall.2016" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#spring.2017
s17.caro.ss.gam <- gam(survival ~ s(euclid), data = subset(caro, sampling.session == "spring.2017"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp")
plot(s17.caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(s17.caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Spring 2017: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(caro$euclid), max(caro$euclid))
)
rug(subset(caro, sampling.session == "spring.2017")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session == "spring.2017" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#fall.2017
f17.caro.ss.gam <- gam(survival ~ s(euclid), data = subset(caro, sampling.session == "fall.2017"), 
                       family = binomial(link = "logit"), 
                       method = "GCV.Cp")
plot(f17.caro.ss.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(f17.caro.ss.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Euclidean distance to overall mean phenotype", 
     main = "Fall 2017: Anolis carolinensis",
     col = "brown",las = 1,xaxt="n", ylim = c(0.0,0.6),
     xlim = c(min(caro$euclid), max(caro$euclid))
)
rug(subset(caro, sampling.session == "fall.2017")$euclid, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session == "fall.2017" & survival == "1")$euclid, side=3, col = "brown", lwd=0.01)

#all sorts of madness going on there

###15b(iii). Estimate selection----

#make blank table
caro.table=data.frame(time=c("fall.2015",
                             "spring.2016", "fall.2016",
                             "spring.2017", "fall.2017",
                             "cumulative",
                             "cumulative*"))
caro.table
caro.table$n <- NA
caro.table$Beta <- NA
caro.table$S.E. <- NA
caro.table$z <- NA
caro.table$p <- NA
caro.table

#sample sizes
caro.table[1,2] <- length(subset(caro, sampling.session == "fall.2015")$sampling.session)
caro.table[2,2] <- length(subset(caro, sampling.session == "spring.2016")$sampling.session)
caro.table[3,2] <- length(subset(caro, sampling.session == "fall.2016")$sampling.session)
caro.table[4,2] <- length(subset(caro, sampling.session == "spring.2017")$sampling.session)
caro.table[5,2] <- length(subset(caro, sampling.session == "fall.2017")$sampling.session)
caro.table[6:7,2] <- length(caro$sampling.session)
caro.table

# test for stabilizing selection
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = caro, 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)
su.caro.ss.simple.glm
#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = caro)
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#time as random effect (linear)
caro.ss.ran.glm <- glmer(survival ~ scale(euclid) + 
                           (1|sampling.session), 
                         data = caro, 
                         family = binomial(link = "logit"))
caro.ss.ran.glm.drop <- glmer(survival ~ 1 + 
                                (1|sampling.session), 
                              data = caro, 
                              family = binomial(link = "logit"))
s.caro.glm.lin.ran = anova(caro.ss.ran.glm,caro.ss.ran.glm.drop)
s.caro.glm.lin.ran

#OLS with random effect
caro.ss.ran.ols <- lmer(survival/mean(survival) ~ scale(euclid) + (1|sampling.session), 
                        data = caro)
s.caro.ss.ran.ols <- summary(caro.ss.ran.ols)

#insert simple model results first (no random effects)
caro.table[6,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[6,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[6,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[6,6] <- round(s.caro.ss.simple.glm[1,3],3)

#now random effects models
caro.table[7,3] <- round(s.caro.ss.ran.ols$coefficients[2,1],3)
caro.table[7,4] <- round(s.caro.ss.ran.ols$coefficients[2,2],3)
caro.table[7,5] <- round(s.caro.glm.lin.ran[2,6],3)
caro.table[7,6] <- round(s.caro.glm.lin.ran[2,8],3)
caro.table

#fall.2015
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2015"), 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)

#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2015"))
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#insert simple model results first (no random effects)
caro.table[1,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[1,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[1,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[1,6] <- round(s.caro.ss.simple.glm[1,3],3)
caro.table

#spring.2016
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = subset(caro, sampling.session == "spring.2016"), 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)

#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = subset(caro, sampling.session == "spring.2016"))
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#insert simple model results first (no random effects)
caro.table[2,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[2,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[2,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[2,6] <- round(s.caro.ss.simple.glm[1,3],3)
caro.table

#fall.2016
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2016"), 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)

#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2016"))
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#insert simple model results first (no random effects)
caro.table[3,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[3,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[3,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[3,6] <- round(s.caro.ss.simple.glm[1,3],3)
caro.table

#spring.2017
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = subset(caro, sampling.session == "spring.2017"), 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)

#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = subset(caro, sampling.session == "spring.2017"))
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#insert simple model results first (no random effects)
caro.table[4,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[4,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[4,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[4,6] <- round(s.caro.ss.simple.glm[1,3],3)
caro.table

#fall.2017
#simple GLM
caro.ss.simple.glm <- glm(survival ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2017"), 
                          family = binomial(link = "logit"))
s.caro.ss.simple.glm = Anova(caro.ss.simple.glm, type="III")
su.caro.ss.simple.glm = summary(caro.ss.simple.glm)

#OLS
caro.ss.simple.ols <- glm(survival/mean(survival) ~ scale(euclid), 
                          data = subset(caro, sampling.session == "fall.2017"))
s.caro.ss.simple.ols <- summary(caro.ss.simple.ols)

#insert simple model results first (no random effects)
caro.table[5,3] <- round(s.caro.ss.simple.ols$coefficients[2,1],3)
caro.table[5,4] <- round(s.caro.ss.simple.ols$coefficients[2,2],3)
caro.table[5,5] <- round(su.caro.ss.simple.glm$coefficients[2,3],3)
caro.table[5,6] <- round(s.caro.ss.simple.glm[1,3],3)
caro.table

###15b(iv). carolinensis stabilizing selection----
caro.table
carolinensis.euclidean.selection <- caro.table
carolinensis.euclidean.selection

#15c. Euclidean tables ----

sagrei.euclidean.selection
carolinensis.euclidean.selection

#woohoo! this is really cool
#strong evidence of stabilizing selection towards a mean phenotype
#but at the micro-temporal level this selection fluctuates in strength


#16. Major axis transect ----

## 16a. carolinensis----
###16a(i) Transect plot----
#pdf(file = "carolinensis_majoraxis.pdf")

par(mfrow = c(1,3))
#2D surface plot
image.plot(caro.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="identify peak and trough")
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

caro$rel.fit = caro.fit$fitted.values
caro.fit.max = caro[which.max(caro$rel.fit),]
points(caro.fit.max$LD1, caro.fit.max$LD3, pch=21, bg="red", cex=2)
caro.fit.min = caro[which.min(caro$rel.fit),]
points(caro.fit.min$LD1, caro.fit.min$LD3, pch=21, bg="blue", cex=2)

#do things
image.plot(caro.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="draw transect line")
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

max.centroid = rbind(caro.fit.max$LD3, caro.fit.max$LD1)
min.centroid = rbind(caro.fit.min$LD3, caro.fit.min$LD1)
f = max.centroid[2,1]
g = max.centroid[1,1]
h = min.centroid[2,1]
j = min.centroid[1,1]

x <- caro$LD1
y <- caro$LD3

ys = c(f, g)
xs = c(h, j)
d = rbind(ys, xs)
d = as.data.frame(d)
d
line <- lm(V2~V1, d)

#draw an abline connecting the two species centroids
abline(line, col = 'red', lwd=2, lty=2)
#plot the two centroids
points(caro.fit.max$LD1, caro.fit.max$LD3, pch=21, bg="red", cex=2)
points(caro.fit.min$LD1, caro.fit.min$LD3, pch=21, bg="blue", cex=2)
#check data match
points(d[,1], d[,2], bg = "black", pch = 21, cex=2)
#looks good

# function draws line between two random points
drawline <- function(V2, V1){
  line <- lm(V2~V1, d)
  abline(line, col = "red", lwd=2, lty=2)
}

# function to measure distance from 'measure_point' to line
dist_to_line <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 

# function to measure distance to 'ref_point'
dist_to_ref <- function(r,t) {
  u <- (r[1]-t[1])^2
  v <- (r[2]-t[2])^2
  j <- sqrt(u+v)
}

line
m = summary(line)
m$coefficients
aa = m$coefficients[1,1]
bb = m$coefficients[2,1]
cc = min(caro$LD1)-2
f = aa+bb*cc
ref_point = c(cc,f)
line_point = c(d[2,1],d[2,2])
line_point


image.plot(caro.out, xlab = "LD1: Body size and foot length", 
           ylab = "LD3: Limb dimensions", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="set reference point")
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

#draw reference point for measuring
points(line_point[1],line_point[2], bg = "blue", pch = 21, cex=2)
points(ref_point[1],ref_point[2], bg = "lightblue", pch = 21, cex=2)
caro_point = c(d[1,1],d[1,2])
drawline(c(ref_point[1],line_point[1]),c(ref_point[2],line_point[2]))
#ok getting there

measure_point <- 0
i = 1
a <- NA

length = length(caro$svl)

image.plot(caro.out, xlab = "LD1: Body size and foot length", 
           ylab = "LD3: Limb dimensions", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="plot all individuals")
contour(caro.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

#visualize position of all individuals relative to the transect
for(i in 1:length){
  measure_point <- c(x[i],y[i])
  points(measure_point[1],measure_point[2], bg = "gray75", pch = 21)
  d <- dist_to_line(measure_point, ref_point, line_point)
  j <- dist_to_ref(measure_point,ref_point)
  dist <- sqrt((j^2 - d^2))
  a[i]<- dist
}
abline(line, col = 'red', lwd=2)

data_plotting <- data.frame(a,1)
#plot transect data
plot(data_plotting[,1],data_plotting[,2], 
     xlim = c(min(data_plotting[,1]-2),
              max(data_plotting[,1])),
     yaxt='n', xaxt = 'n',
     ylab="",
     xlab = "position along major axis transect")
abline(h=1, col="red")

#combine major.axis.transect data to main dataframe
caro$major.axis.transect <- a

caro.gam <- gam(survival ~ s(major.axis.transect), data = caro, 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", main = "",
     col = "forestgreen",las = 1,xaxt="n", ylim = c(0,0.4)
)
rug(caro$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#models
caro.glm <- glmer(survival ~ major.axis.transect+I(major.axis.transect^2) + (1|sampling.session), 
                  data = caro, 
                  family = binomial(link = "logit"))
caro.glm = Anova(caro.glm)
caro.glm

caro.ols <- lmer(survival/mean(survival) ~ major.axis.transect+I(major.axis.transect^2) + (1|sampling.session), 
                 data = caro)
summary(caro.ols)

###16.a(i) Surface plots----

par(mfrow=c(2,3))
image.plot(caro.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Cumulative")
contour(caro.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)


#surfaces through time
caro.f15 = droplevels(subset(caro, sampling.session == "fall.2015"))
x = cbind(caro.f15$LD1, caro.f15$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f15$survival)
caro.15.fit = Tps(x2, survival)
caro.f15.out = predictSurface(caro.15.fit)
image.plot(caro.f15.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2015")
contour(caro.f15.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)

#0.35
caro.s16 = droplevels(subset(caro, sampling.session == "spring.2016"))
x = cbind(caro.s16$LD1, caro.s16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.s16$survival)
fit = Tps(x2, survival)
caro.s16.out = predictSurface(fit)
image.plot(caro.s16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Spring 2016")
contour(caro.s16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)

#0.35
caro.f16 = droplevels(subset(caro, sampling.session == "fall.2016"))
x = cbind(caro.f16$LD1, caro.f16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f16$survival)
fit = Tps(x2, survival)
caro.f16.out = predictSurface(fit)
image.plot(caro.f16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2016")
contour(caro.f16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)
#0.6

caro.s17 = droplevels(subset(caro, sampling.session == "spring.2017"))
x = cbind(caro.s17$LD1, caro.s17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.s17$survival)
fit = Tps(x2, survival)
caro.s17.out = predictSurface(fit)
image.plot(caro.s17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Spring 2017")
contour(caro.s17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)
#0.25

caro.f17 = droplevels(subset(caro, sampling.session == "fall.2017"))
x = cbind(caro.f17$LD1, caro.f17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(caro.f17$survival)
fit = Tps(x2, survival)
caro.f17.out = predictSurface(fit)
image.plot(caro.f17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(caro$LD1)-1, max(caro$LD1)+1),
           ylim = c(min(caro$LD3)-1, max(caro$LD3)+1),
           main="Fall 2017")
contour(caro.f17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(caro.hull)
abline(line, col = 'red', lwd=2)
#0.35

###16a(ii.) GAMs through time----
par(mfrow=c(2,3))
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = caro, 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "Cumuative",
     col = "forestgreen",las = 1,xaxt="n", 
     xlim=c(min(caro$major.axis.transect), max(caro$major.axis.transect)),
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(caro$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#fall.2015
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = subset(caro, sampling.session=="fall.2015"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2015",
     col = "forestgreen",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2015")$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2015" & survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#spring.2016
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = subset(caro, sampling.session=="spring.2016"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "spring.2016",
     col = "forestgreen",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="spring.2016")$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, sampling.session=="spring.2016" & survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#fall.2016
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = subset(caro, sampling.session=="fall.2016"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2016",
     col = "forestgreen",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2016")$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2016" & survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#spring.2017
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = subset(caro, sampling.session=="spring.2017"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "spring.2017",
     col = "forestgreen",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="spring.2017")$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, sampling.session=="spring.2017" & survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#fall.2017
caro.gam <- gam(survival ~ s(major.axis.transect), 
                data = subset(caro, sampling.session=="fall.2017"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2017",
     col = "forestgreen",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2017")$major.axis.transect, side=1, col = "forestgreen", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2017" & survival == "1")$major.axis.transect, side=3, col = "forestgreen", lwd=0.01)

#dev.off()

###16a(iii) Models----
#### . GLMs ----

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = caro, 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#time as random effect (linear)
caro.glm.lin.ran <- glmer(survival ~ scale(major.axis.transect) + 
                            (1|sampling.session), 
                          data = caro, 
                          family = binomial(link = "logit"))
caro.glm.lin.ran.drop <- glmer(survival ~ 1 + 
                                 (1|sampling.session), 
                               data = caro, 
                               family = binomial(link = "logit"))
s.caro.glm.lin.ran = anova(caro.glm.lin.ran,caro.glm.lin.ran.drop)
s.caro.glm.lin.ran

#standard model (including nonlinear term)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+
                       I(scale(major.axis.transect)^2), 
                     data = caro, 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#time as random effect (including nonlinear term)
caro.glm.full.ran <- glmer(survival ~ scale(major.axis.transect)+
                             I(scale(major.axis.transect)^2)+
                             (1|sampling.session), 
                           data = caro, 
                           family = binomial(link = "logit"))
#drop term of interest
caro.glm.full.ran.drop <- glmer(survival ~ scale(major.axis.transect)+
                                  #I(scale(major.axis.transect)^2)+
                                  (1|sampling.session), 
                                data = caro, 
                                family = binomial(link = "logit"))
s.caro.glm.full.ran = anova(caro.glm.full.ran, caro.glm.full.ran.drop)
s.caro.glm.full.ran

#### . LMs ----

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = caro)
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#time as random effect (linear)
caro.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(major.axis.transect) + (1|sampling.session), 
                         data = caro)
s.caro.ols.lin.ran = summary(caro.ols.lin.ran)
s.caro.ols.lin.ran

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+
                      +I(scale(major.axis.transect)^2), 
                    data = caro)
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full
#time as random effect (including nonlinear and correlational terms)
caro.ols.full.ran <- lmer(survival ~ scale(major.axis.transect)+
                            I(scale(major.axis.transect)^2)+
                            (1|sampling.session), 
                          data = caro)
s.caro.ols.full.ran = summary(caro.ols.full.ran)
s.caro.ols.full.ran$coefficients

#### . Results table----
#let's combine all of this together
tab.caro
carolinensis.table = tab.caro
carolinensis.table$major.axis.transect.lin <- NA
carolinensis.table$major.axis.transect.lin.SE <- NA
carolinensis.table$major.axis.transect.lin.p <- NA
carolinensis.table$major.axis.transect.quad <- NA
carolinensis.table$major.axis.transect.quad.SE <- NA
carolinensis.table$major.axis.transect.quad.p <- NA
carolinensis.table

s.caro.glm.full
#add in the standard model results first
carolinensis.table[6,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[6,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[6,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[6,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[6,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[6,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

s.caro.ols.full.ran$coefficients
#and now the random effects models
carolinensis.table[7,2] <- round(s.caro.ols.lin.ran$coefficients[2,1],3)
carolinensis.table[7,3] <- round(s.caro.ols.lin.ran$coefficients[2,2],3)
carolinensis.table[7,4] <- round(s.caro.glm.lin.ran[2,8],3)
carolinensis.table[7,5] <- round(s.caro.ols.full.ran$coefficients[3,1],3)
carolinensis.table[7,6] <- round(s.caro.ols.full.ran$coefficients[3,2],3)
carolinensis.table[7,7] <- round(s.caro.glm.full.ran[2,8],3)
carolinensis.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = subset(caro, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                     data = subset(caro, sampling.session == "fall.2015"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[1,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[1,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[1,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[1,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[1,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[1,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#spring 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = subset(caro, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                     data = subset(caro, sampling.session == "spring.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[2,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[2,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[2,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[2,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[2,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[2,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#fall 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = subset(caro, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                     data = subset(caro, sampling.session == "fall.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[3,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[3,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[3,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[3,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[3,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[3,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#spring 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = subset(caro, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                     data = subset(caro, sampling.session == "spring.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[4,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[4,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[4,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[4,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[4,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[4,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#fall 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                    data = subset(caro, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                     data = subset(caro, sampling.session == "fall.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                   data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[5,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[5,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[5,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[5,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[5,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[5,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#### . Final table----
carolinensis.table
#we need to double the quadratic coefficients (and associated SE's)
carolinensis.table$major.axis.transect.quad = carolinensis.table$major.axis.transect.quad*2
carolinensis.table$major.axis.transect.quad.SE = carolinensis.table$major.axis.transect.quad.SE*2
carolinensis.table
carolinensis.major.axis.transect.selection = carolinensis.table

##16b. sagrei ----

#plot if needed
#pdf("sagrei_majoraxis.pdf")

###16b(i) Transect plot----
par(mfrow = c(1,3))
#2D surface plot
image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="Identify fitness peak and trough")
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

sag$rel.fit = sag.fit$fitted.values
sag.fit.max = sag[which.max(sag$rel.fit),]
points(sag.fit.max$LD1, sag.fit.max$LD3, pch=21, bg="red", cex=2)
sag.fit.min = sag[which.min(sag$rel.fit),]
points(sag.fit.min$LD1, sag.fit.min$LD3, pch=21, bg="blue", cex=2)

#do things
image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="Draw optimal fitness transect line")
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

max.centroid = rbind(sag.fit.max$LD3, sag.fit.max$LD1)
min.centroid = rbind(sag.fit.min$LD3, sag.fit.min$LD1)
f = max.centroid[2,1]
g = max.centroid[1,1]
h = min.centroid[2,1]
j = min.centroid[1,1]

x <- sag$LD1
y <- sag$LD3

ys = c(f, g)
xs = c(h, j)
d = rbind(ys, xs)
d = as.data.frame(d)
d
line <- lm(V2~V1, d)

#draw an abline connecting the two species centroids
abline(line, col = 'red', lwd=2, lty=2)
#plot the two centroids
points(sag.fit.max$LD1, sag.fit.max$LD3, pch=21, bg="red", cex=2)
points(sag.fit.min$LD1, sag.fit.min$LD3, pch=21, bg="blue", cex=2)
#check data match
#points(d[,1], d[,2], bg = "black", pch = 21, cex=2)
#looks good

# function draws line between two random points
drawline <- function(V2, V1){
  line <- lm(V2~V1, d)
  abline(line, col = "red", lwd=2, lty=2)
}

# function to measure distance from 'measure_point' to line
dist_to_line <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 

# function to measure distance to 'ref_point'
dist_to_ref <- function(r,t) {
  u <- (r[1]-t[1])^2
  v <- (r[2]-t[2])^2
  j <- sqrt(u+v)
}

line
m = summary(line)
m$coefficients
aa = m$coefficients[1,1]
bb = m$coefficients[2,1]
cc = min(sag$LD1)-2
f = aa+bb*cc
ref_point = c(cc,f)
line_point = c(d[2,1],d[2,2])
line_point


image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="Set reference point")
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

#draw reference point for measuring
points(line_point[1],line_point[2], bg = "blue", pch = 21, cex=2)
points(ref_point[1],ref_point[2], bg = "lightblue", pch = 21, cex=2)
points(sag.fit.max$LD1, sag.fit.max$LD3, pch=21, bg="red", cex=2)
sag_point = c(d[1,1],d[1,2])
drawline(c(ref_point[1],line_point[1]),c(ref_point[2],line_point[2]))
#ok getting there

measure_point <- 0
i = 1
a <- NA

length = length(sag$svl)

image.plot(sag.out, xlab = "LD Axis 1", 
           ylab = "LD Axis 3", 
           xlim = c(min(res$LD1)-1, max(res$LD1)+1),
           ylim = c(min(res$LD3)-1, max(res$LD3)+1),
           main="Plot distribution of all individuals")
contour(sag.out, add=T, col = "gray50", drawlabels=FALSE, nlevels=10)
lines(comm.hull, lwd = 2, col = "black")

#visualize position of all individuals relative to the major.axis.transect
for(i in 1:length){
  measure_point <- c(x[i],y[i])
  points(measure_point[1],measure_point[2], bg = "gray75", pch = 21)
  d <- dist_to_line(measure_point, ref_point, line_point)
  j <- dist_to_ref(measure_point,ref_point)
  dist <- sqrt((j^2 - d^2))
  a[i]<- dist
}
abline(line, col = 'red', lwd=2, lty=2)

#data_plotting <- data.frame(a,1)
#plot major.axis.transect data
#plot(data_plotting[,1],data_plotting[,2], 
#     xlim = c(min(data_plotting[,1]-2),
#              max(data_plotting[,1])),
#     yaxt='n', xaxt = 'n',
#     ylab="",
#     xlab = "Position along optimal fitness transect")
#abline(h=1, col="red", lty=2)
#head(a)

#combine major.axis.transect data to main dataframe
sag$major.axis.transect <- a

sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = sag, 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Optimal fitness transect", main = "Selection surface across
optimal fitness transect",
     col = "brown",las = 1,xaxt="n", ylim = c(0,0.4)
)
rug(sag$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)
#dev.off()

#models
sag.glm <- glmer(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2) + (1|sampling.session), 
                 data = sag, 
                 family = binomial(link = "logit"))
check_model(sag.glm)
sag.glm = Anova(sag.glm)
sag.glm

sag.ols <- lmer(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2) + (1|sampling.session), 
                data = sag)
summary(sag.ols)

###16b(i) Surfaces ----

par(mfrow=c(2,3))
image.plot(sag.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Cumulative")
contour(sag.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)


#surfaces through time
sag.f15 = droplevels(subset(sag, sampling.session == "fall.2015"))
x = cbind(sag.f15$LD1, sag.f15$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f15$survival)
sag.15.fit = Tps(x2, survival)
sag.f15.out = predictSurface(sag.15.fit)
image.plot(sag.f15.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,
           xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2015")
contour(sag.f15.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)
#0.35

sag.s16 = droplevels(subset(sag, sampling.session == "spring.2016"))
x = cbind(sag.s16$LD1, sag.s16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.s16$survival)
fit = Tps(x2, survival)
sag.s16.out = predictSurface(fit)
image.plot(sag.s16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Spring 2016")
contour(sag.s16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)
#0.35

sag.f16 = droplevels(subset(sag, sampling.session == "fall.2016"))
x = cbind(sag.f16$LD1, sag.f16$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f16$survival)
fit = Tps(x2, survival)
sag.f16.out = predictSurface(fit)
image.plot(sag.f16.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2016")
contour(sag.f16.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)
#0.6

sag.s17 = droplevels(subset(sag, sampling.session == "spring.2017"))
x = cbind(sag.s17$LD1, sag.s17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.s17$survival)
fit = Tps(x2, survival)
sag.s17.out = predictSurface(fit)
image.plot(sag.s17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Spring 2017")
contour(sag.s17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)
#0.25

sag.f17 = droplevels(subset(sag, sampling.session == "fall.2017"))
x = cbind(sag.f17$LD1, sag.f17$LD3)
x2 = as.matrix(x, ncol=2)
survival = c(sag.f17$survival)
fit = Tps(x2, survival)
sag.f17.out = predictSurface(fit)
image.plot(sag.f17.out, 
           zlim2 = c(0,0.65),
           add.legend=F,
           xlab = "", ylab = "", zlab = "",
           box=F
           ,xlim = c(min(sag$LD1)-1, max(sag$LD1)+1),
           ylim = c(min(sag$LD3)-1, max(sag$LD3)+1),
           main="Fall 2017")
contour(sag.f17.out, add=T, col = "gray25", nlevels=5,
        drawlabels = TRUE)
lines(sag.hull)
abline(line, col = 'red', lwd=2)
#0.35


###16b(ii.) GAMs through time----
par(mfrow=c(2,3))
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = sag, 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "Cumuative",
     col = "brown",las = 1,xaxt="n", 
     xlim=c(min(sag$major.axis.transect), max(sag$major.axis.transect)),
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(sag$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#fall.2015
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = subset(sag, sampling.session=="fall.2015"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2015",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2015")$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2015" & survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#spring.2016
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = subset(sag, sampling.session=="spring.2016"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "spring.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="spring.2016")$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="spring.2016" & survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#fall.2016
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = subset(sag, sampling.session=="fall.2016"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2016")$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2016" & survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#spring.2017
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = subset(sag, sampling.session=="spring.2017"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "spring.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="spring.2017")$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="spring.2017" & survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#fall.2017
sag.gam <- gam(survival ~ s(major.axis.transect), 
               data = subset(sag, sampling.session=="fall.2017"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis major.axis.transect", 
     main = "fall.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2017")$major.axis.transect, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2017" & survival == "1")$major.axis.transect, side=3, col = "brown", lwd=0.01)

#dev.off()
###16b(iii) Models----
#### . GLMs ----

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = sag, 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#time as random effect (linear)
sag.glm.lin.ran <- glmer(survival ~ scale(major.axis.transect) + 
                           (1|sampling.session), 
                         data = sag, 
                         family = binomial(link = "logit"))
sag.glm.lin.ran.drop <- glmer(survival ~ 1 + 
                                (1|sampling.session), 
                              data = sag, 
                              family = binomial(link = "logit"))
s.sag.glm.lin.ran = anova(sag.glm.lin.ran,sag.glm.lin.ran.drop)
s.sag.glm.lin.ran

#standard model (including nonlinear term)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = sag, 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#time as random effect (including nonlinear term)
sag.glm.full.ran <- glmer(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2)+
                            (1|sampling.session), 
                          data = sag, 
                          family = binomial(link = "logit"))
#drop term of interest
sag.glm.full.ran.drop <- glmer(survival ~ scale(major.axis.transect)+
                                 (1|sampling.session), 
                               data = sag, 
                               family = binomial(link = "logit"))
s.sag.glm.full.ran = anova(sag.glm.full.ran, sag.glm.full.ran.drop)
s.sag.glm.full.ran

#### . LMs ----

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = sag)
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#time as random effect (linear)
sag.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(major.axis.transect) + (1|sampling.session), 
                        data = sag)
s.sag.ols.lin.ran = summary(sag.ols.lin.ran)
s.sag.ols.lin.ran

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+
                     +I(scale(major.axis.transect)^2), 
                   data = sag)
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full
#time as random effect (including nonlinear and correlational terms)
sag.ols.full.ran <- lmer(survival ~ scale(major.axis.transect)+
                           I(scale(major.axis.transect)^2)+
                           (1|sampling.session), 
                         data = sag)
s.sag.ols.full.ran = summary(sag.ols.full.ran)
.sag.ols.full.ran$coefficients

#### . Results table----
#let's combine all of this together
tab.sag
sagrei.table = tab.sag
sagrei.table$major.axis.transect.lin <- NA
sagrei.table$major.axis.transect.lin.SE <- NA
sagrei.table$major.axis.transect.lin.p <- NA
sagrei.table$major.axis.transect.quad <- NA
sagrei.table$major.axis.transect.quad.SE <- NA
sagrei.table$major.axis.transect.quad.p <- NA
sagrei.table

s.sag.glm.full
#add in the standard model results first
sagrei.table[6,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[6,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[6,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[6,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[6,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[6,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

s.sag.ols.full.ran$coefficients
#and now the random effects models
sagrei.table[7,2] <- round(s.sag.ols.lin.ran$coefficients[2,1],3)
sagrei.table[7,3] <- round(s.sag.ols.lin.ran$coefficients[2,2],3)
sagrei.table[7,4] <- round(s.sag.glm.lin.ran[2,8],3)
sagrei.table[7,5] <- round(s.sag.ols.full.ran$coefficients[3,1],3)
sagrei.table[7,6] <- round(s.sag.ols.full.ran$coefficients[3,2],3)
sagrei.table[7,7] <- round(s.sag.glm.full.ran[2,8],3)
sagrei.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = subset(sag, sampling.session == "fall.2015"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(sag, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                   data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[1,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[1,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[1,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[1,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[1,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[1,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#spring 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = subset(sag, sampling.session == "spring.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(sag, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                   data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[2,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[2,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[2,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[2,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[2,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[2,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#fall 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = subset(sag, sampling.session == "fall.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(sag, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                   data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[3,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[3,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[3,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[3,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[3,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[3,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#spring 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = subset(sag, sampling.session == "spring.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(sag, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                   data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[4,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[4,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[4,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[4,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[4,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[4,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#fall 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(major.axis.transect), 
                   data = subset(sag, sampling.session == "fall.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                    data = subset(sag, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(major.axis.transect), 
                  data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(major.axis.transect)+I(scale(major.axis.transect)^2), 
                   data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[5,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[5,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[5,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[5,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[5,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[5,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#### . Final table----
sagrei.table
#we need to double the quadratic coefficients (and associated SE's)
sagrei.table$major.axis.transect.quad = sagrei.table$major.axis.transect.quad*2
sagrei.table$major.axis.transect.quad.SE = sagrei.table$major.axis.transect.quad.SE*2
sagrei.table
sagrei.major.axis.transect.selection = sagrei.table



#17. Projection pursuit regression ----

## 17a. Anolis sagrei----

sag.mod.ppr <- ppr(survival ~ LD1 + LD3, nterms = 1, max.terms = 5, data = sag, sm.method = 'gcvspline', gcvpen = 1)
summary(sag.mod.ppr)

# Call:
# ppr(formula = survival ~ LD1 + LD3, data = sag, nterms = 1, max.terms = 5, 
# sm.method = "gcvspline", gcvpen = 1)

# Goodness of fit:
# 1 terms  2 terms  3 terms  4 terms  5 terms 
# 255.9598 255.9237   0.0000   0.0000   0.0000 

# Projection direction vectors ('alpha'):
# LD1       LD3 
# 0.6082112 0.7937752 

# Coefficients of ridge terms ('beta'):
# term 1 
# 0.03573184 

# Equivalent df for ridge terms:
# term 1 
# 3.14 

par(mfrow=c(2,2))
plot(sag.mod.ppr)
plot(update(sag.mod.ppr, bass = 5), main = "update(..., bass = 5)")
plot(update(sag.mod.ppr, sm.method = "gcv", gcvpen = 2),
     main = "update(..., sm.method=\"gcv\", gcvpen=2)")

sag.projected_z <- (sag.mod.ppr$alpha[1]*sag$LD1) + (sag.mod.ppr$alpha[2]*sag$LD3)
sag$projected_z <- (sag.projected_z - mean(sag.projected_z))/sd(sag.projected_z)
hist(sag$projected_z)

image.plot(sag.out)
sag.mod.ppr$yb; sag.mod.ppr$ys
summary(sag.mod.ppr)

sag.glm01 <- glm(survival ~ projected_z + I(projected_z^2), 
                 data = sag, family = binomial(link = 'logit'))
Anova(sag.glm01, type="III")

#gam
sag.gam <- gam(survival ~ s(projected_z), 
               data = sag, 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis transect", 
     main = "Cumulative: PPR",
     col = "brown",las = 1,
     #xaxt="n", 
     xlim=c(min(sag$projected_z), max(sag$projected_z))
     #,ylim = c(0,max(sag.gam$fitted.values)+0.1)
)

# Analysis of Deviance Table (Type II tests)
##17b. GAMs through time----
par(mfrow=c(2,3))
sag.gam <- gam(survival ~ s(projected_z), 
               data = sag, 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "Cumuative",
     col = "brown",las = 1,xaxt="n", 
     xlim=c(min(sag$projected_z), max(sag$projected_z)),
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(sag$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2015
sag.gam <- gam(survival ~ s(projected_z), 
               data = subset(sag, sampling.session=="fall.2015"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2015",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2015")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2015" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#spring.2016
sag.gam <- gam(survival ~ s(projected_z), 
               data = subset(sag, sampling.session=="spring.2016"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "spring.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="spring.2016")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="spring.2016" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2016
sag.gam <- gam(survival ~ s(projected_z), 
               data = subset(sag, sampling.session=="fall.2016"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2016")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2016" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#spring.2017
sag.gam <- gam(survival ~ s(projected_z), 
               data = subset(sag, sampling.session=="spring.2017"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "spring.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="spring.2017")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="spring.2017" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2017
sag.gam <- gam(survival ~ s(projected_z), 
               data = subset(sag, sampling.session=="fall.2017"), 
               family = binomial(link = "logit"), 
               method = "GCV.Cp")

plot(sag.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(sag.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(sag.gam$fitted.values)+0.1)
)
rug(subset(sag, sampling.session=="fall.2017")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(sag, sampling.session=="fall.2017" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

##17c. Models----
###17c(ii) GLMs ----

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = sag, 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#time as random effect (linear)
sag.glm.lin.ran <- glmer(survival ~ scale(projected_z) + 
                           (1|sampling.session), 
                         data = sag, 
                         family = binomial(link = "logit"))
sag.glm.lin.ran.drop <- glmer(survival ~ 1 + 
                                (1|sampling.session), 
                              data = sag, 
                              family = binomial(link = "logit"))
s.sag.glm.lin.ran = anova(sag.glm.lin.ran,sag.glm.lin.ran.drop)
s.sag.glm.lin.ran

#standard model (including nonlinear term)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = sag, 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#time as random effect (including nonlinear term)
sag.glm.full.ran <- glmer(survival ~ scale(projected_z)+I(scale(projected_z)^2)+
                            (1|sampling.session), 
                          data = sag, 
                          family = binomial(link = "logit"))
#drop term of interest
sag.glm.full.ran.drop <- glmer(survival ~ scale(projected_z)+
                                 (1|sampling.session), 
                               data = sag, 
                               family = binomial(link = "logit"))
s.sag.glm.full.ran = anova(sag.glm.full.ran, sag.glm.full.ran.drop)
s.sag.glm.full.ran

###17c(iii) LMs ----

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = sag)
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#time as random effect (linear)
sag.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(projected_z) + (1|sampling.session), 
                        data = sag)
s.sag.ols.lin.ran = summary(sag.ols.lin.ran)
s.sag.ols.lin.ran

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+
                     +I(scale(projected_z)^2), 
                   data = sag)
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full
#time as random effect (including nonlinear and correlational terms)
sag.ols.full.ran <- lmer(survival ~ scale(projected_z)+
                           I(scale(projected_z)^2)+
                           (1|sampling.session), 
                         data = sag)
s.sag.ols.full.ran = summary(sag.ols.full.ran)
.sag.ols.full.ran$coefficients

###17c(iv) Results table----
#let's combine all of this together
tab.sag
sagrei.table = tab.sag
sagrei.table$ppr.lin <- NA
sagrei.table$ppr.lin.SE <- NA
sagrei.table$ppr.lin.p <- NA
sagrei.table$ppr.quad <- NA
sagrei.table$ppr.quad.SE <- NA
sagrei.table$ppr.quad.p <- NA
sagrei.table

s.sag.glm.full
#add in the standard model results first
sagrei.table[6,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[6,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[6,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[6,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[6,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[6,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

s.sag.ols.full.ran$coefficients
#and now the random effects models
sagrei.table[7,2] <- round(s.sag.ols.lin.ran$coefficients[2,1],3)
sagrei.table[7,3] <- round(s.sag.ols.lin.ran$coefficients[2,2],3)
sagrei.table[7,4] <- round(s.sag.glm.lin.ran[2,8],3)
sagrei.table[7,5] <- round(s.sag.ols.full.ran$coefficients[3,1],3)
sagrei.table[7,6] <- round(s.sag.ols.full.ran$coefficients[3,2],3)
sagrei.table[7,7] <- round(s.sag.glm.full.ran[2,8],3)
sagrei.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = subset(sag, sampling.session == "fall.2015"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(sag, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                   data = subset(sag, sampling.session == "fall.2015"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[1,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[1,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[1,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[1,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[1,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[1,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#spring 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = subset(sag, sampling.session == "spring.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(sag, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                   data = subset(sag, sampling.session == "spring.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[2,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[2,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[2,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[2,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[2,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[2,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#fall 2016

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = subset(sag, sampling.session == "fall.2016"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(sag, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                   data = subset(sag, sampling.session == "fall.2016"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[3,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[3,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[3,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[3,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[3,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[3,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#spring 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = subset(sag, sampling.session == "spring.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(sag, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                   data = subset(sag, sampling.session == "spring.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[4,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[4,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[4,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[4,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[4,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[4,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

#fall 2017

#standard model (linear)
sag.glm.lin <- glm(survival ~ scale(projected_z), 
                   data = subset(sag, sampling.session == "fall.2017"), 
                   family = binomial(link = "logit"))
s.sag.glm.lin = Anova(sag.glm.lin, type="III")
s.sag.glm.lin

#standard model (including nonlinear and correlational terms)
sag.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(sag, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.sag.glm.full = Anova(sag.glm.full, type="III")
s.sag.glm.full

#standard model (linear)
sag.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                  data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.lin = summary(sag.ols.lin)
s.sag.ols.lin

#standard model (including nonlinear and correlational terms)
sag.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                   data = subset(sag, sampling.session == "fall.2017"))
s.sag.ols.full = summary(sag.ols.full)
s.sag.ols.full

#add in the standard model results first
sagrei.table[5,2] <- round(s.sag.ols.lin$coefficients[2,1],3)
sagrei.table[5,3] <- round(s.sag.ols.lin$coefficients[2,2],3)
sagrei.table[5,4] <- round(s.sag.glm.lin[1,3],3)
sagrei.table[5,5] <- round(s.sag.ols.full$coefficients[3,1],3)
sagrei.table[5,6] <- round(s.sag.ols.full$coefficients[3,2],3)
sagrei.table[5,7] <- round(s.sag.glm.full[2,3],3)
sagrei.table

###17c(v) PPR final table----
sagrei.table
#we need to double the quadratic coefficients (and associated SE's)
sagrei.table$ppr.quad = sagrei.table$ppr.quad*2
sagrei.table$ppr.quad.SE = sagrei.table$ppr.quad.SE*2
sagrei.table
sagrei.ppr.selection = sagrei.table


## 17d. Anolis carolinensis----

caro.mod.ppr <- ppr(survival ~ LD1 + LD3, nterms = 1, max.terms = 5, data = caro, sm.method = 'gcvspline', gcvpen = 1)
summary(caro.mod.ppr)

# Call:
# ppr(formula = survival ~ LD1 + LD3, data = caro, nterms = 1, max.terms = 5, 
# sm.method = "gcvspline", gcvpen = 1)

# Goodness of fit:
# 1 terms  2 terms  3 terms  4 terms  5 terms 
# 255.9598 255.9237   0.0000   0.0000   0.0000 

# Projection direction vectors ('alpha'):
# LD1       LD3 
# 0.6082112 0.7937752 

# Coefficients of ridge terms ('beta'):
# term 1 
# 0.03573184 

# Equivalent df for ridge terms:
# term 1 
# 3.14 

par(mfrow=c(2,2))
plot(caro.mod.ppr)
plot(update(caro.mod.ppr, bass = 5), main = "update(..., bass = 5)")
plot(update(caro.mod.ppr, sm.method = "gcv", gcvpen = 2),
     main = "update(..., sm.method=\"gcv\", gcvpen=2)")

caro.projected_z <- (caro.mod.ppr$alpha[1]*caro$LD1) + (caro.mod.ppr$alpha[2]*caro$LD3)
caro$projected_z <- (caro.projected_z - mean(caro.projected_z))/sd(caro.projected_z)
#hist(caro$projected_z)

image.plot(caro.out)
caro.mod.ppr$yb; caro.mod.ppr$ys
abline(caro.mod.ppr$ys, caro.mod.ppr$yb)
summary(caro.mod.ppr)

caro.glm01 <- glm(survival ~ projected_z + I(projected_z^2), 
                  data = caro, family = binomial(link = 'logit'))
Anova(caro.glm01, type="III")

#gam
caro.gam <- gam(survival ~ s(projected_z), 
                data = caro, 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis transect", 
     main = "Cumulative: PPR",
     col = "brown",las = 1,
     #xaxt="n", 
     xlim=c(min(caro$projected_z), max(caro$projected_z))
     #,ylim = c(0,max(caro.gam$fitted.values)+0.1)
)

# Analysis of Deviance Table (Type II tests)
##17e. GAMs through time----
par(mfrow=c(2,3))
caro.gam <- gam(survival ~ s(projected_z), 
                data = caro, 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "Cumuative",
     col = "brown",las = 1,xaxt="n", 
     xlim=c(min(caro$projected_z), max(caro$projected_z)),
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(caro$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2015
caro.gam <- gam(survival ~ s(projected_z), 
                data = subset(caro, sampling.session=="fall.2015"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2015",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2015")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2015" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#spring.2016
caro.gam <- gam(survival ~ s(projected_z), 
                data = subset(caro, sampling.session=="spring.2016"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "spring.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="spring.2016")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session=="spring.2016" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2016
caro.gam <- gam(survival ~ s(projected_z), 
                data = subset(caro, sampling.session=="fall.2016"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2016",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2016")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2016" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#spring.2017
caro.gam <- gam(survival ~ s(projected_z), 
                data = subset(caro, sampling.session=="spring.2017"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "spring.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="spring.2017")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session=="spring.2017" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

#fall.2017
caro.gam <- gam(survival ~ s(projected_z), 
                data = subset(caro, sampling.session=="fall.2017"), 
                family = binomial(link = "logit"), 
                method = "GCV.Cp")

plot(caro.gam, n = 500, se = 1, seWithMean = TRUE, rug = FALSE, 
     shift = mean(predict(caro.gam)),
     trans = function(x){exp(x)/(1+exp(x))}, ylab = "Survival prob.", 
     xlab = "Major axis projected_z", 
     main = "fall.2017",
     col = "brown",las = 1,xaxt="n", 
     ylim = c(0,max(caro.gam$fitted.values)+0.1)
)
rug(subset(caro, sampling.session=="fall.2017")$projected_z, side=1, col = "brown", lwd=0.01)
rug(subset(caro, sampling.session=="fall.2017" & survival == "1")$projected_z, side=3, col = "brown", lwd=0.01)

##17f. Models----
###17f(ii) GLMs ----

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = caro, 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#time as random effect (linear)
caro.glm.lin.ran <- glmer(survival ~ scale(projected_z) + 
                            (1|sampling.session), 
                          data = caro, 
                          family = binomial(link = "logit"))
caro.glm.lin.ran.drop <- glmer(survival ~ 1 + 
                                 (1|sampling.session), 
                               data = caro, 
                               family = binomial(link = "logit"))
s.caro.glm.lin.ran = anova(caro.glm.lin.ran,caro.glm.lin.ran.drop)
s.caro.glm.lin.ran

#standard model (including nonlinear term)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = caro, 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#time as random effect (including nonlinear term)
caro.glm.full.ran <- glmer(survival ~ scale(projected_z)+I(scale(projected_z)^2)+
                             (1|sampling.session), 
                           data = caro, 
                           family = binomial(link = "logit"))
#drop term of interest
caro.glm.full.ran.drop <- glmer(survival ~ scale(projected_z)+
                                  (1|sampling.session), 
                                data = caro, 
                                family = binomial(link = "logit"))
s.caro.glm.full.ran = anova(caro.glm.full.ran, caro.glm.full.ran.drop)
s.caro.glm.full.ran

###17f(iii) LMs ----

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = caro)
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#time as random effect (linear)
caro.ols.lin.ran <- lmer(survival/mean(survival) ~ scale(projected_z) + (1|sampling.session), 
                         data = caro)
s.caro.ols.lin.ran = summary(caro.ols.lin.ran)
s.caro.ols.lin.ran

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+
                      +I(scale(projected_z)^2), 
                    data = caro)
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full
#time as random effect (including nonlinear and correlational terms)
caro.ols.full.ran <- lmer(survival ~ scale(projected_z)+
                            I(scale(projected_z)^2)+
                            (1|sampling.session), 
                          data = caro)
s.caro.ols.full.ran = summary(caro.ols.full.ran)
.caro.ols.full.ran$coefficients

###17f(iv) Results table----
#let's combine all of this together
tab.caro
carolinensis.table = tab.caro
carolinensis.table$ppr.lin <- NA
carolinensis.table$ppr.lin.SE <- NA
carolinensis.table$ppr.lin.p <- NA
carolinensis.table$ppr.quad <- NA
carolinensis.table$ppr.quad.SE <- NA
carolinensis.table$ppr.quad.p <- NA
carolinensis.table

s.caro.glm.full
#add in the standard model results first
carolinensis.table[6,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[6,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[6,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[6,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[6,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[6,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

s.caro.ols.full.ran$coefficients
#and now the random effects models
carolinensis.table[7,2] <- round(s.caro.ols.lin.ran$coefficients[2,1],3)
carolinensis.table[7,3] <- round(s.caro.ols.lin.ran$coefficients[2,2],3)
carolinensis.table[7,4] <- round(s.caro.glm.lin.ran[2,8],3)
carolinensis.table[7,5] <- round(s.caro.ols.full.ran$coefficients[3,1],3)
carolinensis.table[7,6] <- round(s.caro.ols.full.ran$coefficients[3,2],3)
carolinensis.table[7,7] <- round(s.caro.glm.full.ran[2,8],3)
carolinensis.table

#and now let's start the individual sampling sessions

#fall 2015

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = subset(caro, sampling.session == "fall.2015"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = subset(caro, sampling.session == "fall.2015"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(caro, sampling.session == "fall.2015"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[1,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[1,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[1,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[1,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[1,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[1,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#spring 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = subset(caro, sampling.session == "spring.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = subset(caro, sampling.session == "spring.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(caro, sampling.session == "spring.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[2,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[2,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[2,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[2,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[2,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[2,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#fall 2016

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = subset(caro, sampling.session == "fall.2016"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = subset(caro, sampling.session == "fall.2016"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(caro, sampling.session == "fall.2016"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[3,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[3,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[3,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[3,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[3,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[3,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#spring 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = subset(caro, sampling.session == "spring.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = subset(caro, sampling.session == "spring.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(caro, sampling.session == "spring.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[4,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[4,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[4,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[4,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[4,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[4,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

#fall 2017

#standard model (linear)
caro.glm.lin <- glm(survival ~ scale(projected_z), 
                    data = subset(caro, sampling.session == "fall.2017"), 
                    family = binomial(link = "logit"))
s.caro.glm.lin = Anova(caro.glm.lin, type="III")
s.caro.glm.lin

#standard model (including nonlinear and correlational terms)
caro.glm.full <- glm(survival ~ scale(projected_z)+I(scale(projected_z)^2), 
                     data = subset(caro, sampling.session == "fall.2017"), 
                     family = binomial(link = "logit"))
s.caro.glm.full = Anova(caro.glm.full, type="III")
s.caro.glm.full

#standard model (linear)
caro.ols.lin <- lm(survival/mean(survival) ~ scale(projected_z), 
                   data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.lin = summary(caro.ols.lin)
s.caro.ols.lin

#standard model (including nonlinear and correlational terms)
caro.ols.full <- lm(survival/mean(survival) ~ scale(projected_z)+I(scale(projected_z)^2), 
                    data = subset(caro, sampling.session == "fall.2017"))
s.caro.ols.full = summary(caro.ols.full)
s.caro.ols.full

#add in the standard model results first
carolinensis.table[5,2] <- round(s.caro.ols.lin$coefficients[2,1],3)
carolinensis.table[5,3] <- round(s.caro.ols.lin$coefficients[2,2],3)
carolinensis.table[5,4] <- round(s.caro.glm.lin[1,3],3)
carolinensis.table[5,5] <- round(s.caro.ols.full$coefficients[3,1],3)
carolinensis.table[5,6] <- round(s.caro.ols.full$coefficients[3,2],3)
carolinensis.table[5,7] <- round(s.caro.glm.full[2,3],3)
carolinensis.table

###17f(v) PPR final table----
carolinensis.table
#we need to double the quadratic coefficients (and associated SE's)
carolinensis.table$ppr.quad = carolinensis.table$ppr.quad*2
carolinensis.table$ppr.quad.SE = carolinensis.table$ppr.quad.SE*2
carolinensis.table
carolinensis.ppr.selection = carolinensis.table


#18. Final tables----

## 18a. Anolis sagrei ----

### Major transect----
sagrei.major.axis.transect.selection
write.table(sagrei.major.axis.transect.selection, "clipboard", sep="\t", row.names=TRUE)
###PPR----
sagrei.ppr.selection
write.table(sagrei.ppr.selection, "clipboard", sep="\t", row.names=TRUE)
###Euclidean----
sagrei.euclidean.selection
write.table(sagrei.euclidean.selection, "clipboard", sep="\t", row.names=TRUE)

## 18b. Anolis carolinensis----

### Major transect----
carolinensis.major.axis.transect.selection
write.table(carolinensis.major.axis.transect.selection, "clipboard", sep="\t", row.names=TRUE)
###PPR----
carolinensis.ppr.selection
write.table(carolinensis.ppr.selection, "clipboard", sep="\t", row.names=TRUE)
###Euclidean----
carolinensis.euclidean.selection
write.table(carolinensis.euclidean.selection, "clipboard", sep="\t", row.names=TRUE)
