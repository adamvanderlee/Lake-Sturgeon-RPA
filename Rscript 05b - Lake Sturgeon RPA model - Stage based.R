################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Lake Sturgeon
# Recovery Potential Assessment - Estimates of Allowable harm, recovery targets
# and critical habitat
#-------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())  ### clear the workspace

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load Libraries
library(pacman)
p_load(popbio)
p_load(ggplot2)
p_load(ggmap)
p_load(RgoogleMaps)
p_load(car)
p_load(nlme)
p_load(parallel)

#-------------------------------------------------------------------------------
# Functions

# Files to source 
source("Rscript 01 - Parameter functions.r")  # sub routines for parameter fitting
source("Rscript 02 - Plot functions.r")       # sub routines for plotting
source("Rscript 03 - Matrix functions.r")     # sub routines for matrices
source("Rscript 04b - Simulation functions - stage based.r") # sub routines for model simulation

#-------------------------------------------------------------------------------
# ggplot theme

theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 11, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 10, family = "sans", colour = "black"),
        axis.text.y = element_text(size = 10, family = "sans", hjust = 0.6,
                                   angle = 90, colour = "black"),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        strip.text = element_text(size = 11, family = "sans", face = "bold"))

#-------------------------------------------------------------------------------
# Data

# life_history data

# Length data
VBGF <- read.csv("data/VBGF.csv", header = TRUE)
Fortin_VB <- read.csv("data/Fortin 1996 - Length-Lat.csv", header = TRUE)
Power_VB <- read.csv("data/Power 1997 - Length-Lat.csv", header = TRUE)

# L-W data
L_W <- read.csv("data/LW-regression.csv", header = TRUE)
L_W$a <- L_W$a / 1000 # Convert units to kg

# Maturity
MO_data <- read.csv("data/Bruch 2008 - Maturity.csv", header = TRUE)

# Mortality
M_data <- read.csv("data/Mortality.csv", header = TRUE)

# Fecundity
Fec_data <- read.csv("data/COSEWIC - Fecundity.csv", header = TRUE)

# Population data
pop_data <- read.csv("data/COSEWIC 2017 - Pop Status.csv", header = TRUE)

# Data from Tom Pratt at Goulais Bay - not included on Git hub
# GB_data <- read.csv("data/Goulais Bay.csv", header = TRUE)

#-------------------------------------------------------------------------------
# COSEWIC POPULATION DATA

with(pop_data, aggregate(Estimate, by = list(DU), length))

with(pop_data, aggregate(Estimate, by = list(DU), mean, na.rm = T))
with(pop_data, aggregate(Estimate, by = list(DU), geom_mean, na.rm = T))
with(pop_data, aggregate(Estimate, by = list(DU), sum, na.rm = T))

lapply(1:4, function(x){
  data <- pop_data[pop_data$DU == x, ]
  apply(data[, 14:20], 2, sum, na.rm = T)
})

lapply(1:4, function(x){
  data <- pop_data[pop_data$DU == x, ]
  apply(data[, 14:20], 2, sum, na.rm = T)/dim(data)[1]
})

lapply(1:4, function(x){
  data <- pop_data[pop_data$DU == x, ]
  table(data$Trajectory)
})

lapply(1:4, function(x){
  data <- pop_data[pop_data$DU == x, ]
  table(data$Recuitment)
})

#===============================================================================
# Map: with position of simulated population per DU

# Postional data
# Lat and Lon of each simulated population
pos_data <- as.data.frame(list(
  Lat = c(57.3, 49.9, 49.1, 50.5, 42.5, 45.4),
  Lon = c(102, 111.7, 95, 79.1, 82.7, 73.8),
  name = c("DU1: Western Hudson Bay", "DU2: South Saskatchewan River", "DU2: Lake of the Woods",
           "DU3: Southern Hudson Bay - James Bay", "DU4: Lake St. Clair", "DU4: St. Lawrence River")
))

map <- get_googlemap(center = c(-94.5, 52), zoom = 4, maptype = "terrain",
                     style = 'element:labels|visibility:off')

p <- ggmap(map)
p <- p + geom_point(data = pos_data, aes(x = -Lon, y = Lat), colour = "black", size = 1)
p <- p + geom_text(data = pos_data, aes(x = -Lon, y = Lat, label = name), colour = "black", nudge_y = -.5, size = 2.5)
#p <- p + geom_point(data = VBGF, aes(x = -Longitude, y = Latitude, colour = as.factor(DU)), size = 2)
p <- p + scale_x_continuous(limits = c(-120, -69))
p <- p + scale_y_continuous(limits = c(42, 60))
p <- p + labs(x = "Longitude", y = "Latitude")

#png("Figures/DU - Location.png", width = 5, height = 3, unit = "in", res = 300)
p
#dev.off()

# Clean up
rm(map, p)

#===============================================================================
# Parameters

DU <- c("DU1", "DU2a", "DU2b", "DU3", "DU4a", "DU4b")

# Initalize list to store parameter info
lh_data <- lapply(DU, function(x) list())   # life history data; parameters to est. growth, fecundity, and survival
names(lh_data) <- DU

Status <- c("Endangered", "Endangered", "Endangered", "Special Concern", "Threatened", "Threatened")
lh_data <- lapply(1:6, function(x) {lh_data[[x]]$Status = Status[x]; lh_data[[x]]})
lh_data <- lapply(lh_data, function(x) {x$Tmax = 100; x})

#-------------------------------------------------------------------------------
# GROWTH
#-------------------------------------------------------------------------------

# Hatch length 
H0.range <- c(6.5, 14) # COSEWIC 2017
H0 <- mean(H0.range)

# Save to lh_data
lh_data <- lapply(lh_data, function(x) { x$H0 = H0; x })

# Clean up
rm(H0, H0.range)

#-------------------------------------------------------------------------------
# VBGF fits
# Model

if(FALSE){
  
  H0 <- lh_data[[1]]$H0
  
  # Load Length-at-age data
  Barth_LD <- read.csv("data/Barth 2011 - Length at age.csv", header = TRUE)
  Block_LD <- read.csv("data/Block 2001 - Length at age.csv", header = TRUE)
  Block_LD <- Block_LD[Block_LD$Source != "Fortin et al 1993",]
  Bruch_LD <- read.csv("data/Bruch 2008 - Length at age.csv", header = TRUE)
  Fortin_LD <- read.csv("data/Fortin et al. 1992 - Length at age.csv", header = TRUE)
  Sunde_LD <- read.csv("data/Sunde 1961 - Length at age.csv", header = TRUE)
  Dieterman_LD <- read.csv("data/Dieterman 2010 - Length at age.csv", header = TRUE)
  
  # Combine data
  length_data <- as.data.frame(list(
    Source = c(rep("Barth", nrow(Barth_LD)),
               rep("Block", nrow(Block_LD)),
               rep("Bruch", nrow(Bruch_LD)),
               rep("Fortin", nrow(Fortin_LD)),
               rep("Sunde", nrow(Sunde_LD)),
               rep("Dieterman", nrow(Dieterman_LD))),
    Location = c(as.character(Barth_LD$Location),
                 as.character(Block_LD$Location),
                 as.character(Bruch_LD$Location),
                 as.character(Fortin_LD$Lake),
                 as.character(Sunde_LD$Location),
                 rep("Kettle River, MN", nrow(Dieterman_LD))),
    DU = c(Barth_LD$DU,
           Block_LD$DU,
           Bruch_LD$DU,
           Fortin_LD$DU,
           Sunde_LD$DU,
           Dieterman_LD$DU),
    Age = c(Barth_LD$Age,
            Block_LD$Age,
            Bruch_LD$Age,
            Fortin_LD$Age,
            Sunde_LD$Age,
            Dieterman_LD$Age),
    TL = c(Barth_LD$TL_mean_mm,
           Block_LD$TL_mm,
           Bruch_LD$TL_mean_mm,
           Fortin_LD$TL_mean_mm,
           Sunde_LD$TL_mm,
           Dieterman_LD$TL_mm)
  ))
  
  DUs <- c(1, 2, 2, 3, 4, 4)
  pred_data <- do.call(rbind, lapply(1:length(lh_data), function(x) {
    as.data.frame(list(
      DU = DUs[x],
      Curve = DU[x],
      Age = 0:100,
      TL =  with(lh_data[[x]], VB.growth(0:100, Linf, K, L0 = H0))
    ))
  }))
  
  
  p <- ggplot(length_data, aes(Age, TL))
  p <- p + geom_point(aes(colour = as.factor(DU)), size = 1, shape = 1)
  p <- p + geom_path(data = pred_data, aes(x = Age, y = TL, colour = as.factor(DU), group = Curve),
                     size = 1.0)
  p <- p + labs(x = "Age", y = "Length (mm)", colour = "DU") 
  p <- p + theme_me
  
  png("Figures/VB by DU with data.png", width = 4, height = 3, unit = "in", res = 300)
  p
  dev.off()
  
  # Barth VB fit
  Barth.L0.all.w <- nls(TL_mean_mm ~ Linf - (Linf - H0) * exp(-K * Age), data = Barth_LD,
                        start = c(Linf = 1000, K = 0.05), weights = n)
  
  Barth.t0.all.w <- nls(TL_mean_mm ~ Linf * (1 -  exp(-K * (Age-t0))), data = Barth_LD,
                        start = c(Linf = 1000, K = 0.05, t0 = 1), weights = n)
  with(Barth_LD, plot(Age, TL_mean_mm))
  xv = 0:28
  lines(xv, predict(Barth.L0.all.w, list(Age = xv)))
  lines(xv, predict(Barth.t0.all.w, list(Age = xv)), col = 2)
  
  # Block VB fit
  Block.L0.s <- nls(TL_mm ~ Linf - (Linf - H0) * exp(-K * Age), data = Block_LD,
                    start = c(Linf = 1500, K = 0.05))
  
  VBGF.f <- deriv(~ Linf ~ Linf - (Linf - H0) * exp(-K * Age), 
                  namevec = c('Linf', 'K'), function.arg=c('Age', 'Linf', 'K'))
  
  Block.L0.lme <- nlme(TL_mm ~ VBGF.f(Age, Linf, K), data = Block_LD, 
                       fixed = Linf + K ~ 1, random = Linf + K ~ 1|Location/Source,
                       na.action = na.omit,
                       start = list(fixed = coef(Block.L0.s)))
  
  with(Block_LD, plot(Age, TL_mm))
  
  Location <- unique(Block_LD$Location)
  Source <- unique(Block_LD$Source)
  counter <- 0
  Block_list <- NULL
  names <- NULL
  for(i in Source){
    for(j in Location) {
      subData <- Block_LD[Block_LD$Location == j & Block_LD$Source == i, ]
      if(dim(subData)[1] > 0){
        counter <- counter + 1
        Block_list[counter] <- list(subData)
        names <- c(names, paste0(i,"/",j))
      }
    }
  }
  names(Block_list) <- names
  
  #Clean up
  rm(Location, Source, counter, names, i, j, subData)
  
  # Model
  Block.L0.ind <- lapply(Block_list, function(x) nls(TL_mm ~ Linf - (Linf - H0) * exp(-K * Age), data = x,
                                                     start = c(Linf = 1500, K = 0.05)))
  Block.L0.coefs <- t(data.frame(lapply(Block.L0.ind, function(x) coef(x))))
  Block.L0.coefs <- as.data.frame(Block.L0.coefs[order(row.names(Block.L0.coefs)), ])
  Block.L0.coefs$t0 <- as.vector(apply(Block.L0.coefs, 1, function(x) agebyVB(0, x[[1]], x[[2]], L0 = H0)))
  #write.table(Block.L0.coefs, "Block VB.csv", sep = ",")
  
  # Bruch VB fit
  with(Bruch_LD, plot(Age, TL_mean_mm))
  
  Bruch.L0.w <- nls(TL_mean_mm ~ Linf - (Linf - H0) * exp(-K * Age), data = Bruch_LD,
                    start = c(Linf = 1500, K = 0.05), weights = N)
  lines(0:100, predict(Bruch.L0.w, list(Age = 0:100)))
  agebyVB(0, coef(Bruch.L0.w)[[1]], coef(Bruch.L0.w)[[2]], L0 = H0)
  
  # Fortin VB fit
  with(Fortin_LD, plot(Age, TL_mean_mm, col = Lake))
  Fortin_list <- lapply(unique(Fortin_LD$Lake), function(x) Fortin_LD[Fortin_LD$Lake == x, ])
  names(Fortin_list) <- unique(Fortin_LD$Lake)
  
  Fortin.L0.ind <- lapply(Fortin_list, function(x) nls(TL_mean_mm ~ Linf - (Linf - H0) * exp(-K * Age), 
                                                       data = x,
                                                       start = c(Linf = 1500, K = 0.05)))
  Fortin.L0.coefs <- t(data.frame(lapply(Fortin.L0.ind, function(x) coef(x))))
  lapply(1:3, function(x) {
    lines(1:35, lapply(Fortin.L0.ind, function(x) predict(x, list(Age = 1:35)))[[x]], col = x)
  })
  Fortin.L0.coefs <- as.data.frame(Fortin.L0.coefs[order(row.names(Fortin.L0.coefs)), ])
  Fortin.L0.coefs$t0 <- as.vector(apply(Fortin.L0.coefs, 1, function(x) agebyVB(0, x[[1]], x[[2]], L0 = H0)))
  
  # Sunde 1961
  with(Sunde_LD, plot(Age, TL_mm, col = Sex))
  
  Sunde.L0.w <- nls(TL_mm ~ Linf - (Linf - H0) * exp(-K * Age), data = Sunde_LD[Sunde_LD$Sex == "F", ],
                    start = c(Linf = 1500, K = 0.05), weights = n)
  lines(0:100, predict(Sunde.L0.w, list(Age = 0:100)))
  agebyVB(0, coef(Sunde.L0.w)[[1]], coef(Sunde.L0.w)[[2]], L0 = H0)
  
  # Dieterman 2010
  with(Dieterman_LD, plot(Age, TL_mm))
  
  Dieterman.L0.w <- nls(TL_mm ~ Linf - (Linf - H0) * exp(-K * Age), data =Dieterman_LD,
                        start = c(Linf = 1500, K = 0.05), weights = N)
  lines(0:100, predict(Dieterman.L0.w, list(Age = 0:100)))
  agebyVB(0, coef(Dieterman.L0.w)[[1]], coef(Dieterman.L0.w)[[2]], L0 = H0)
  
  # Clean up
  rm(xv, Barth.L0.all.w, Barth.t0.all.w, Barth_LD, Block.L0.s, VBGF.f, Block.L0.lme,
     Block_LD, Block.L0.ind, Block_list, Block.L0.coefs, Bruch_LD, Bruch.L0.w,
     Fortin_LD, Fortin_list, Fortin.L0.ind, Fortin.L0.coefs, Sunde_LD, Sunde.L0.w,
     Dieterman_LD, Dieterman.L0.w)
  
}

#-------------------------------------------------------------------------------
# Length and Latitude

meanLat <- round(mean(VBGF$Latitude), 2)
meanLon <- round(mean(VBGF$Longitude), 2)

VBGF$cLat <-  VBGF$Latitude - meanLat
VBGF$cLon <-  VBGF$Longitude - meanLon
L25_mod <- lm(L25 ~ cLat * cLon, data = VBGF)
summary(L25_mod)
vif(L25_mod)

lat <- seq(43, 55, 0.1); lon <- quantile(VBGF$Longitude, 0.5)
pred_data <- as.data.frame(list(
  Latitude = rep(lat, 2),
  L25 = c(predict(L25_mod, data.frame(cLat = lat - meanLat, cLon = lon - meanLon)),
          2569.6 - 49.1 * lat + 11.4 * lon),
  Model = c(rep("Fitted", length(lat)), rep("Fortin", length(lat)))
))

# Plot - 
p <- ggplot(VBGF, aes(Latitude, L25))
p <- p + geom_point(aes(colour = as.factor(DU)), size = 1)
p <- p + geom_path(data = pred_data, aes(x = Latitude, y = L25, linetype = Model), size = 1)
p <- p + scale_x_continuous(breaks = seq(43, 55, 2))
p <- p + scale_y_continuous(breaks = seq(700, 1900, 300))
p <- p + labs(x = "Latitude", y = "Length at Age-25 (mm)", colour = "DU") 
p <- p + theme_me +
  #theme(legend.position = "top")
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1), 
        legend.direction = "horizontal", legend.spacing.y = unit(-0.25, "cm"))


#png("Figures/Latitude - Length.png", width = 4, height = 3, unit = "in", res = 300)
p
#dev.off()

lat <- quantile(VBGF$Latitude, 0.5); lon <- seq(71, 111, 0.1)
pred_data <- as.data.frame(list(
  Longitude = rep(lon , 2),
  L25 = c(predict(L25_mod, data.frame(cLat = lat - meanLat, cLon = lon - meanLon)),
          2569.6 - 49.1 * lat + 11.4 * lon),
  Model = c(rep("Fitted", length(lon)), rep("Fortin", length(lon)))
))


predict(L25_mod, data.frame(cLat = 46.7 - meanLat, cLon = 84.45 - meanLon))
        
        
# Plot - 
q <- ggplot(VBGF, aes(Longitude, L25))
q <- q + geom_point(aes(colour = as.factor(DU)), size = 1)
q <- q + geom_path(data = pred_data, aes(x = Longitude, y = L25, linetype = Model), size = 1)
q <- q + scale_x_continuous(breaks = seq(71, 111, 5))
q <- q + scale_y_continuous(breaks = seq(700, 1900, 300))
q <- q + labs(x = "Longitude", y = "Length at Age-25 (mm)") 
q <- q + theme_me + theme(legend.position="none")

#png("Figures/Latitude - Length.png", width = 4, height = 3, unit = "in", res = 300)
q
#dev.off()

#png("Figures/Latitude and Longitude - Length.png", width = 6, height = 2.5, unit = "in", res = 300)
multiplot(p, q, cols = 2)
#dev.off()

# Test difference in L25 between DUs
DU_mod <- aov(L25~as.factor(DU), data = VBGF)
TukeyHSD(DU_mod)

# Plot - 
p <- ggplot(VBGF, aes(as.factor(DU), L25))
p <- p + geom_boxplot(fill = "#4271AE", alpha = 0.75)
p <- p + annotate("text", x = 2, y = 1175, label = "a")
p <- p + annotate("text", x = 1, y = 1350, label = "b")
p <- p + annotate("text", x = 3, y = 1685, label = "b")
#p <- p + coord_trans(y = "log")
p <- p + labs(x = "DU", y = "Length at Age-25 (mm)")               # Axis names
p <- p + theme_me 

#png("Figures/DU - Length.png", width = 4, height = 3, unit = "in", res = 300)
p
#dev.off()

# clean up
rm(lat, pred_data, DU_mod, p)

#-------------------------------------------------------------------------------
# VBGF Parameters - K ~ Linf

# Isolate VBGF data
VB_data <- VBGF[is.finite(VBGF$Linf) & is.finite(VBGF$K), ]
VB_data <- VB_data[VB_data$Source != "Magnin 1977 in Block 2001/my fit", ] # remove outlier

# Model K as a function of Linf
K_mod <- lm((K) ~ log(Linf) + log(L25), data = VB_data)
summary(K_mod)

Linf <- 1105:2090
pred_data <- as.data.frame(list(
  Linf = Linf,
  K = (predict(K_mod, data.frame(Linf = Linf, L25 = mean(VB_data$L25))))
))

# Plot - 
p <- ggplot(VB_data, aes(Linf, K))
p <- p + geom_point(aes(colour = as.factor(DU)), size = 1)
p <- p + geom_path(data = pred_data, aes(x = Linf, y = K), size = 1)
#p <- p + coord_trans(x = "log")
p <- p + labs(x = expression(bold(L[infinity])), y = "k", colour = "DU")               # Axis names
p <- p + theme_me 

#png("Figures/K - Linf.png", width = 3.5, height = 2.5, unit = "in", res = 300)
p
#dev.off()

# Clean uo
rm(Linf, pred_data, p, VB_data)

#-------------------------------------------------------------------------------
# VBGF 
# Predict VBGF for surgeon population location
# Predict L25, then optimize Linf by finding L25 with K related to Linf from above

# Expression
L25.ex <- expression( Linf - (Linf - H0) * exp(-K.ex * 25) )
K.ex <- expression( (coef(K_mod)[[1]] + coef(K_mod)[[2]] * log(Linf) + coef(K_mod)[[3]] * log(L25)) )

# Optimization function
Linf_optim_f <- function(Linf, target.L25, H0) {
  
  L.25 <- eval(as.expression(sapply(L25.ex, proc)), list(Linf = Linf, H0 = H0, L25 = target.L25))
  (L.25 - target.L25)^2
}

# Predict Length at age 25
L25.est <- (predict(L25_mod, data.frame(cLat = pos_data$Lat - meanLat, cLon = pos_data$Lon - meanLon)))

# Linf - Optimize for
Linf.est <- Map(function(L.25) {
  optimize(Linf_optim_f, c(500, 3500), tol = 1e-16, target.L25 = L.25, 
           H0 = lh_data[[1]]$H0)$minimum
}, L.25 = L25.est)

optimize(Linf_optim_f, c(500, 3500), tol = 1e-16, target.L25 = 1204.658, 
         H0 = lh_data[[1]]$H0)$minimum

# Predict K
K.est <- (predict(K_mod, data.frame(Linf = unlist(Linf.est), L25 = unlist(L25.est))))
predict(K_mod, data.frame(Linf = 1711.976, L25 = 1204.658))

# Save VBGF Parameters
lh_data <- lapply(1:6, function(x) { lh_data[[x]]$Linf = Linf.est[[x]]; lh_data[[x]] }) # Linf
lh_data <- lapply(1:6, function(x) { lh_data[[x]]$K = K.est[[x]]; lh_data[[x]] })       # K
lh_data <- lapply(lh_data, function(x) { x$cv.Li = 0.15; x})                            # Growth increment variability (constant proportion of mean) (Lorenzen 2016)

# Clean up
rm(L25.ex, K.ex, Linf_optim_f, L25.est, Linf.est, K.est)

#-------------------------------------------------------------------------------
# Length-Weight Relationship
#-------------------------------------------------------------------------------
# Fit relationship between length-weight regression coefs
# Then fin mean b and predict a

# Model relationship between a and b
LW_mod <- lm(log10(a) ~ log10(b), data = L_W)

# Predict values
b <- seq(2.5, 3.7, 0.05)
pred_data <- as.data.frame(list(
  a = 10^(predict(LW_mod, data.frame(b = b))),
  b = b
))

# Plot 
p <- ggplot(L_W, aes(b, a * 1e6))
p <- p + geom_point(aes(colour = as.factor(DU)), size = 1)
p <- p + geom_path(data = pred_data, aes(x = b, y = a * 1e6), size = 1)
p <- p + scale_y_continuous(trans = "log10", breaks = c(1e-4, 1e-3, 1e-2, 1e-1))
p <- p + scale_x_continuous(trans = "log10", breaks = seq(2.5, 3.5, 0.2))
p <- p + labs(y = bquote(a~(x~10^-6)), x = "b", colour = "DU")               # Axis names
#p <- p + labs(y = expression(paste(bold("a", 1e^-6))), x = "b", colour = "DU")               # Axis names

p <- p + theme_me 

#png("Figures/LW lm.png", width = 3.75, height = 2.75, unit = "in", res = 300)
p
#dev.off()

lh_data <- lapply(lh_data, function(x) { x$b.LW = mean(L_W$b); x})                             
lh_data <- lapply(lh_data, function(x) { x$a.LW = 10^(predict(LW_mod, data.frame(b = mean(L_W$b))))[[1]]; x}) 

# clean up
rm(b, pred_data, p, LW_mod)

# Examine lenght-weight by habitat type


subdata <- L_W[L_W$Habitat != "Both",]
subdata$W1000 <- with(subdata, a * 1000 ^ b)
subdata$Habitat <- factor(subdata$Habitat, levels = unique(subdata$Habitat))

anova(lm(W1000 ~ Habitat, data = subdata))
t.test(W1000~Habitat, data = subdata)

p <- ggplot(subdata, aes(Habitat, W1000))
p <- p + geom_violin(fill = "grey", alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75))
p <- p + labs(x = "Habitat Type", y = "Weight at 1000 mm (kg)")
p <- p + theme_me 

#png("Figures/LW  by habitat.png", width = 3.75, height = 2.75, unit = "in", res = 300)
p
#dev.off()

# no impact of habitat type weight at length

L <- seq(0, 2000, 10)
pred_W <- do.call(rbind, lapply(1:nrow(L_W), function(x){
  data <- L_W[x,]
  as.data.frame(list(
    L = L,
    W = with(data, a*L^b),
    DU = data$DU,
    Habitat = data$Habitat,
    row = x
  ))
}))

mean_W <- as.data.frame(list(
  L = L,
  W = lh_data[[1]]$a.LW * L ^ lh_data[[1]]$b.LW
))

# Plot 
p <- ggplot(pred_W, aes(L, W))
p <- p + geom_path(aes(colour = as.factor(Habitat), group = row), size = 1)
p <- p + geom_path(data = mean_W, aes(x = L, y = W), col = "black", size = 2)
p <- p + labs(x = "Length (mm)", y = "Weight (kg)", colour = "Habitat")
p <- p + theme_me + theme(legend.position = "top") 

#png("Figures/LW each relationship.png", width = 3.75, height =  3, unit = "in", res = 300)
p
#dev.off()


#-------------------------------------------------------------------------------
# Reproduction
#-------------------------------------------------------------------------------

# Periodicity - 2-7 years (COSEWIC 2017)
lh_data <- lapply(lh_data, function(x) { x$P =  mean(c(2,7)); x})                             

# Sex ratio - offspring sex
lh_data <- lapply(lh_data, function(x) { x$R =  0.5; x})    

# Maturity Ogive - Bruch 2008 - Lake Winnebago WI
MO_mod = nls(MO ~ 1.0/(1+exp(-b*(Age-c))), start = list(b = 0.35, c = 25), 
             data = MO_data[MO_data$Source == "Bruch 2008", ])
lh_data <- lapply(lh_data, function(x) { x$b.MO =   coef(MO_mod)[[1]]; x})    
lh_data <- lapply(lh_data, function(x) { x$c.MO =   round(coef(MO_mod)[[2]]); x})    

# Predict values
Age <- 5:50
pred_data <- as.data.frame(list(
  Age = Age,
  MO = predict(MO_mod, data.frame(Age = Age))
))

# Plot 
p <- ggplot(MO_data, aes(Age, MO))
p <- p + geom_point(aes(colour = Location), size = 1)
p <- p + geom_path(data = pred_data, aes(x = Age, y = MO), size = 1)
p <- p + labs(x = "Age", y = "Proportion Mature")               # Axis names
p <- p + theme_me + theme(legend.position = c(0.99,0.01), legend.justification = c(1,0))

#png("Figures/Maturity.png", width = 3, height = 2.75, unit = "in", res = 300)
p
#dev.off()

#clean up
rm(Age, pred_data, p)

# Stage Boundaries
lh_data <- lapply(lh_data, function(x) { x$j1.Tmax = floor((x$c.MO + 1) /2); x}) 
lh_data <- lapply(lh_data, function(x) { x$j2.Tmax = x$c.MO; x}) 
lh_data <- lapply(lh_data, function(x) { x$a1.Tmax = floor((x$Tmax + x$c.MO)/ 2); x}) 

# Lenth at maturity
Lmat <- lapply(lh_data, function(x){
  with(x, VB.growth(26, Linf, K, L0 = H0))
})

# Sex ratio
# Proportion male at age
# COSEWIC 2006
R_data <- as.data.frame(list(
  Age = c(0, 25, 40, 80),   # Age
  PM = c(0.5, 1/3, 1/7, 0) # Proportion Male
))

R <- NULL
R <- with(R_data, {
for(t in 0:80){
  if(t < Age[2]){
    s <- (PM[2] - PM[1])/(Age[2] - Age[1])
    R <- c(R, (t-Age[1])*s + PM[1])
  } else if(t < Age[3]){
    s <- (PM[3] - PM[2])/(Age[3] - Age[2])
    R <- c(R, (t-Age[2])*s + PM[2])
  } else if(t < Age[4]){
    s <- (PM[4] - PM[3])/(Age[4] - Age[3])
    R <- c(R, (t-Age[3])*s + PM[3])
  } else if(t >= Age[4]){
    R <- c(R, 0)
  }
}
R
})

# Predict values
Age <- 0:80
pred_data <- as.data.frame(list(
  Age = 0:80,
  R = R
))

# Plot 
p <- ggplot(R_data, aes(Age, PM))
p <- p + geom_point(aes(), size = 1.5)
p <- p + geom_path(data = pred_data, aes(x = Age, y = R), size = 1)
p <- p + labs(x = "Age", y = "Proportion Male")               # Axis names
p <- p + theme_me 

#png("Figures/Sex Ratio.png", width = 3, height = 2.75, unit = "in", res = 300)
p
#dev.off()

rm(R_data, Age, pred_data, p)

#-------------------------------------------------------------------------------
# Fecundity Estimates 
#-------------------------------------------------------------------------------
# 
# F1 <- 18683.6 * W - 183527 # Bruch 2008
# F2 <- 13270.96 * W - 36937 # COSEWIC 2017
# F3 <- 16640 * W - 150683   # Bruch et al. 2006
# F4 <- 12000 * W            # Velez-espino and Koops 2010

# COSEWIC 2017
Fec_mod <- lm(log(Eggs) ~ log(W_kg), data = Fec_data)
summary(Fec_mod)

# Predict values
W <- seq(5, 85, 0.25)
pred_data <- as.data.frame(list(
  W = W,
  Fec = exp(coef(Fec_mod)[[1]]) * W ^ coef(Fec_mod)[[2]]
))

# Plot 
p <- ggplot(Fec_data, aes(W_kg, Eggs))
p <- p + geom_point(aes(), size = 1)
p <- p + geom_path(data = pred_data, aes(x = W, y = Fec), size = 1)
p <- p + labs(y = "Egg Count", x = "Weight (kg)")               # Axis names
p <- p + theme_me 

#png("Figures/Fecundity.png", width = 4, height = 3, unit = "in", res = 300)
p
#dev.off()

# Parameters - COSEWIC 2017
lh_data <- lapply(lh_data, function(x) { x$a.F = exp(coef(Fec_mod)[[1]]); x})    
lh_data <- lapply(lh_data, function(x) { x$b.F = coef(Fec_mod)[[2]]; x})    

# clean up
rm(p, W, pred_data, t, L, W, F_pred)

#-------------------------------------------------------------------------------
# Life History Plot
#-------------------------------------------------------------------------------

# LENGTH
t <- 0:80
VB_pred <- as.data.frame(list(
  Age = rep(t, 6),
  Y = unlist(lapply(lh_data, function(x) x$Linf - (x$Linf - x$H0) * exp(-x$K * t))),
  DU = c(rep(1, length(t)), rep("2a", length(t)), rep("2b", length(t)),
         rep(3, length(t)), rep("4a", length(t)), rep("4b", length(t))),
  Var = "Length (mm)"
))

# WEIGHT
t <- 0:80
L <- lapply(lh_data, function(x) x$Linf - (x$Linf - x$H0) * exp(-x$K * t))

W_pred <- as.data.frame(list(
  Age = rep(t, 6),
  Y = unlist(lapply(1:6, function(x) lh_data[[x]]$a.LW * L[[x]] ^ lh_data[[x]]$b.LW)),
  DU = c(rep(1, length(t)), rep("2a", length(t)), rep("2b", length(t)),
         rep(3, length(t)), rep("4a", length(t)), rep("4b", length(t))),
  Var = "Weight (kg)"
))

# FECUNDITY
t <- 19:80
L <- lapply(lh_data, function(x) x$Linf - (x$Linf - x$H0) * exp(-x$K * t))
W <- lapply(1:6, function(x) lh_data[[x]]$a.LW * L[[x]] ^ lh_data[[x]]$b.LW)

F_pred <- as.data.frame(list(
  Age = rep(t, 6),
  Y = unlist(lapply(1:6, function(x) lh_data[[x]]$a.F * W[[x]] ^ lh_data[[x]]$b.F)),
  DU = c(rep(1, length(t)), rep("2a", length(t)), rep("2b", length(t)),
         rep(3, length(t)), rep("4a", length(t)), rep("4b", length(t))),
  Var = "Fecundity (eggs/female)"
))

# PLOT
p <- ggplot(rbind(VB_pred, W_pred, F_pred), aes(Age, Y))
p <- p + geom_path(aes(colour = as.factor(DU)), size = 1)
p <- p + facet_wrap(~Var, scales = "free_y", ncol = 1)
p <- p + labs(y = "", x = "Age", colour = "DU")               # Axis names
p <- p + theme_me 

#png("Figures/Life History.png", width = 5.5, height = 6.5, unit = "in", res = 300)
p
#dev.off()

# Clean up
rm(p, F_pred, W_pred, VB_pred, t, L, W)

#-------------------------------------------------------------------------------
# SURVIVAL
#-------------------------------------------------------------------------------
# Mortality estimates are avaialbe from lake sturgeon population
# Either use M values as a guide to estiamte area-specifinc M from a literature
# relationship for take teh mid-point of M estiamtes and apply accross locations
# Estiamted of juvenile survival are also avaiable - use as a guide for estiamtes
# form Lorrencen (2000) length base estimator. 
# YOY estimated are availble - ue a guide when optimizing for YOY survival

# Plot - box plot comparing M by estimation method and DU
p <- ggplot(M_data[is.na(M_data$Juvenile) & is.na(M_data$F), ], aes(as.factor(Method), M))
p <- p + geom_boxplot(aes(fill = as.factor(DU)), alpha = 0.75)
#p <- p + coord_trans(y = "log")
p <- p + labs(x = NULL, y = "Intantaneous Mortality", fill = "DU")               # Axis names
p <- p + theme_me 

#png("Figures/DU - Mortality.png", width = 3.5, height = 2.5, unit = "in", res = 300)
p
#dev.off()

with(M_data[is.na(M_data$Juvenile) & is.na(M_data$F), ], aggregate(M, by = list(Method), mean))
with(M_data[is.na(M_data$Juvenile) & is.na(M_data$F), ], aggregate(M, by = list(Method), sd))

# Median annual mortality from adult non-commercially fished population

#lh_data <- lapply(1:6, function(x) { lh_data[[x]]$M = lh_data[[x]]$K * 1.5; lh_data[[x]]}) # Beverton and Holt
lh_data <- lapply(1:6, function(x) { lh_data[[x]]$M = 4.1 * lh_data[[x]]$K ^ 0.7 * lh_data[[x]]$Linf ^ -0.3; lh_data[[x]]}) # Then et al. 2015
#lh_data <- lapply(1:6, function(x) { lh_data[[x]]$M = with(M_data[is.na(M_data$Juvenile) & is.na(M_data$F), ], median(M)); lh_data[[x]]}) 
lh_data <- lapply(lh_data, function(x) { x$cv.M = 0.2; x}) 
#lh_data <- lapply(lh_data, function(x) { x$syoy = NA; x})    

# Estimate m0 - mortality at unit length 
t <- round(lh_data[[1]]$c.MO)                                               # Age at 50% maturity - from maturity ogive
La <- lapply(lh_data, function(x) x$Linf - (x$Linf - x$H0) * exp(-x$K * t)) # Length at age-of-maturity
m0.list <- lapply(1:6, function(x) La[[x]] * lh_data[[x]]$M)                # Mortality at unit length

#lh_data <- lapply(1:6, function(x) c(lh_data[[x]], m0 = m0.list[x]))

# PLot Measured survival at age and estiamted survival at age (Lorenzen 2000)
t = 1:65
L.list <- lapply(lh_data, function(x) x$Linf - (x$Linf - x$H0) * exp(-x$K * t))

s <- lapply(1:6, function(x){
       data = lh_data[[x]]
       L = with(data, VB.growth(t, Linf, K, L0 = H0))
       ifelse(t < data$c.MO, (L[1:(length(L)-1)] / L[2:(length(L))] * exp(-data$K)) ^ (m0.list[[x]] / data$K / data$Linf),
              exp(-data$M))
    })
Ms <- lapply(s, function(x) -log(x))
M.tot <- lapply(Ms, sum)

pred_data <- as.data.frame(list(
  age = rep(t, 6),
  M = unlist(Ms),
  DU = sort(rep(DU, (length(t))))
))

# Plot - Instantaneous mortality
col <- ggplot_colour(3)[c(1,3)]
p <- ggplot(M_data[is.finite(M_data$Age_min) & is.na(M_data$F), ])
p <- p + geom_path(data = pred_data, aes(x = age, y = M, linetype = as.factor(DU)))
p <- p + geom_segment(aes(x = Age_min, xend = Age_max, y = M, yend = M, colour = as.factor(DU)), 
                      arrow = arrow(angle = 90, ends = "both", length = unit(0.075, "inches")))
#p <- p + geom_point(aes(x = (Age_min + Age_max)/2, y = M, colour = as.factor(DU)))
p <- p + scale_x_continuous(breaks = seq(10,60, 10))
#p <- p + scale_colour_manual(values = col)
p <- p + labs(x = "Age", y = "Instantaeous Mortality", colour = "DU", linetype = "DU")               # Axis names
p <- p + theme_me + theme(legend.position = "none")

# Plot - annual surivval
col <- ggplot_colour(3)[c(1,3)]
q <- ggplot(M_data[is.finite(M_data$Age_min) & is.na(M_data$F), ])
q <- q + geom_path(data = pred_data, aes(x = age, y = exp(-M), linetype = as.factor(DU)))
q <- q + geom_segment(aes(x = Age_min, xend = Age_max, y = exp(-M), yend = exp(-M), colour = as.factor(DU)), 
                      arrow = arrow(angle = 90, ends = "both", length = unit(0.075, "inches")))
#q <- q + geom_point(aes(x = (Age_min + Age_max)/2, y = exp(-M), colour = as.factor(DU)))
q <- q + scale_x_continuous(breaks = seq(10,60, 10))
#p <- p + scale_colour_manual(values = col)
q <- q + labs(x = "Age", y = "Annual Survival", colour = "DU", linetype = "DU")               # Axis names
q <- q + theme_me +
  theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0), legend.box= "horizontal")

#png("Figures/Age vs Mortality.png", width = 8, height = 3, unit = "in", res = 300)
multiplot(p, q, cols = 2)
#dev.off()


# Clean uo
rm(p, t, La, m0.list, s, Ms, M.tot, pred_data, col)


if(FALSE){

M.EtoL <- c(2.813, 4.711, 8.517-4.05, 7.264-2.34)
S.EtoL <- exp(-M.EtoL)
DU <- c(4, 4, 4, 4)

M.Lto0 <- c(4.05, 2.34)
S.Lto0 <- exp(-M.Lto0)
DU <- c(4, 4)

M.Eto0 <- c(8.517, 7.264)
S.Eto0 <- exp(-M.Eto0)
DU <- c(4, 4)

M.0_winter <- -log(0.4)
S.0_winter <- 0.4

M_yoy <- M.Eto0 + M.0_winter
S_yoy <- exp(-M_yoy)

# Recruitment
S = seq(0, 10000, 100)
R = 20.07706*S*exp(-0.000445*S)/2 # For female Age-1 recruite by female spawner
plot(S/668, R)

S = seq(0, 10000, 100)
R = 2.99395*S*exp(-0.000436*S)  # For female Age-21 recruite by female spawner
plot(S, R)


Sa <- exp(-1.65*0.042)
s1a = Sa*1/24
s1b = Sa*14/24

mean(c(s1a, s1b))

s2a = Sa*14/24
s2b = Sa*24/24

mean(c(s2a, s2b))

}

#===============================================================================
# MODEL EXPRESSION
#-------------------------------------------------------------------------------
# R - sex ratio
# MO - Maturity Ogive
# P - Periodicity
# Ft - fecundity
# st - Survival

mat_ogive <- function(t, b, c, names = TRUE) {
  MO <- round(1.0 / (1 + exp(-b * (t - c))), 2)
  if(names == TRUE) names(MO) <- paste0("MO", t)
  MO
}

fec_est <- function(W, a, b, names = TRUE){
  f <- a * W ^ b 
  if(names == TRUE) names(f) <- paste0("f", 0:(length(W) - 1))
  f
}

surv_est <- function(t, data, names = TRUE){
  
  s <- with(data, 
            ifelse(t < 1, s0, 
                   ifelse(t >= c.MO, exp(-M), {
                     m0 <-  VB.growth(c.MO, Linf, K, L0 = H0) * M
                     (VB.growth(t, Linf, K, L0 = H0) / VB.growth(t + 1, Linf, K, L0 = H0) * exp(-K)) ^ 
                       (m0 / K / Linf)
                   }
                   )))
  if(names == TRUE) names(s) <- paste0("s", t)
  s
}


Fertility <- function(data, f, MO, names = TRUE) {
  Ft <- data$R * MO * f / data$P
  if(names == TRUE) names(Ft) <- paste0("F", 0:(length(f) - 1))
  Ft
}

# Stage transition
Gj1 <- expression( ((sj1/lambda) ^ (j1.Tmax - 1) - (sj1/lambda) ^ (j1.Tmax - 1 - 1)) / 
                     ((sj1/lambda) ^ (j1.Tmax - 1) - 1) )
Gj2 <- expression( ((sj2/lambda) ^ (j2.Tmax - j1.Tmax) - (sj2/lambda) ^ (j2.Tmax - j1.Tmax - 1)) / 
                     ((sj2/lambda) ^ (j2.Tmax - j1.Tmax) - 1) )
Ga1 <- expression( ((sa1/lambda) ^ (a1.Tmax - j2.Tmax) - (sa1/lambda) ^ (a1.Tmax - j2.Tmax - 1)) / 
                     ((sa1/lambda) ^ (a1.Tmax - j2.Tmax) - 1) )
Ga2 <- expression( ((sa2/lambda) ^ (Tmax - a1.Tmax) - (sa2/lambda) ^ (Tmax - a1.Tmax - 1)) / 
                     ((sa2/lambda) ^ (Tmax - a1.Tmax) - 1) )


#Gj1 <- expression( 1 / (j1.Tmax - 1) )
#Gj2 <- expression( 1 / (j2.Tmax - j1.Tmax) )
#Ga1 <- expression( 1 / (a1.Tmax - j2.Tmax) )
#Ga2 <- expression( 1 / (Tmax - a1.Tmax) )

#===============================================================================
# Define the matrix structure 
#-------------------------------------------------------------------------------
# Projection matrix structure of a single population
# Matrix structure - Vital rates

A.LS <- expression(
    0,        0,         Fj*sj2*(1-Gj2) + F1*sj2*Gj2, F1*sa1*(1-Ga1) + F2*sa1*Ga1,  F2*sa2*(1-Ga2),
  syoy, sj1 * (1 - Gj1),                0,                           0,                  0,
    0,     sj1 * Gj1,             sj2 * (1 - Gj2),                   0,                  0,
    0,        0,                    sj2 * Gj2,                 sa1 * (1 - Ga1),          0,
    0,        0,                        0,                       sa1 * Ga1,        sa2 * (1 - Ga2))
  
# Matrix struction - lower level function with M instead of m0 (m0 = M * L3) 
mx.LS <- as.expression(sapply(A.LS, proc))  # proc() substitutes sub expressions in (see MatrixFunctions.r)

#-------------------------------------------------------------------------------
# Minimum Lambda vital rates
#-------------------------------------------------------------------------------

# Params used to reduce lambda to lambda min
harm <- c("equil", "M0", "Mj", "Ma", "Mf")
l_min <- lapply(harm, function(x) list())
names(l_min) <- harm

l_min$equil <- lapply(lh_data, function(x) {
  list(M0 = 1,  # YOY mortality increase
       Mj = 1, # Young juvenile mortality increase
       Ma = 1, # Young adult mortality increase
       Mf = 1)  # Fecundity reduction
})

#clean up
rm(harm)

# function to generate mean projection matrces
pmx_create <- function(data) {
  
  Lt <- VB.growth(0:data$Tmax, data$Linf, data$K, L0 = data$H0, names = TRUE)
  Wt <-  LW_fun(Lt, data)
  ss <- surv_est(0:data$Tmax, c(data))
  MO <-  mat_ogive(0:data$Tmax, b = data$b.MO, c = data$c.MO, names = TRUE)
  fs <- with(data, fec_est(Wt, a = a.F, b.F, names = TRUE))
  Fs <- Fertility(data, fs, MO)
  
  s_stage <- with(data, list(
    syoy = ss[1],
    sj1 = exp(-mean(-log(ss[2:j1.Tmax]))),
    sj2 = exp(-mean(-log(ss[(j1.Tmax+1):j2.Tmax]))),
    sa1 = exp(-mean(-log(ss[(j2.Tmax+1):a1.Tmax]))),
    sa2 = exp(-mean(-log(ss[(a1.Tmax+1):(Tmax+1)])))
  ))
  
  F_stage <- {
    
    # Age-specific mortality rate
    M <- -log(ss)
    
    # stage stage dist'n within stage - use to weight fecundity est.
    wj <- exp(-cumsum(M[(data$j1.Tmax + 1):data$j2.Tmax]))
    w1 <- exp(-cumsum(M[(data$j2.Tmax + 1):data$a1.Tmax]))
    w2 <- exp(-cumsum(M[(data$a1.Tmax + 1):data$Tmax]))
    
    # aggregate Fertility to stage estimat
    list(
      Fj = sum(Fs[(data$j1.Tmax + 1):data$c.MO] * wj) / sum(wj), # Pre adult - first half of maturity ogive
      F1 = sum(Fs[(data$c.MO + 1):data$a1.Tmax] * w1) / sum(w1), # young adult - Tmat to half adult stage
      F2 = sum(Fs[(data$a1.Tmax+1):(data$Tmax)] * w2) / sum(w2)  # old adult - ent young adult to Tmax
    )
  }
  
  # Mean annual stage transition
  G_stage <- c(Ex_means(prefix = "Gj", suffix = 1:2, data = c(data, s_stage, lambda = 1)),
               Ex_means(prefix = "Ga", suffix = 1:2, data = c(data, s_stage, lambda = 1)))
  
  # If lambda min then reduce vital rates by x
  if(any(with(data, c(M0, Mj, Ma, Mf)) < 1)){
    if(data$M0 < 1) {
      s_stage <- {
        s_stage$syoy = s_stage$syoy * data$M0 
        s_stage
      }

    } else if(data$Mj < 1) { 
      s_stage <- {
          s_stage$sj1 = s_stage$sj1 * data$Mj 
          s_stage$sj2 = s_stage$sj2 * data$Mj 
          s_stage
      }
    } else if(data$Ma < 1) { 
      s_stage <- {
          s_stage$sa1 = s_stage$sa1 * data$Ma 
          s_stage$sa2 = s_stage$sa2 * data$Ma 
          s_stage
      }
    } else if(data$Mf < 1) {  # Fecundity reduction
      F_stage <- {
          F_stage$Fj = F_stage$Fj * data$Mf 
          F_stage$F1 = F_stage$F1 * data$Mf 
          F_stage$F2 = F_stage$F2 * data$Mf 
          F_stage
      }
    }
  }
  
  # evaluate projection matrix
  pmxs <- pmx_eval(A.LS, c(data, F_stage, s_stage, G_stage))
  
  # output
  pmxs
  
}

#-------------------------------------------------------------------------------
# Variability
#-------------------------------------------------------------------------------

# Initial length at age
L_mean <- lapply(lh_data, function(x) VB.growth(0:x$Tmax, x$Linf, x$K, L0 = x$H0, names = TRUE))

# Mean annual growth increments by age
Li_means <- lapply(1:6, function(x) L_mean[[x]][2:(lh_data[[x]]$Tmax + 1)] - L_mean[[x]][1:lh_data[[x]]$Tmax])

# Mean annual Survival 
s_means <- lapply(lh_data, function(x) surv_est(0:(x$Tmax - 1), c(x, s0 = NA)))

# Mean annual Instantanesou Mortality
M_means <- lapply(s_means, function(x){
  Ms <- -log(x)
  names(Ms) <- paste0("M", 0:(length(x) - 1))
  Ms
})

# names for  vectors
names_Li <- lapply(1:6, function(x) paste0("Li", 0:(lh_data[[x]]$Tmax - 1))) # Lengths
names_s <- lapply(1:6, function(x) paste0("s", 0:(lh_data[[x]]$Tmax - 1)))   # survival

# Correlation structure
cor_Li <- lapply(lh_data, function(x) matrix(1, ncol = x$Tmax, nrow = x$Tmax))
#cor_Li <- lapply(lh_data, function(x) {mat <- diag(x$Tmax); mat[mat == 0] <- 0.5; mat })
#cor_Li <- lapply(lh_data, function(x) AR1(n = x$Tmax, rho = 0.75))
cor_s <- lapply(1:6, function(x) {
  mat <- matrix(1, ncol = length(s_means[[x]]), nrow = length(s_means[[x]]))
  mat[, 1] <- c(1, rep(0, (length(s_means[[x]]) - 1)))
  mat[1, ] <- c(1, rep(0, (length(s_means[[x]]) - 1)))
  mat
})

# Distributions for each random variables
dist_Li <- lapply(lh_data, function(x) rep("beta_stretch", x$Tmax))
dist_s <- lapply(1:6, function(x) rep("beta_stretch", length(s_means[[x]])))

# parameter values for random variables distributions
params_Li <- lapply(1:6, function(x) { # Growth increments
  lapply(1:lh_data[[x]]$Tmax, function(i) beta_stretch_val(Li_means[[x]][[i]], sd = Li_means[[x]][[i]] * lh_data[[x]]$cv.Li))
})

params_s <- lapply(1:6, function(x) { # Survival
  lapply(1:(lh_data[[x]]$Tmax), function(i) beta_stretch_val(M_means[[x]][[i]], sd = (M_means[[x]][[i]] * lh_data[[x]]$cv.M)))
})

#-------------------------------------------------------------------------------
# YOY SURVIVAL - lambda = 1
#-------------------------------------------------------------------------------
# Optimize for  estiamtes of YOY survival that produce matrices with geometric mean
# lambda values of 1, mean and max lambda

# Deterministic
# no variability - to give pmx with lambda = 1
M_optim_f <- function(syoy, target.lambda, data){
  
  pmx <- pmx_create(c(data, s0 = syoy))
  lambda.est <- lambda(pmx)
  
  # est distance from trarget lambda
  (lambda.est - target.lambda)^2
}
M_optim <- lapply(1:6, function(x){
  optimize(M_optim_f, c(0, 1e-1), tol = 1e-16, target.lambda = 1, 
           data = c(lh_data[[x]], l_min$equil[[x]], lambda = 1))$minimum
})

lh_data <- lapply(1:6, function(x) { lh_data[[x]]$syoy.1d = M_optim[[x]]; lh_data[[x]]})    

# clean up
rm(M_optim_f, M_optim)

#-------------------------------------------------------------------------------
# Stochastic syoy

if(FALSE){
# Function to optimize
M_optim_f = function (s0, target.lambda, numReps, years, pmxs, data){
  
  # mean projection matrix
  pmx_mean <- pmx_create(c(data, s0 = s0))
  
  # initial population size
  N0 <- initPop(pmx_mean, 100, "adult", p.rep = c(0, 0, 0, 1, 1))
  
  # Mean annual survival rates by age
  M0 <- -log(s0)
  
  # Vital rates
  params.s0 <- beta_stretch_val(M0, sd = (M0 * data$cv.M))
  
  # Generate random YOY survival rate
  s0.rand <- exp(-with(params.s0, rbeta(years * numReps, shape1 = a, shape2 = b) * (max - min) + min))
  s0.rand <- matrix(s0.rand, nrow = numReps, ncol = years)
  
  # insert s.yoy into projection matrices
  pmx_list <- lapply(1:numReps, function (n){ lapply(1:years, function(y) {
    pmxs[[n]][[y]][2,1] = s0.rand[n, y]; pmxs[[n]][[y]]
  })})
  
  # est. geom_mean of lambda of pmx.list
  lambdas <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$longterm) 
  #lambdas <- sapply(pmx_list, lambda_mean)
  mean_lambda <- mean(lambdas)
  
  # est distance from trarget lambda
  (mean_lambda - target.lambda)^2
}

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores) 
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))

M_optim <- parLapply(cl, 1:6, function(x){
  #cat("Location", x, "Of 6", "\n")
  
  years = 100
  numReps = 500
  
  # Generate random variables
  pmx_rand <- pmx_build_runs(mx = A.LS,     # Projection matrix expression
                 numReps = numReps,         # number of simulation replicates
                 years = years,             # number of years to simulate
                 lambda = 1,                # Pop growth rate
                 cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                 names.Li = names_Li[[x]],  # names for growth increments
                 dist.Li = dist_Li[[x]],    # distributions names for random length increments
                 params.Li = params_Li[[x]],# distributions parameters fore random length increments
                 cor.s = cor_s[[x]],        # Correlation struction for vital rates
                 names.s = names_s[[x]],    # names for vital rates
                 dist.s = dist_s[[x]],      # distributions names for random  vital rates
                 params.s = params_s[[x]],  # distributions parameters fore random vital rate
                 L.mean = L_mean[[x]],      # Initial (mean) length at age
                 data = c(lh_data[[x]], l_min$equil[[x]], s0 = NA)) # Life history data 
  
  pmx_init <- pmx_rand$pmx
  
  rm(pmx_rand)
  gc()
  
  # Solve for syoy with lambda = 1
  equil <- optimize(M_optim_f, c(0, 1e-1), tol = 1e-16, target.lambda = 1,
                    numReps = numReps,         # number of simulation replicates
                    years = years,             # number of years to simulate
                    pmxs = pmx_init,           # random projection matrcies
                    data = c(lh_data[[x]], l_min$equil[[x]])) # Life history data 
  
  list(equil = equil$minimum)
})

stopCluster(cl)

# Assign syoy to lh_data
lh_data <- lapply(1:6, function(x) { 
  lh_data[[x]]$syoy.1 = M_optim[[x]]$equil   # lambda = 1
  lh_data[[x]]
}) 

rm(years, numReps, M_optim, M_optim_f)
}

# Set YOY survival
Ms.1 <- c(0.0001186, 0.0000326, 0.0000049, 0.0000127, 0.0000145, 0.0000391)

lh_data <- lapply(1:6, function(x) { lh_data[[x]]$syoy.1 = Ms.1[x]; lh_data[[x]]})    

# Mean annual Mortality
params_s.1 <- lapply(1:6, function(x) { # Survival
  
  # Mean annual Survival 
  s_means <- surv_est(0:(lh_data[[x]]$Tmax - 1), c(lh_data[[x]], s0 = lh_data[[x]]$syoy.1))
  
  # Mean annual Instantanesou Mortality
  M_means <- {
    Ms <- -log(s_means)
    names(Ms) <- paste0("M", 0:(length(s_means) - 1))
    Ms
  }
  
  lapply(M_means, function(i) beta_stretch_val(i, sd = (i * lh_data[[x]]$cv.M)))
})

#-------------------------------------------------------------------------------
# Projection matrix
#-------------------------------------------------------------------------------
# Evaluate porjection matrix

pmx <- lapply(1:6, function(x) pmx_create(c(lh_data[[x]], s0 = lh_data[[x]]$syoy.1, l_min$equil[[x]])))

# Summary data
lambda.est <- lapply(pmx, lambda)                      # Lambda - pop growth rate (estimated M)
generation_time <- lapply(pmx, generation.time)        # Generation time
stable_stage <- lapply(pmx, stable.stage)              # Stabel stage structure
reproductive_value <- lapply(pmx, reproductive.value)  # Reproductive value

# Deterministic matrices
pmx.d <- lapply(1:6, function(x) pmx_create(c(lh_data[[x]], s0 = lh_data[[x]]$syoy.1d, l_min$equil[[x]])))
lapply(pmx.d, lambda)

#-------------------------------------------------------------------------------
# POPULATION GROWTH RATE (Lambda)
#-------------------------------------------------------------------------------
# 

# Dieterman et al. (2010)
lambda <- c(NA,  2.07, 1.24, 0.88, 1.40, NA, 1.11, 0.87, 1.16, 1.09, 0.78, 0.47, NA)
geom_mean(lambda, na.rm = TRUE) # ~ 1.035 - good representation of long term growth rate of recovering pop
r <- log(lambda)

# McDougall et al. (2017)
Na1 <- c(772, 927, 1342, 1252, 1644, 1973, 1988, 1867)
Na2 <- c(5794, 9919, 8061, 6600, 5384, 4391, 5005)
geom_mean(Na1[2:length(Na1)]/Na1[1:(length(Na1)-1)])
geom_mean(Na2[2:length(Na2)]/Na2[1:(length(Na2)-1)])

#year <- c(1993:1998, 2000:2006)
#N <- c(170, 154, 163, 169, 131, 205, 208, 268, 237, 271, 294, 210, 119)

# Haxton et al. (2014)
r <- c( 0.16, 0.13, 0.084, 0.146, 0.1, 0.118, 0.125)
exp((r))

# Bruch 2008
# r <- 0.016, 0.086, 0.49
# Model: Chapman, Jolly Seber, SCAA

# Maximum Lambda - Randall and Minn (2000) - estiamte
lambda.max <- lapply(lh_data, function(x){
  Lmat <- with(x, VB.growth(c.MO, Linf = Linf, K = K, L0 = H0)) # LEngth at maturity
  Wmat <- x$a.LW * Lmat ^ x$b.LW                                # Weight at maturity  
  exp(2.64 * (Wmat*1000) ^ -0.35)                               # Lambda max
})

# Minimum lambda - COSEWIC: 70% decline ofer 10 years or 3 generations for endangered status
decline <- list(DU1 = 0.98,
                DU2a = 0.9,
                DU2b = 0.9,
                DU3 = 0.25,
                DU4a = 0.99,
                DU4b = 0.99)

lambda.min <- lapply(1:6, function(x) {
  
  time <- 3 * generation_time[[x]]
  if(time < 10) { time <- 10 }

  if(time > 100){ time <- 100 }
  
  l = (1 - decline[[x]]) ^ (1 / time)
  
  #if(lh_data[[x]]$Status == "Endangered") {            # 70% decline over 10 years
  #  l = 0.3 ^ (1 / time)
  #} else if(lh_data[[x]]$Status == "Threatened"){      # 50% decline over 10 years
  #  l = 0.5 ^ (1 / time)
  #} else if(lh_data[[x]]$Status == "Special Concern"){ # 25% decline over 10 years
  #  l = 0.75 ^ (1 / time)
  #}
  l
})

# Mean Lambda - geometric mean of min, max and 1 lambdas
lambda.mean <- lapply(1:6, function(x){
  geom_mean(lambda, na.rm = TRUE) 
})

# Clean up
rm(lambda, r)

# COSEWIC STATUS - 6 month interim report
lambda.exp <- list(
  DU1 = (1 - 0.2) ^ (1 / 10),                                  # DU 1:  ~20% decline over 10 year
  DU2a = mean(c((1 + 0.1) ^ (1 / 10), (1 + 0.5) ^ (1 / 10))),  # DU 2: 10:50% increase projeced over next 10 years
  DU2b = mean(c((1 + 0.1) ^ (1 / 10), (1 + 0.5) ^ (1 / 10))),  # DU 2: 10:50% increase projeced over next 10 years
  DU3 = NA,                                                    # DU 3: Unknown
  DU4a = mean(c((1 + 0.05) ^ (1 / 10), (1 + 0.2) ^ (1 / 10))), # DU 4:5-20% increase projected over next 10 years
  DU4b = mean(c((1 + 0.05) ^ (1 / 10), (1 + 0.2) ^ (1 / 10)))  # DU 4:5-20% increase projected over next 10 years
)

lambda.1 <- lapply(1:6, function(x){
  1
})

#-------------------------------------------------------------------------------
# YOY SURVIVAL - mean and max lambda
#-------------------------------------------------------------------------------
# Optimize for  estiamtes of YOY survival that produce matrices with geometric mean
# lambda values of 1, mean and max lambda

if(FALSE){
  # Function to optimize
  M_optim_f = function (s0, target.lambda, numReps, years, pmxs, data){
    
    # mean projection matrix
    pmx_mean <- pmx_create(c(data, s0 = s0))
    
    # initial population size
    N0 <- initPop(pmx_mean, 100, "adult", p.rep = c(0, 0, 0, 1, 1))
    
    # Mean annual survival rates by age
    M0 <- -log(s0)
    
    # Vital rates
    params.s0 <- beta_stretch_val(M0, sd = (M0 * data$cv.M))
    
    # Generate random YOY survival rate
    s0.rand <- exp(-with(params.s0, rbeta(years * numReps, shape1 = a, shape2 = b) * (max - min) + min))
    s0.rand <- matrix(s0.rand, nrow = numReps, ncol = years)
    
    # insert s.yoy into projection matrices
    pmx_list <- lapply(1:numReps, function (n){
      lapply(1:years, function(y) {pmxs[[n]][[y]][2,1] = s0.rand[n, y]; pmxs[[n]][[y]]})
    })
    
    # est. geom_mean of lambda of pmx.list
    lambdas <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$longterm) 
    #lambdas <- sapply(pmx_list, lambda_mean)
    mean_lambda <- mean(lambdas)
    
    
    # est distance from trarget lambda
    (mean_lambda - target.lambda)^2
  }
  
  no_cores <- detectCores() - 1
  
  cl <- makeCluster(no_cores) 
  clusterExport(cl, ls())
  clusterEvalQ(cl, library(popbio))
  
  M_optim <- parLapply(cl, 1:6, function(x){
    
    years = 100
    numReps = 500
    
    # Generate random variables
    pmx_rand <- pmx_build_runs(mx = A.LS,     # Projection matrix expression
                               numReps = numReps,         # number of simulation replicates
                               years = years,             # number of years to simulate
                               lambda = lambda.max[[x]],   # pop growth rates
                               cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                               names.Li = names_Li[[x]],  # names for growth increments
                               dist.Li = dist_Li[[x]],    # distributions names for random length increments
                               params.Li = params_Li[[x]],# distributions parameters fore random length increments
                               cor.s = cor_s[[x]],        # Correlation struction for vital rates
                               names.s = names_s[[x]],    # names for vital rates
                               dist.s = dist_s[[x]],      # distributions names for random  vital rates
                               params.s = params_s[[x]],  # distributions parameters fore random vital rate
                               L.mean = L_mean[[x]],      # Initial (mean) length at age
                               data = c(lh_data[[x]], l_min$equil[[x]], s0 = NA), # Life history data 
                               sensitivity = FALSE)
    
    pmx_init <- pmx_rand$pmx
    
    rm(pmx_rand)
    gc()
    
    # Solve for syoy with lambda = max
    max <- optimize(M_optim_f, c(0, 1e-1), tol = 1e-16, target.lambda = lambda.max[[x]],
                    numReps = numReps,         # number of simulation replicates
                    years = years,             # number of years to simulate
                    pmxs = pmx_init,           # random projection matrcies
                    data =  c(lh_data[[x]], l_min$equil[[x]])) # Life history data 
    
    # Generate random variables
    pmx_rand <- pmx_build_runs(mx = A.LS,     # Projection matrix expression
                               numReps = numReps,         # number of simulation replicates
                               years = years,             # number of years to simulate
                               lambda = lambda.mean[[x]],  # pop growth rates
                               cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                               names.Li = names_Li[[x]],  # names for growth increments
                               dist.Li = dist_Li[[x]],    # distributions names for random length increments
                               params.Li = params_Li[[x]],# distributions parameters fore random length increments
                               cor.s = cor_s[[x]],        # Correlation struction for vital rates
                               names.s = names_s[[x]],    # names for vital rates
                               dist.s = dist_s[[x]],      # distributions names for random  vital rates
                               params.s = params_s[[x]],  # distributions parameters fore random vital rate
                               L.mean = L_mean[[x]],      # Initial (mean) length at age
                               data = c(lh_data[[x]], l_min$equil[[x]], s0 = NA), # Life history data 
                               sensitivity = FALSE)
    
    pmx_init <- pmx_rand$pmx
    
    rm(pmx_rand)
    gc()
    
    # Solve for syoy with lambda = mean
    mean <- optimize(M_optim_f, c(0, 1e-1), tol = 1e-16, target.lambda = lambda.mean[[x]],
                     numReps = numReps,         # number of simulation replicates
                     years = years,             # number of years to simulate
                     pmxs = pmx_init,           # random projection matrcies
                     data = c(lh_data[[x]], l_min$equil[[x]])) # Life history data 
    
    list(mean = mean$minimum, max = max$minimum)
  })
  
  stopCluster(cl)
  
  # Assign syoy to lh_data
  lh_data <- lapply(1:6, function(x) { 
    lh_data[[x]]$syoy.mean = M_optim[[x]]$mean # lambda = mean
    lh_data[[x]]$syoy.max = M_optim[[x]]$max   # lmabda = max
    lh_data[[x]]
  }) 
  
  rm(years, numReps, M_optim, M_optim_f)
}

# Set YOY survival
Ms.mean <- c(0.0005621, 0.0001758, 0.0000380, 0.0000973, 0.0000804, 0.0002073)
lh_data <- lapply(1:6, function(x) { lh_data[[x]]$syoy.mean = Ms.mean[x]; lh_data[[x]]})    

Ms.max <- c(0.0092573, 0.0018055, 0.0008689, 0.0054470, 0.0006173, 0.0022097)
lh_data <- lapply(1:6, function(x) { lh_data[[x]]$syoy.max = Ms.max[x]; lh_data[[x]]})    

# Mean annual Mortality
params_s.mean <- lapply(1:6, function(x) { # Survival
  
  # Mean annual Survival 
  s_means <- surv_est(0:(lh_data[[x]]$Tmax - 1), c(lh_data[[x]], s0 = lh_data[[x]]$syoy.mean))
  
  # Mean annual Instantanesou Mortality
  M_means <- {
    Ms <- -log(s_means)
    names(Ms) <- paste0("M", 0:(length(s_means) - 1))
    Ms
  }
  
  lapply(M_means, function(i) beta_stretch_val(i, sd = (i * lh_data[[x]]$cv.M)))
})

params_s.max <- lapply(1:6, function(x) { # Survival
  
  # Mean annual Survival 
  s_means <- surv_est(0:(lh_data[[x]]$Tmax - 1), c(lh_data[[x]], s0 = lh_data[[x]]$syoy.max))
  
  # Mean annual Instantanesou Mortality
  M_means <- {
    Ms <- -log(s_means)
    names(Ms) <- paste0("M", 0:(length(s_means) - 1))
    Ms
  }
  
  lapply(M_means, function(i) beta_stretch_val(i, sd = (i * lh_data[[x]]$cv.M)))
})

rm(params_s)

pmx.mean <- lapply(1:6, function(x) pmx_create(c(lh_data[[x]], s0 = lh_data[[x]]$syoy.mean, l_min$equil[[x]])))
pmx.max <- lapply(1:6, function(x) pmx_create(c(lh_data[[x]], s0 = lh_data[[x]]$syoy.max, l_min$equil[[x]])))

#-------------------------------------------------------------------------------
# minimum pop growth adjustments
#-------------------------------------------------------------------------------
# Optimize for reduction in vital rates to produce matrices with geometric mean
# lambdas of lambda min;
# Separtes reductions for yoy, young juv, olf juv, young adt, and old adt survival
# and fecundity 

if(FALSE){
  # Function to optimize
  M_optim_f = function (x, par, target.lambda, mx, numReps, years, ss, Fs, Gs, data){
    
    if(par == "M0") {        # Young juvenile mortality increase
      ss <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          ss[[n]][[y]]$syoy = ss[[n]][[y]]$syoy * x 
          ss[[n]][[y]]
          })
      })
      l.min <- list(M0 = x, # YOY mortality increase
                    Mj = 1, # Young juvenile mortality increase
                    Ma = 1, # Young adult mortality increase
                    Mf = 1) # Fecundity reduction
    } else if(par == "Mj") { # Young juvenile mortality increase
      ss <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          ss[[n]][[y]]$sj1 = ss[[n]][[y]]$sj1 * x
          ss[[n]][[y]]$sj2 = ss[[n]][[y]]$sj2 * x
          ss[[n]][[y]]
          })
      })
      l.min <- list(M0 = 1, # YOY mortality increase
                    Mj = x, # Young juvenile mortality increase
                    Ma = 1, # Young adult mortality increase
                    Mf = 1) # Fecundity reduction
    } else if(par == "Ma") { # Young adult mortality increase
      ss <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          ss[[n]][[y]]$sa1 = ss[[n]][[y]]$sa1 * x
          ss[[n]][[y]]$sa2 = ss[[n]][[y]]$sa2 * x
          ss[[n]][[y]]
        })
      })
      l.min <- list(M0 = 1, # YOY mortality increase
                    Mj = 1, # Young juvenile mortality increase
                    Ma = x, # Young adult mortality increase
                    Mf = 1) # Fecundity reduction
    } else if(par == "Mf") {  # Fecundity reduction
      Fs <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          Fs[[n]][[y]]$Fj = Fs[[n]][[y]]$Fj * x 
          Fs[[n]][[y]]$F1 = Fs[[n]][[y]]$F1 * x 
          Fs[[n]][[y]]$F2 = Fs[[n]][[y]]$F2 * x 
          Fs[[n]][[y]]
          })
      })
      l.min <- list(M0 = 1, # YOY mortality increase
                    Mj = 1, # Young juvenile mortality increase
                    Ma = 1, # Young adult mortality increase
                    Mf = x) # Fecundity reduction
    }
    
    # mean projection matrix
    pmx_mean <- pmx_create(c(data, s0 = data$syoy.1, l.min))
    
    # initial population size
    N0 <- initPop(pmx_mean, 100, "adult", p.rep = c(0, 0, 0, 1, 1))
    
    # stochastic projection matrices
    pmx_list <- lapply(1:numReps, function(x){
      lapply(1:years, function(y) pmx_eval(mx, 
                                           vars = c(data, ss[[x]][[y]], Fs[[x]][[y]], Gs[[x]][[y]])))
    })
    
    # est. geom_mean of lambda of pmx.list
    lambdas <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$longterm) 
    #lambdas <- sapply(pmx_list, lambda_mean)
    mean_lambda <- mean(lambdas)
    
    # est distance from trarget lambda
    (mean_lambda - target.lambda)^2
  }
  
  
  no_cores <- detectCores() - 1
  
  cl <- makeCluster(no_cores) 
  clusterExport(cl, ls())
  clusterEvalQ(cl, library(popbio))
  
  M_optim <- parLapply(cl, 1:6, function(x){
    
    years = 100
    numReps = 500
    
    # Generate random variables
    pmx_rand <- pmx_build_runs(mx = A.LS,                # Projection matrix expression
                               numReps = numReps,         # number of simulation replicates
                               years = years,             # number of years to simulate
                               lambda = lambda.min[[x]],  # pop growth rate
                               cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                               names.Li = names_Li[[x]],  # names for growth increments
                               dist.Li = dist_Li[[x]],    # distributions names for random length increments
                               params.Li = params_Li[[x]],# distributions parameters fore random length increments
                               cor.s = cor_s[[x]],        # Correlation struction for vital rates
                               names.s = names_s[[x]],    # names for vital rates
                               dist.s = dist_s[[x]],      # distributions names for random  vital rates
                               params.s = params_s.1[[x]],# distributions parameters fore random vital rate
                               L.mean = L_mean[[x]],      # Initial (mean) length at age
                               data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.1), # Life history data 
                               sensitivity = TRUE)
    
    st_init <- pmx_rand$ss
    Ft_init <- pmx_rand$Fs
    Gt_init <- pmx_rand$Gs
    
    rm(pmx_rand)
    
    # Solve for M0 with lambda = min
    cat("   M0", "\n")
    M0 <- optimize(M_optim_f, c(0, 1), tol = 1e-16, par = "M0", target.lambda = lambda.min[[x]],
                   mx = A.LS,                # Projection matrix expression
                   numReps = numReps,         # number of simulation replicates
                   years = years,             # number of years to simulate
                   ss = st_init,              # random survival rates
                   Fs = Ft_init,              # Random fertility
                   Gs = Gt_init,              # random stage transition rates
                   data = lh_data[[x]])       # Life history data 
    
    # Solve for Mj with lambda = min
    cat("   Mj", "\n")
    Mj <- optimize(M_optim_f, c(0, 1), tol = 1e-16, par = "Mj", target.lambda = lambda.min[[x]],
                   mx = A.LS,                # Projection matrix expression
                   numReps = numReps,         # number of simulation replicates
                   years = years,             # number of years to simulate
                   ss = st_init,              # random survival rates
                   Fs = Ft_init,              # Random fertility
                   Gs = Gt_init,              # random stage transition rates
                   data = lh_data[[x]])       # Life history data 
    
    # Solve for Ma with lambda = min
    cat("   Ma", "\n")
    Ma <- optimize(M_optim_f, c(0, 1), tol = 1e-16, par = "Ma", target.lambda = lambda.min[[x]],
                   mx = A.LS,                # Projection matrix expression
                   numReps = numReps,         # number of simulation replicates
                   years = years,             # number of years to simulate
                   ss = st_init,              # random survival rates
                   Fs = Ft_init,              # Random fertility
                   Gs = Gt_init,              # random stage transition rates
                   data = lh_data[[x]])       # Life history data 
    
    # Solve for Mf with lambda = min
    cat("   Mf", "\n")
    Mf <- optimize(M_optim_f, c(0, 1), tol = 1e-16, par = "Mf", target.lambda = lambda.min[[x]],
                   mx = A.LS,                # Projection matrix expression
                   numReps = numReps,         # number of simulation replicates
                   years = years,             # number of years to simulate
                   ss = st_init,              # random survival rates
                   Fs = Ft_init,              # Random fertility
                   Gs = Gt_init,              # random stage transition rates
                   data = lh_data[[x]])       # Life history data 
    
    list(M0 = M0$minimum, Mj = Mj$minimum,  Ma = Ma$minimum, Mf = Mf$minimum)
  })
  
  stopCluster(cl)
  
}

# Decresing vital rates giving min lambda
M0s <- c(0.1697, 0.4190, 0.3173, 0.8834, 0.1327, 0.1456) # YOY survival
Mjs <- c(0.8298, 0.9251, 0.9061, 0.9910, 0.8092, 0.8165) # Juvenile survival
Mas <- c(0.7294, 0.9125, 0.9353, 0.9938, 0.7562, 0.7397) # Adult survival
Mfs <- c(0.1696, 0.4189, 0.3172, 0.8833, 0.1327, 0.1455) # Fertility

# Params used to reduce lambda to lambda min
l_min$M0 <- lapply(M0s, function(x) {
  list(M0 = x, # YOY mortality increase
       Mj = 1, # Young juvenile mortality increase
       Ma = 1, # Young adult mortality increase
       Mf = 1) # Fecundity reduction
})

l_min$Mj <- lapply(Mjs, function(x) {
  list(M0 = 1, # YOY mortality increase
       Mj = x, # Young juvenile mortality increase
       Ma = 1, # Young adult mortality increase
       Mf = 1) # Fecundity reduction
})

l_min$Ma <- lapply(Mas, function(x) {
  list(M0 = 1, # YOY mortality increase
       Mj = 1, # Young juvenile mortality increase
       Ma = x, # Young adult mortality increase
       Mf = 1) # Fecundity reduction
})

l_min$Mf <- lapply(Mfs, function(x) {
  list(M0 = 1, # YOY mortality increase
       Mj = 1, # Young juvenile mortality increase
       Ma = 1, # Young adult mortality increase
       Mf = x) # Fecundity reduction
})

# Clean u
rm(M0s, Mjs, Mas, Mfs)

pmx.min <- lapply(2:5, function(y) lapply(1:6, function(x) {
  pmx_create(c(lh_data[[x]], s0 = lh_data[[x]]$syoy.1, l_min[[y]][[x]]))
}))
names(pmx.min) <- c("M0", "Mj", "Ma", "Mf")

#-------------------------------------------------------------------------------
# Stable Stage Distributions
#-------------------------------------------------------------------------------

if(FALSE){
# simulation criteria
numReps = 500
years = 100

# Est. random matrices for each cor str
no_cores <- detectCores() - 1

cl <- makeCluster(no_cores) 
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))

stable_stage.est <- parLapply(cl, 1:6, function(x){
  pmx_rand <- pmx_build_runs(mx = A.LS,       # Projection matrix expression
                             numReps = numReps,         # number of simulation replicates
                             years = years,             # number of years to simulate
                             lambda = lambda.1[[x]], # pop growth rate                                                ****
                             cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                             names.Li = names_Li[[x]],  # names for growth increments
                             dist.Li = dist_Li[[x]],    # distributions names for random length increments
                             params.Li = params_Li[[x]],# distributions parameters fore random length increments
                             cor.s = cor_s[[x]],        # Correlation struction for vital rates
                             names.s = names_s[[x]],    # names for vital rates
                             dist.s = dist_s[[x]],      # distributions names for random  vital rates
                             params.s = params_s.1[[x]],# distributions parameters fore random vital rate             ****
                             L.mean = L_mean[[x]],      # Initial (mean) length at age
                             data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.1), # Life history data  ****
                             sensitivity = FALSE)

  # extract results  
  pmx_list <- pmx_rand$pmx  # matrices
  
  # clean up
  rm(pmx_rand) 
  
  # Calc stable stage for each matric
  ss <- lapply(pmx_list, function(x) lapply(x, stable.stage))
  
  # MGeometric mean SS within replicate
  ss.mean <- lapply(ss, function(x) apply(matrix(unlist(x), ncol = 5, byrow = T),2,geom_mean))
  
  # Mean SS accross replicates
  stable_stage_est <- apply(matrix(unlist(ss.mean), ncol = 5, byrow = T),2,mean)
  stable_stage_est
})
stopCluster(cl)
}

stable_stage_est <- list(
  c(9.979382e-01, 6.376061e-04, 9.236090e-05, 7.963214e-05, 6.022005e-06),
  c(9.991370e-01, 1.873928e-04, 3.308316e-05, 4.081866e-05, 5.174556e-06),
  c(9.996922e-01, 3.218751e-05, 8.715567e-06, 2.325243e-05, 1.116471e-05),
  c(9.994175e-01, 7.957906e-05, 2.052591e-05, 4.757157e-05, 1.999389e-05),
  c(9.994948e-01, 8.655402e-05, 1.712313e-05, 2.601713e-05, 4.500543e-06),
  c(9.990459e-01, 2.208576e-04, 3.802090e-05, 4.447168e-05, 5.237405e-06)
)
names(stable_stage_est) <- DU

#-------------------------------------------------------------------------------
# Model Variability
#-------------------------------------------------------------------------------

# simulation criteria
numReps = 1000
years = 100
x = 2 # DU Location for random variability plotting

# Est. random matrices for each cor str
system.time({
  pmx_rand <- pmx_build_runs(mx = A.LS,       # Projection matrix expression
                             numReps = numReps,         # number of simulation replicates
                             years = years,             # number of years to simulate
                             lambda = lambda.1[[x]], # pop growth rate                                                ****
                             cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                             names.Li = names_Li[[x]],  # names for growth increments
                             dist.Li = dist_Li[[x]],    # distributions names for random length increments
                             params.Li = params_Li[[x]],# distributions parameters fore random length increments
                             cor.s = cor_s[[x]],        # Correlation struction for vital rates
                             names.s = names_s[[x]],    # names for vital rates
                             dist.s = dist_s[[x]],      # distributions names for random  vital rates
                             params.s = params_s.1[[x]],# distributions parameters fore random vital rate             ****
                             L.mean = L_mean[[x]],      # Initial (mean) length at age
                             data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.1), # Life history data  ****
                             sensitivity = FALSE)
})

# extract results  
pmx_list <- pmx_rand$pmx  # matrices
Lt_list <- pmx_rand$Lts   # vital rates
st_list <- pmx_rand$ss    # vital rates
Ft_list <- pmx_rand$Fs    # vital rates
Gt_list <- pmx_rand$Gs    # vital rates

rm(pmx_rand) 

# extract pop growth
lambda_year.d <- sapply(unlist(pmx_list, recursive = F), lambda) # individual years
lambda_rep.d <- sapply(pmx_list, lambda_mean)                    # mean deterministic lambda within replicate

lambda_year.e <- (lapply(1:length(pmx_list), function(i) lambda.real(N0 = initPop(pmx[[x]], 100, "adult", c(0, 0, 0, 1, 1)), pmx_list[[i]])$annual))
lambda_rep.e <- sapply(1:length(pmx_list), function(i) lambda.real(N0 = initPop(pmx[[x]], 100, "adult", c(0, 0, 0, 1, 1)), pmx_list[[i]])$longterm)

mean(lambda_rep.d) # mean lambda - overall
mean(lambda_rep.e)   # mean lambda - overall
lambda.min[[x]]

sd(lambda_rep.e)   # sd of lambda among replicates
sd(log(unlist(lambda_year.e)))   # sd of lambda among replicates
quantile(unlist(lambda_year.e), prob = c(0.005, 0.995)) # 99% CI for annual lambdas

# Annual pop growth rates
plot_data <- as.data.frame(list(
  lambda = c(unlist(lambda_year.e), unlist(lambda_rep.e)),
  DU = c(sort(rep(DU[x], years * numReps)),
          sort(rep(DU[x], numReps))),
  strip = c(rep("Annual", years * numReps), rep("Long Term", numReps))
))
plot_data <- plot_data[plot_data$lambda < 1.35 & plot_data$lambda > 0.85,]

p <- ggplot(plot_data, aes(x = lambda)) 
p <- p + geom_density(aes(fill = DU), alpha = 0.5)
p <- p + facet_wrap(~strip, scales = "free")
p <- p + labs(x="Lambda", y = "Density", fill = "DU")  
p <- p + guides(fill = FALSE)
p <- p + theme_me +
  theme(axis.text.y = element_blank())

#png("Figures/Lambda dist.png", width = 7, height = 3.5, unit = "in", res = 300)
p
#dev.off()

# Length 
plot_data <- as.data.frame(list(
  Length = as.vector(unlist(Lt_list)),
  Age = 0:100,
  Var = "Length"
))

p <- ggplot(plot_data[plot_data$Age %in% c(1, seq(10, 100, 10)), ], aes(x = Length)) 
p <- p + geom_density(aes(fill = as.factor(Age)), alpha = 0.75)
p <- p + facet_wrap(~Var)
p <- p + labs(x = "Length (mm)", y = "Density", fill = "Age")
p <- p + scale_x_continuous(breaks = seq(0, 2000, 400))

p <- p + guides(fill = guide_legend(ncol=2))
p <- p + theme_me + theme(legend.position = "right", axis.text.y = element_blank())
p 

# Fecundity 
plot_data <- as.data.frame(list(
  Fec = as.vector(unlist(Ft_list)),
  Age = c("Old Juvenile", "Young Adult", "Old Adult"),
  Var = "Fertility"
))
plot_data$Age <- factor(plot_data$Age, levels = c("Old Juvenile", "Young Adult", "Old Adult"))

q <- ggplot(plot_data, aes(x = Fec)) 
q <- q + geom_density(aes(fill = Age), alpha = 0.75)
q <- q + scale_x_continuous(trans = "log10", breaks = c(100, 1000, 10000),
                            minor_breaks = c(seq(200, 900, 100), seq(2000, 9000, 1000),
                                             11000))
q <- q + facet_wrap(~Var)
q <- q + labs(x = "Fertility", y = "Density", fill = "Stage")
q <- q + theme_me + theme(legend.position = "right", axis.text.y = element_blank())
q 

# Mortality
plot_data <- as.data.frame(list(
  M = -log(unlist(st_list)),
  Age = c("YOY", "Young Juvenile", "Old Juvenile", "Young Adult", "Old Adult"),
  Var = "Instantaneous Mortality"
))
plot_data$Age <- factor(plot_data$Age, levels = c("YOY", "Young Juvenile", "Old Juvenile", "Young Adult", "Old Adult"))


o <- ggplot(plot_data, aes(x = M)) 
o <- o + geom_density(aes(fill = as.factor(Age)), alpha = 0.75)
o <- o + scale_x_continuous(trans = "log10", breaks =  c(0.01, 0.05, 0.1, 0.5, 1, 5, 10), 
                            minor_breaks = c(seq(0.02, 0.09, 0.01),seq(0.2, 0.9, 0.1),
                                             seq(2, 9, 1), 20))
o <- o + facet_wrap(~Var)
o <- o + labs(x = "Instantaneous Mortality", y = "Density", fill = "Stage")
o <- o + theme_me + theme(legend.position = "right", axis.text.y = element_blank())
o  

# All together
#png("Figures/Variable dist.png", width = 8, height = 5, unit = "in", res = 300)
multiplot(p, o)
#dev.off()

#Clean up
rm(numReps, years, pmx_list, plot_data, p, o, q, x, Lt_list, Ft_list, st_list, Gt_list,
   lambda_rep.e, lambda_year.e, lambda_rep.d, lambda_year.d)

#-------------------------------------------------------------------------------
# Stochastic Elasticity
#-------------------------------------------------------------------------------
# Estimate elasticities values for vital rate. repeat for different pop growth 
# rates and for different correlation structure among random vital rates

# Simluation runs
numReps <- 500
years <- 100

# List of parameter values for vital rate distributions
# use to produce matrices with different mean pop growth rates
params_s <- list(
  M.min.M0 = params_s.1,  # Minimum pop growth rate 
  M.min.Mj = params_s.1,  # Minimum pop growth rate 
  M.min.Ma = params_s.1,  # Minimum pop growth rate 
  M.min.Mf = params_s.1,  # Minimum pop growth rate 
  M.1 = params_s.1,       # Equilibrium
  M.mean = params_s.mean, # mean pop growth rate 
  M.max = params_s.max    # maximu pop growth rate 
)


# Elasticity for different pop growth rates
# Apply elasticity_lambdas over list of vital rate parameters
system.time({
  no_cores <- detectCores()
  
  cl <- makeCluster(no_cores)
  clusterExport(cl, ls())
  
  elasts <- parLapply(cl, 1:length(params_s), function(par){
    
    lapply(1:6, function(x){
      
      if(par == 1) {
        l_mins <- l_min$M0
        lambda <- lambda.min[[x]]
        s.yoy <- lh_data[[x]]$syoy.1
      } else if (par == 2) {
        l_mins <- l_min$Mj
        lambda <- lambda.min[[x]]
        s.yoy <-lh_data[[x]]$syoy.1
      } else if (par == 3) {
        l_mins <- l_min$Ma
        lambda <- lambda.min[[x]]
        s.yoy <- lh_data[[x]]$syoy.1
      } else if (par == 4) {
        l_mins <- l_min$Mf
        lambda <- lambda.min[[x]]
        s.yoy <- lh_data[[x]]$syoy.1
      } else if (par == 5) {
        l_mins <- l_min$equil
        lambda <- lambda.1[[x]]
        s.yoy <- lh_data[[x]]$syoy.1
      } else if (par == 6) {
        l_mins <- l_min$equil
        lambda <- lambda.mean[[x]]
        s.yoy <- lh_data[[x]]$syoy.mean
      } else if (par == 7) {
        l_mins <- l_min$equil
        lambda <- lambda.max[[x]]
        s.yoy <- lh_data[[x]]$syoy.max
      }
      
      elasticity_lambdas(mx = A.LS,                      # Projection matrix expression
                         numReps = numReps,              # number of simulation replicates
                         years = years,                  # number of years to simulate
                         lambda = lambda,                # pop growth rate
                         cor.Li = cor_Li[[x]],           # Correlation struction for growth increments
                         names.Li = names_Li[[x]],       # names for growth increments
                         dist.Li = dist_Li[[x]],         # distributions names for random length increments
                         params.Li = params_Li[[x]],     # distributions parameters fore random length increments
                         cor.s = cor_s[[x]],             # Correlation struction for vital rates
                         names.s = names_s[[x]],         # names for vital rates
                         dist.s = dist_s[[x]],           # distributions names for random  vital rates
                         params.s = params_s[[par]][[x]],# distributions parameters fore random vital rate
                         L.mean = L_mean[[x]],           # Initial (mean) length at age
                         data = c(lh_data[[x]], l_mins[[x]], s0 = s.yoy))# Life history data
    })
  })
  stopCluster(cl)
})

# Save sensitivity analysis results
#save(elasts, file = "Results 01 - elasticity.Rdata")

# Clean up
rm(elasts, years, numReps)

#-------------------------------------------------------------------------------
# Extract and organize MVP simulatin results

load("Results 01 - elasticity.Rdata")

elasts <- lapply(elasts, function(x) lapply(x, function(y) lapply(y, function(z) {
  m <- matrix(unlist(z), nrow = years, byrow = T)
  colnames(m) = names(z[[1]])
  m
})))

# sums age specific VR to stage-specific  
elast_sum <- function(x){
  
  syoy <- x[ ,1]
  sj1 <- x[ ,2]
  sj2 <- x[ ,3]
  sa1 <- x[ ,4]
  sa2 <- x[ ,5]
  sj <- apply(x[, 2:3], 1, sum)
  sa <- apply(x[, 4:5], 1, sum)
  s <- apply(x[, 1:5], 1, sum)
  s1plus <- apply(x[, 2:5], 1, sum)

  f <- apply(x[, 6:8], 1, sum)
  
  e <- cbind(syoy, sj1, sj2, sa1, sa2, f, sj, sa, s, s1plus)
  e
}

# Sum elasticities to stage specific results
elasts_comp <- lapply(elasts, function(x) lapply(x, function(y) lapply(y, elast_sum)))

# calc mean elasticity values over years, within simluation
# also est. confidence intrvals

# stage-specific vital rates
elasts_comp_rep <- lapply(elasts_comp, function(x) lapply(x, elasticity_mean))

# estimate means and CI over replicates

# takes mean of elasticity values over replicates
elast_mean <- function(x, loc, vrNum, names) {
  
  e <- apply(matrix(unlist(x[[loc]]), ncol = vrNum, byrow = TRUE), 2, mean)
  names(e) <- names
  e
}

# names to assign to vectors
names_stage <- c("syoy", "sj1", "sj2", "sa1", "sa2", "F", "sj",  "sa", "s", "s1plus")       # stage-specific results

# stage-specific with independent VR
elasts_comp_mean <- lapply(elasts_comp_rep, function(x) lapply(x, elast_mean, 1, vrNum = 10, names = names_stage))
elasts_comp_LCI <- lapply(elasts_comp_rep, function(x) lapply(x, elast_mean, 2, vrNum = 10, names = names_stage))
elasts_comp_UCI <- lapply(elasts_comp_rep, function(x) lapply(x, elast_mean, 3, vrNum = 10, names = names_stage))

# Clean up
rm(names_stage,  elasts_comp, elasts_comp_rep, elasts, elast_mean, elast_sum)

#-------------------------------------------------------------------------------
# Plot elsticity results

# Stage based elasticities
names_stage <- c("syoy", "sj1", "sj2", "sa1", "sa2", "F", "sj",  "sa", "s", "s1plus")  
lambda <- c("min s0", "min sj", "min sa", "min f","1", "mean", "max")
plotData <- as.data.frame(list(  
  DU = rep(sort(rep(DU, length(names_stage))), length(lambda)),
  VR = rep(rep(names_stage, length(DU)), length(lambda)),
  Lambda = lambda[sort(rep(rep(1:length(lambda), length(names_stage)), length(DU)))],
  Elast = unlist(elasts_comp_mean),
  LCI = unlist(elasts_comp_LCI), 
  UCI = unlist(elasts_comp_UCI)
))
plotData <- plotData[!(plotData$VR %in% c("sj", "sa", "s","s1plus")), ]

plotData$Lambda <- factor(plotData$Lambda, levels = lambda,
                          labels = lambda)
plotData$VR <- factor(plotData$VR, levels = c("syoy", "sj1", "sj2", "sa1", "sa2", "F"),
                      labels = c("syoy", "sj1", "sj2", "sa1", "sa2", "F"))
rm(names_stage, lambda)

labels <- c(expression(paste(italic(sigma[YOY]))),  # Survival - YOY
            expression(paste(italic(sigma[paste(j,"1")]))),  # Survival - juvenile
            expression(paste(italic(sigma[paste(j,"2")]))),
            expression(paste(italic(sigma[paste(a,"1")]))),
            expression(paste(italic(sigma[paste(a,"2")]))),
            expression(paste(italic("F")))) # Fecundity)  # Survival - adult

guides <- c(expression(paste(lambda[min]," - ", italic(sigma[YOY]))), # min - YOY
            expression(paste(lambda[min]," - ", italic(sigma[j]))),   # min - juv
            expression(paste(lambda[min]," - ", italic(sigma[a]))),   # min - adt
            expression(paste(lambda[min]," - ",italic(F))),           # min - F
            expression(paste(lambda[1])),                             # 1
            expression(paste(lambda[mean])),                          # mean
            expression(paste(lambda[max])))                           # max

p <- ggplot(plotData, aes(x = VR, y = Elast, fill = Lambda)) 
p <- p + geom_bar(position = position_dodge(), stat="identity")
p <- p + geom_errorbar(aes(ymin = LCI, ymax = UCI),
                       width = .25,                    # Width of the error bars
                       position = position_dodge(0.9))
p <- p + facet_wrap(~DU, ncol = 2)
p <- p + scale_x_discrete(labels = labels)
p <- p + labs(x = "Vital Rates", y = "Elasticities")  
p <- p + scale_fill_discrete(name = NULL, labels = guides, guide = guide_legend(nrow = 1))
p <- p + theme_me + theme(legend.position = "top", legend.text.align = 0)

#png("Figures/Sensitivity - Stage.png", width = 7, height = 8, unit = "in", res = 300)
p
#dev.off()

# clean up
rm(p, plotData, labels)

#-------------------------------------------------------------------------------
# CHRONIC ALLOWABLE HARM
#-------------------------------------------------------------------------------
# Estimate allowable harm and recovery effort using vital rate elasticity
# Allowable harm - the proportional decrease in vital rate that would not cause
# lambda to decrease to < 1
# A_harm = 

# est. Allowable harm
A_harm.mean <- lapply(1:6, function(x) {
  allowable_harm(elasts_comp_UCI[[6]][[x]], Lambda = lambda.mean[[x]])
})

A_harm.max <- lapply(1:6, function(x) {
  allowable_harm(elasts_comp_UCI[[7]][[x]], Lambda = lambda.max[[x]])
})

# est. recovery action
Rec <- lapply(1:4, function(y) lapply(1:6, function(x) {
  allowable_harm(elasts_comp_mean[[y]][[x]], Lambda = lambda.min[[x]])
}))
  
# clean up
rm(elasts_comp_mean, elasts_comp_LCI, elasts_comp_UCI)

#-------------------------------------------------------------------------------
# Simlation of impacts of chronic harm to population growth rate
# Compare the distributions of lambdas values for populaitons that are not and are
# experiencing harm - start with mean pop growth rates (~1.14)

# simluation criteria
numReps <- 500
years <- 100

# Calc mean values
s_means <- lapply(lh_data, function(x) surv_est(0:(x$Tmax - 1), c(x, s0 = x$syoy.mean)))

# Amount of harm
h <- seq(-0.00, -0.99, -0.02)
nullMat <- matrix(0, ncol = 100, nrow = length(h))
h.y <- nullMat; h.y[, 1] <- h
h.j <- nullMat; h.j[, with(lh_data[[1]], 2:j2.Tmax)] <- h
h.a <- nullMat; h.a[, with(lh_data[[1]], (j2.Tmax+1):Tmax)] <- h
h.1plus <- cbind(rep(0, length(h)), matrix(rep(h, 99), ncol = 99))
h.all <- matrix(rep(h, 100), ncol = 100)

# Population vector
N0 <- lapply(1:6, function(x) initPop(pmx.mean[[x]], 100, "adult", c(0, 0, 0, 1, 1)))

# Simulate effects of harm
# compares un-impacted and impacted pop growth rates
# harm only applies to survival rate

# YOY harm

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
harm_lambda.y <- parLapply(cl, 1:6, function(y) {
  cat("Location", y, "Of 6", "\n")
  apply(h.y, 1, function(x) {
    A_harm_impact(s_means = s_means[[y]],    # mean unharmed survival rates
                  harm = x,                  # proportional change to survival
                  mx = A.LS,                 # Projection matrix expression
                  numReps = numReps,         # number of simulation replicates
                  years = years,             # number of years to simulate
                  N0 = N0[[y]],              # inital pop size vector
                  lambda = lambda.mean[[y]], # Pop growth rate
                  cor.Li = cor_Li[[y]],      # Correlation struction for growth increments
                  names.Li = names_Li[[y]],  # names for growth increments
                  dist.Li = dist_Li[[y]],    # distributions names for random length increments
                  params.Li = params_Li[[y]],# distributions parameters fore random length increments
                  cor.s = cor_s[[y]],        # Correlation struction for vital rates
                  names.s = names_s[[y]],    # names for vital rates
                  dist.s = dist_s[[y]],      # distributions names for random  vital rates
                  L.mean = L_mean[[y]],      # Initial (mean) length at age
                  data = c(lh_data[[y]], l_min$equil[[y]], s0 = lh_data[[y]]$syoy.mean)) # Life history data
  })
})
stopCluster(cl)

# Save sensitivity analysis results
save(harm_lambda.y, file = "Results 01b.1 - Harm impact - YOY.Rdata")

# clean up
rm(harm_lambda.y)

# Juv Harm 

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
harm_lambda.j <- parLapply(cl, 1:6, function(y) {
  cat("Location", y, "Of 6", "\n")
  apply(h.j, 1, function(x) {
    A_harm_impact(s_means[[y]],              # mean unharmed survival rates
                  harm = x,                  # proportional change to survival
                  mx = A.LS,                 # Projection matrix expression
                  numReps = numReps,         # number of simulation replicates
                  years = years,             # number of years to simulate
                  N0 = N0[[y]],              # inital pop size vector
                  lambda = lambda.mean[[y]], # Pop growth rate
                  cor.Li = cor_Li[[y]],      # Correlation struction for growth increments
                  names.Li = names_Li[[y]],  # names for growth increments
                  dist.Li = dist_Li[[y]],    # distributions names for random length increments
                  params.Li = params_Li[[y]],# distributions parameters fore random length increments
                  cor.s = cor_s[[y]],        # Correlation struction for vital rates
                  names.s = names_s[[y]],    # names for vital rates
                  dist.s = dist_s[[y]],      # distributions names for random  vital rates
                  L.mean = L_mean[[y]],      # Initial (mean) length at age
                  data = c(lh_data[[y]], l_min$equil[[y]], s0 = lh_data[[y]]$syoy.mean)) # Life history data
  })
})
stopCluster(cl)

# Save sensitivity analysis results
save(harm_lambda.j, file = "Results 01b.2 - Harm impact - Juv.Rdata")

# clean up
rm(harm_lambda.j)

# Adt Harm 

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
harm_lambda.a <- parLapply(cl, 1:6, function(y) {
  cat("Location", y, "Of 6", "\n")
  apply(h.a, 1, function(x) {
    A_harm_impact(s_means[[y]],              # mean unharmed survival rates
                  harm = x,                  # proportional change to survival
                  mx = A.LS,                 # Projection matrix expression
                  numReps = numReps,         # number of simulation replicates
                  years = years,             # number of years to simulate
                  N0 = N0[[y]],              # inital pop size vector
                  lambda = lambda.mean[[y]], # Pop growth rate
                  cor.Li = cor_Li[[y]],      # Correlation struction for growth increments
                  names.Li = names_Li[[y]],  # names for growth increments
                  dist.Li = dist_Li[[y]],    # distributions names for random length increments
                  params.Li = params_Li[[y]],# distributions parameters fore random length increments
                  cor.s = cor_s[[y]],        # Correlation struction for vital rates
                  names.s = names_s[[y]],    # names for vital rates
                  dist.s = dist_s[[y]],      # distributions names for random  vital rates
                  L.mean = L_mean[[y]],      # Initial (mean) length at age
                  data = c(lh_data[[y]], l_min$equil[[y]], s0 = lh_data[[y]]$syoy.mean)) # Life history data
  })
})
stopCluster(cl)

# Save sensitivity analysis results
save(harm_lambda.a, file = "Results 01b.3 - Harm impact - Adt.Rdata")

# clean up
rm(harm_lambda.a)

# 1+ Harm 

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
harm_lambda.1plus <- parLapply(cl, 1:6, function(y) {
  cat("Location", y, "Of 6", "\n")
  apply(h.1plus, 1, function(x) {
    A_harm_impact(s_means[[y]],              # mean unharmed survival rates
                  harm = x,                  # proportional change to survival
                  mx = A.LS,                 # Projection matrix expression
                  numReps = numReps,         # number of simulation replicates
                  years = years,             # number of years to simulate
                  N0 = N0[[y]],              # inital pop size vector
                  lambda = lambda.mean[[y]], # Pop growth rate
                  cor.Li = cor_Li[[y]],      # Correlation struction for growth increments
                  names.Li = names_Li[[y]],  # names for growth increments
                  dist.Li = dist_Li[[y]],    # distributions names for random length increments
                  params.Li = params_Li[[y]],# distributions parameters fore random length increments
                  cor.s = cor_s[[y]],        # Correlation struction for vital rates
                  names.s = names_s[[y]],    # names for vital rates
                  dist.s = dist_s[[y]],      # distributions names for random  vital rates
                  L.mean = L_mean[[y]],      # Initial (mean) length at age
                  data = c(lh_data[[y]], l_min$equil[[y]], s0 = lh_data[[y]]$syoy.mean)) # Life history data
  })
})
stopCluster(cl)

# Save sensitivity analysis results
save(harm_lambda.1plus, file = "Results 01b.4 - Harm impact - 1+.Rdata")

# clean up
rm(harm_lambda.1plus)

# All Harm 

no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
harm_lambda.all <- parLapply(cl, 1:6, function(y) {
  cat("Location", y, "Of 6", "\n")
  apply(h.all, 1, function(x) {
    A_harm_impact(s_means[[y]],              # mean unharmed survival rates
                  harm = x,                  # proportional change to survival
                  mx = A.LS,                 # Projection matrix expression
                  numReps = numReps,         # number of simulation replicates
                  years = years,             # number of years to simulate
                  N0 = N0[[y]],              # inital pop size vector
                  lambda = lambda.mean[[y]], # Pop growth rate
                  cor.Li = cor_Li[[y]],      # Correlation struction for growth increments
                  names.Li = names_Li[[y]],  # names for growth increments
                  dist.Li = dist_Li[[y]],    # distributions names for random length increments
                  params.Li = params_Li[[y]],# distributions parameters fore random length increments
                  cor.s = cor_s[[y]],        # Correlation struction for vital rates
                  names.s = names_s[[y]],    # names for vital rates
                  dist.s = dist_s[[y]],      # distributions names for random  vital rates
                  L.mean = L_mean[[y]],      # Initial (mean) length at age
                  data = c(lh_data[[y]], l_min$equil[[y]], s0 = lh_data[[y]]$syoy.mean)) # Life history data
  })
})
stopCluster(cl)

# Save sensitivity analysis results
save(harm_lambda.all, file = "Results 01b.5 - Harm impact - All.Rdata")

# clean up
rm(harm_lambda.all)

#-------------------------------------------------------------------------------
# Allowable harm simulation results

load("Results 01b.1 - Harm impact - YOY.Rdata")
load("Results 01b.2 - Harm impact - Juv.Rdata")
load("Results 01b.3 - Harm impact - Adt.Rdata")
load("Results 01b.4 - Harm impact - 1+.Rdata")
load("Results 01b.5 - Harm impact - All.Rdata")

# Probability of population delcine

# Annual pop growth
prob_decline_year <- lapply(1:6, function(x) {
  data <- as.data.frame(list(
    H = rep(-h, 5), 
    P = c(sapply(harm_lambda.y[[x]], function(y) length(which(y$lambda_year < 1)) / (numReps * years)),
          sapply(harm_lambda.j[[x]], function(y) length(which(y$lambda_year < 1)) / (numReps * years)),
          sapply(harm_lambda.a[[x]], function(y) length(which(y$lambda_year < 1)) / (numReps * years)),
          sapply(harm_lambda.1plus[[x]], function(y) length(which(y$lambda_year < 1)) / (numReps * years)),
          sapply(harm_lambda.all[[x]], function(y) length(which(y$lambda_year < 1)) / (numReps * years))),
    Stage = c(rep("YOY", length(h)), 
              rep("Juvenile", length(h)), 
              rep("Adult", length(h)),
              rep("Age-1+", length(h)),
              rep("All Ages", length(h))),
    Type = "Annual",
    DU = DU[[x]]
  ))
  data$Stage <- factor(data$Stage, levels = unique(data$Stage))
  data
})
prob_decline_year <- do.call(rbind, prob_decline_year)

# mean pop growth
prob_decline_mean <- lapply(1:6, function(x) {
  data <- as.data.frame(list(
    H = rep( -h, 5), 
    P = c(sapply(harm_lambda.y[[x]], function(y) length(which(y$lambda_mean < 1)) / (numReps)),
          sapply(harm_lambda.j[[x]], function(y) length(which(y$lambda_mean < 1)) / (numReps)),
          sapply(harm_lambda.a[[x]], function(y) length(which(y$lambda_mean < 1)) / (numReps)),
          sapply(harm_lambda.1plus[[x]], function(y) length(which(y$lambda_mean < 1)) / (numReps)),
          sapply(harm_lambda.all[[x]], function(y) length(which(y$lambda_mean < 1)) / (numReps))),
    Stage = c(rep("YOY", length(h)), 
              rep("Juvenile", length(h)), 
              rep("Adult", length(h)),
              rep("Age-1+", length(h)),
              rep("All Ages", length(h))),
    Type = "Long Term",
    DU = DU[[x]]
  ))
  data$Stage <- factor(data$Stage, levels = unique(data$Stage))
  data
})
prob_decline_mean <- do.call(rbind, prob_decline_mean)

# Plot
b <- seq(0, 1, 0.25)
l <- c(expression(frac(0, 100)), expression(frac(25, 100)), expression(frac(50, 100)), 
       expression(frac(75, 100)), expression(frac(100, 100)))
p <- ggplot(prob_decline_year, aes(x = H, y = P)) 
p <- p + geom_path(aes(colour = Stage), size = 1.5)
p <- p + facet_wrap(~ DU, ncol = 3)
p <- p + labs(x = "Harm", y = "Probability of Population decline") 
p <- p + scale_x_continuous(breaks = b, labels = l)
p <- p + theme_me #+ 
  #theme(legend.position = "top")

#png("Figures/Allowable harm - annual decline prob.png", width = 8, height = 5, unit = "in", res = 300)
p
#dev.off()

b <- seq(0, 1, 0.25)
l <- c(expression(frac(0, 100)), expression(frac(25, 100)), expression(frac(50, 100)), 
       expression(frac(75, 100)), expression(frac(100, 100)))
p <- ggplot(prob_decline_mean, aes(x = H, y = P)) 
p <- p + geom_path(aes(colour = Stage), size = 1.5)
p <- p + facet_wrap(~DU, ncol = 3)
p <- p + labs(x = "Harm", y = "Probability of Population decline") 
p <- p + scale_x_continuous(breaks = b, labels = l)
p <- p + theme_me #+ 
#theme(legend.position = "top")

#png("Figures/Allowable harm - 100 year decline prob.png", width = 8, height = 5, unit = "in", res = 300)
p
#dev.off()


# Lambda Distributions
# Annual pop growth rates
lambda_year <- as.data.frame(list(
  lambda = c(harm_lambda.all[[4]][[1]]$lambda_year,
             harm_lambda.all[[4]][[2]]$lambda_year),
  Harm = c(rep("Unharmed", years * numReps),
           rep("Harmed", years * numReps)),
  Type = "Annual"
))
lambda_year <- lambda_year[lambda_year$lambda <= 1.5, ]

# Mean pop growth over years
lambda_long <- as.data.frame(list(
  lambda = c(harm_lambda.all[[4]][[1]]$lambda_mean,
             harm_lambda.all[[4]][[2]]$lambda_mean),
  Harm = c(rep("Unharmed", numReps),
           rep("Harmed", numReps)),
  Type = "Long Term"
))

# plot
p <- ggplot(rbind(lambda_year, lambda_long), aes(x = lambda)) 
p <- p + geom_density(aes(fill = Harm), alpha = 0.5)
p <- p + facet_wrap(~ Type, ncol = 2, scales = "free_y")
#p <- p + scale_x_continuous(limits = c(0.5, 2.5))
p <- p + labs(x = "Lambda", y = "Density")           
p <- p + theme_me + 
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        axis.text.y = element_blank())

#png("Figures/Allowable harm lambda.png", width = 7, height = 3.5, unit = "in", res = 300)
p
#dev.off()

# clean up
rm(years, numReps, h, s_means, lambda_year, lambda_long, p, harm_impact, h.y,
   h.1plus, h.all, unharmed_lambda, harm_lambda.y,harm_lambda.j, harm_lambda.a,
   harm_lambda.1plus, harm_lambda.all, p.uh, prob_decline_mean, prob_decline_year)

#-------------------------------------------------------------------------------
# TRANSIENT ALLOWABLE HARM
#-------------------------------------------------------------------------------
# Estimate the impact of a one time removal of fish on populaiton growth rate
# over 3 generation or 10 years (whatever is shorter)

# Simluation criteria
numReps <- 5000
years <- 10

# Amount of harm
h <- seq(0.00, 0.99, 0.01) # inividual replicates
m <- rep(0, 5)
h_y <- lapply(as.list(h), function(x) {h.m <- m; h.m[1] <- x; h.m})       # harm to YOY
h_j <- lapply(as.list(h), function(x) {h.m <- m; h.m[2:3] <- x; h.m})     # harm to juveniles
h_a <- lapply(as.list(h), function(x) {h.m <- m; h.m[4:5] <- x; h.m})     # harm to adults
h_1plus <- lapply(as.list(h), function(x) {h.m <- m; h.m[2:5] <- x; h.m}) # harm to age-1+
h_all <- lapply(as.list(h), function(x) {h.m <- m; h.m[1:5] <- x; h.m})   # harm to all

# YOY 
no_cores <- detectCores() -1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
T_harm_yoy <- parLapply(cl, 1:6, function(x) {
  cat("Location", x, "Of 6", "\n")
  T_harm_impact(harm = h_y,                       # proportional change to survival
                mx = A.LS,                        # Projection matrix expression
                numReps = numReps,                # number of simulation replicates
                years = years,                    # number of years to simulate
                lambda = lambda.mean[[x]],        # pop grwoth rate
                cor.Li = cor_Li[[x]],             # Correlation struction for growth increments
                names.Li = names_Li[[x]],         # names for growth increments
                dist.Li = dist_Li[[x]],           # distributions names for random length increments
                params.Li = params_Li[[x]],       # distributions parameters fore random length increments
                cor.s = cor_s[[x]],               # Correlation struction for vital rates
                names.s = names_s[[x]],           # names for vital rates
                dist.s = dist_s[[x]],             # distributions names for random  vital rates
                params.s = params_s.mean[[x]],    # distributions parameters fore random vital rate
                L.mean = L_mean[[x]],             # Initial (mean) length at age
                data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.mean))# Life history data 
})

# Juveile
T_harm_juv <- parLapply(cl, 1:6, function(x) {
  cat("Location", x, "Of 6", "\n")
  T_harm_impact(harm = h_j,                       # proportional change to survival
                mx = A.LS,                        # Projection matrix expression
                numReps = numReps,                # number of simulation replicates
                years = years,                    # number of years to simulate
                lambda = lambda.mean[[x]],        # pop grwoth rate
                cor.Li = cor_Li[[x]],             # Correlation struction for growth increments
                names.Li = names_Li[[x]],         # names for growth increments
                dist.Li = dist_Li[[x]],           # distributions names for random length increments
                params.Li = params_Li[[x]],       # distributions parameters fore random length increments
                cor.s = cor_s[[x]],               # Correlation struction for vital rates
                names.s = names_s[[x]],           # names for vital rates
                dist.s = dist_s[[x]],             # distributions names for random  vital rates
                params.s = params_s.mean[[x]],    # distributions parameters fore random vital rate
                L.mean = L_mean[[x]],             # Initial (mean) length at age
                data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.mean))# Life history data 
})

# Adult
T_harm_adt <- parLapply(cl, 1:6, function(x) {
  cat("Location", x, "Of 6", "\n")
  T_harm_impact(harm = h_a,                       # proportional change to survival
                mx = A.LS,                        # Projection matrix expression
                numReps = numReps,                # number of simulation replicates
                years = years,                    # number of years to simulate
                lambda = lambda.mean[[x]],        # pop grwoth rate
                cor.Li = cor_Li[[x]],             # Correlation struction for growth increments
                names.Li = names_Li[[x]],         # names for growth increments
                dist.Li = dist_Li[[x]],           # distributions names for random length increments
                params.Li = params_Li[[x]],       # distributions parameters fore random length increments
                cor.s = cor_s[[x]],               # Correlation struction for vital rates
                names.s = names_s[[x]],           # names for vital rates
                dist.s = dist_s[[x]],             # distributions names for random  vital rates
                params.s = params_s.mean[[x]],    # distributions parameters fore random vital rate
                L.mean = L_mean[[x]],             # Initial (mean) length at age
                data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.mean))# Life history data 
})

# 1plus 
T_harm_1plus <- parLapply(cl, 1:6, function(x) {
  cat("Location", x, "Of 6", "\n")
  T_harm_impact(harm = h_1plus,                   # proportional change to survival
                mx = A.LS,                        # Projection matrix expression
                numReps = numReps,                # number of simulation replicates
                years = years,                    # number of years to simulate
                lambda = lambda.mean[[x]],        # pop grwoth rate
                cor.Li = cor_Li[[x]],             # Correlation struction for growth increments
                names.Li = names_Li[[x]],         # names for growth increments
                dist.Li = dist_Li[[x]],           # distributions names for random length increments
                params.Li = params_Li[[x]],       # distributions parameters fore random length increments
                cor.s = cor_s[[x]],               # Correlation struction for vital rates
                names.s = names_s[[x]],           # names for vital rates
                dist.s = dist_s[[x]],             # distributions names for random  vital rates
                params.s = params_s.mean[[x]],    # distributions parameters fore random vital rate
                L.mean = L_mean[[x]],             # Initial (mean) length at age
                data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.mean))# Life history data 
})

# All 
T_harm_all <- parLapply(cl, 1:6, function(x) {
  cat("Location", x, "Of 6", "\n")
  T_harm_impact(harm = h_all,                     # proportional change to survival
                mx = A.LS,                        # Projection matrix expression
                numReps = numReps,                # number of simulation replicates
                years = years,                    # number of years to simulate
                lambda = lambda.mean[[x]],        # pop grwoth rate
                cor.Li = cor_Li[[x]],             # Correlation struction for growth increments
                names.Li = names_Li[[x]],         # names for growth increments
                dist.Li = dist_Li[[x]],           # distributions names for random length increments
                params.Li = params_Li[[x]],       # distributions parameters fore random length increments
                cor.s = cor_s[[x]],               # Correlation struction for vital rates
                names.s = names_s[[x]],           # names for vital rates
                dist.s = dist_s[[x]],             # distributions names for random  vital rates
                params.s = params_s.mean[[x]],    # distributions parameters fore random vital rate
                L.mean = L_mean[[x]],             # Initial (mean) length at age
                data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.mean))# Life history data 
})

stopCluster(cl)

save(T_harm_yoy, T_harm_juv, T_harm_1plus, T_harm_adt, T_harm_all, file = "Results 02 - transient Harm.Rdata")

# clean up
rm(T_harm_yoy, T_harm_1plus, T_harm_adt, T_harm_all, params_vr, h_yoy,
   h_juv, h_adt, h_all, h)

#-------------------------------------------------------------------------------
# Extract and plot transient harm results

load(file = "Results 02 - transient Harm.Rdata")

# Amount of harm
h <- seq(0.00, 0.99, 0.01) 

T_harm_data.y <- lapply(1:6, function(x) {
  D <- lapply(1:length(h), function(y){
  
    as.data.frame(list(
      h = h[[y]],
      p = length(which(T_harm_yoy[[x]][[y]] < 1))/numReps,
      stage = "YOY",
      DU = DU[[x]]
    ))
  })
do.call(rbind, D)
})
T_harm_data.y <- do.call(rbind, T_harm_data.y)

T_harm_data.j <- lapply(1:6, function(x) {
  D <- lapply(1:length(h), function(y){
    
    as.data.frame(list(
      h = h[[y]],
      p = length(which(T_harm_juv[[x]][[y]] < 1))/numReps,
      stage = "Juvenile",
      DU = DU[[x]]
    ))
  })
  do.call(rbind, D)
})
T_harm_data.j <- do.call(rbind, T_harm_data.j)

T_harm_data.a <- lapply(1:6, function(x) {
  D <- lapply(1:length(h), function(y){
    
    as.data.frame(list(
      h = h[[y]],
      p = length(which(T_harm_adt[[x]][[y]] < 1))/numReps,
      stage = "Adult",
      DU = DU[[x]]
    ))
  })
  do.call(rbind, D)
})
T_harm_data.a <- do.call(rbind, T_harm_data.a)

T_harm_data.1plus <- lapply(1:6, function(x) {
  D <- lapply(1:length(h), function(y){
    
    as.data.frame(list(
      h = h[[y]],
      p = length(which(T_harm_1plus[[x]][[y]] < 1))/numReps,
      stage = "1+",
      DU = DU[[x]]
    ))
  })
  do.call(rbind, D)
})
T_harm_data.1plus <- do.call(rbind, T_harm_data.1plus)

T_harm_data.all <- lapply(1:6, function(x) {
  D <- lapply(1:length(h), function(y){
    
    as.data.frame(list(
      h = h[[y]],
      p = length(which(T_harm_all[[x]][[y]] < 1))/numReps,
      stage = "All Ages",
      DU = DU[[x]]
    ))
  })
  do.call(rbind, D)
})
T_harm_data.all <- do.call(rbind, T_harm_data.all)

# PLOT
b <- seq(0, 1, 0.25)
l <- c(expression(frac(0, 100)), expression(frac(25, 100)), expression(frac(50, 100)), 
       expression(frac(75, 100)), expression(frac(100, 100)))
p <- ggplot(rbind(T_harm_data.y, T_harm_data.j, T_harm_data.a, T_harm_data.1plus, T_harm_data.all), 
            aes(x = h, y = p)) 
p <- p + geom_path(aes(colour = stage), size = 1.5)
p <- p + facet_wrap(~ DU, ncol = 3)
p <- p + labs(x = "Harm", y = "Probability of Population decline", colour = "Stage")      
p <- p + scale_x_continuous(breaks = b, labels = l)
p <- p + theme_me 

png("Figures/Transient harm - decline prob.png", width = 8, height = 5, unit = "in", res = 300)
p
dev.off()

# clean up
rm(T_harm_yoy, T_harm_juv, T_harm_adt, T_harm_1plus, T_harm_all,h, p,
   T_harm_data.y, T_harm_data.j, T_harm_data.a, T_harm_data.1plus, T_harm_data.all)

#===============================================================================
# MVP & MAPV
#===============================================================================
# MVP ANALYSIS
# Generature stochastic matrices; run pop simulation with % fo catastrphe and 
# calc prob of extripation for diff adult pop size; fit relationship; calc MVP
# for defined prob of extriction

# Simulation runs
maxRuns <- 5      # number of independent runs
numReps <- 100    # number of replicated simulations
years <- 1500     # years to run simulation

# MVP parameters
p.cat <- c(0.15)  # prob of catastrophe
e.cat <- 0.5      # effect of catastrophe

# Mean adult pop sizes to loop through
adultPop <- c(seq(30, 90, 10),
              seq(100, 450, 50),
              seq(500, 900, 100), 
              seq(1000, 1750, 250),
              seq(2000, 4500, 500),
              seq(5000, 15000, 1000))


# maturity
p.rep <- c(rep(0, 3), rep(1, 2))

# Run simulations
# loop through inital adult pop sizes, run numreps sims of years maxRuns times, counts number
# of extripations (<= pop.fail adults). them calc the prob of exticntion for a given adult
# pop size

no_cores <- detectCores()  - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
MVP_run <- parLapply(cl, 1:6, function(x) {
  
  cat("\n")
  
  lapply(as.list(adultPop), function(y) {
      
    # indicator
    cat("    pop size:", y, "\n")
    cat("    ", date(), "\n")
    
    # inital population size
    N0 <- initPop(pmx[[x]], y, "adult", p.rep)
    
    rep_proj(mx = A.LS,                       # projection matrices as expression
             N0 = N0,                         # initial pop size
             maxRuns = maxRuns,               # number of runs of simlulation (new matrices est each run) 
             years = years,                   # years to run simulation
             numReps = numReps,               # Number of simulation per run
             lambda = lambda.1[[x]],          # Pop growth rate
             cor.Li = cor_Li[[x]],            # Correlation struction for growth increments
             names.Li = names_Li[[x]],        # names for growth increments
             dist.Li = dist_Li[[x]],          # distributions names for random length increments
             params.Li = params_Li[[x]],      # distributions parameters fore random length increments
             cor.s = cor_s[[x]],              # Correlation struction for vital rates
             names.s = names_s[[x]],          # names for vital rates
             dist.s = dist_s[[x]],            # distributions names for random  vital rates
             params.s = params_s.1[[x]],      # distributions parameters fore random vital rate
             L.mean = L_mean[[x]],            # Initial (mean) length at age  
             p.rep = p.rep,                   # Proportion mature by age class
             p.cat = p.cat,                   # Probability of castastrophe
             e.cat = e.cat,                   # effect on castastrophe on pop size (percent recuction)
             gen.time = generation_time[[x]], # Generation Time
             data = c(lh_data[[x]], l_min$equil[[x]], s0 = lh_data[[x]]$syoy.1)# Life histoty data
      )
  })
})

stopCluster(cl)

save(MVP_run, file = "Results 03 all - MVP population.Rdata")

rm(MVP_run)

#-------------------------------------------------------------------------------
# MVP Results

load("Results 03 all - MVP population.Rdata")

# Extinction Probability
time_frame =  250
PVA.P <- lapply(1:6, function(x){
  A <- lapply(1:length(adultPop), function(a){
    R <- lapply(1:maxRuns, function(y){
      as.data.frame(list(
        N = adultPop[[a]],         # inital pop size
        rep = 1:numReps,           # replicate identifier
        run = y,                   # run identifier
        P = is_ext_list(MVP_run[[x]][[a]][[y]], threshold = 25, years = time_frame), # extinction boolean
        DU = DU[[x]],              # DU   
        cat = "0.15/Generation"    # Catastrophe rate
      ))
    })
    do.call(rbind, R)
  })
  do.call(rbind, A)
})

mod.PVA <- lapply(PVA.P, function(x) glm(P ~ log10(N), data = x, family = binomial))

# Time to Extinction
ext_time <- function(x, threshold){
  ext <- which(x <= threshold)
  if(length(ext) < 1){
    ext = length(x)
  }
  min(ext)
}

ext_time_list <- function(x, threshold){
  sapply(x, ext_time, threshold = threshold)
}

PVA.T <- lapply(1:6, function(x){
  A <- lapply(1:length(adultPop), function(a){
      
    as.data.frame(list(
      N = adultPop[[a]],         # inital pop size
      rep = 1:numReps,           # replicate identifier
      Time = ext_time_list(unlist(MVP_run[[x]][[a]], recursive = F), threshold = 25), # extinction time
      DU = DU[[x]],              # DU   
      cat = "0.15/Generation"    # Catastrophe rate
    ))
  })
  do.call(rbind, A)
})

# clean up
rm(MVP_run)

#-------------------------------------------------------------------------------
# Plot

# Subset to reasonable probability range for fitting NLS curve 
plot_data <- lapply(1:6, function(x){
  A <- lapply(1:length(adultPop), function(a){
    R <- lapply(1:maxRuns, function(y){
      data <- PVA.P[[x]][PVA.P[[x]]$N == adultPop[[a]] & PVA.P[[x]]$run == y,]
      
      as.data.frame(list(
        N = adultPop[[a]],         # inital pop size
        run = y,                   # run identifier
        P = sum(data$P)/numReps, # extinction boolean
        DU = DU[[x]],              # DU   
        cat = "0.15/Generation"    # Catastrophe rate
      ))
    })
    do.call(rbind, R)
  })
  do.call(rbind, A)
})

# adult pop size
Na <- seq(min(adultPop), max(adultPop) + 100, 10)

pred_data <- as.data.frame(do.call(rbind, 
      lapply(1:6, function(x){
        as.data.frame(list(
          N = Na,
          rep = NA,
          P = predict(mod.PVA[[x]], list(N = Na), type = "response"),
          DU = DU[[x]],
          cat = "0.15/Generation"
        ))
      })
))

# Plot - adult pop size vs. prob of extirpation
p <- ggplot(do.call(rbind, plot_data), aes(N, P))
p <- p + geom_hline(aes(yintercept = 0.01), colour = "black", linetype = "dashed")
p <- p + geom_hline(aes(yintercept = 0.05), colour = "black", linetype = "dashed")
p <- p + geom_point(aes(), size = 1.5)
p <- p + geom_path(data = pred_data, aes(N, P))
p <- p + scale_x_continuous(trans = "log10")
p <- p + facet_wrap(~ as.factor(DU), ncol = 3)
p <- p + labs(x ="Adult population size ", 
              y = "Probability of Extirpation") 
p <- p + theme_me +
     theme(strip.text = element_text(size = 11, family = "sans", face = "bold"),
           legend.justification = c(1, 1), legend.position = c(0.99, 0.99))

#png("Figures/MVP.png", width = 8, height = 4,5, unit = "in", res = 300)
p
#dev.off()

# Time to Exticntion PDF/CDF
PVA.T.pdf <- lapply(1:6, function(x){
  A <- lapply(1:length(adultPop), function(a){
    
    data <- PVA.T[[x]][PVA.T[[x]]$N == adultPop[[a]],]
    
    pdf <- NULL
    maxYears <- 500
    for(i in 1:maxYears){
      pdf <- c(pdf, length(which(data$Time == i)) / (numReps * maxRuns))
    }
    cdf <- cumsum(pdf)
    
    as.data.frame(list(
      N = adultPop[[a]],      # inital pop size
      year = 1:maxYears,      # years
      PDF = pdf,              # extinction at year probabilty density function
      CDF = cdf,              # extinction at year cummulativedensity function
      DU = DU[[x]],           # DU   
      cat = "0.15/Generation" # Catastrophe rate
    ))
  })
  do.call(rbind, A)
})

# Plot - adult pop size vs. prob of extirpation


a <- c(50, 100, 250, 500, 1000, 5000)
plot_data <- do.call(rbind, PVA.T.pdf)
plot_data <- plot_data[plot_data$N %in% a,]
p <- ggplot(plot_data , aes(year, CDF))
p <- p + geom_path(aes(colour = as.factor(N)))
p <- p + facet_wrap(~ as.factor(DU), ncol = 3)
p <- p + scale_y_continuous(breaks = seq(0, 1, 0.1))
p <- p + labs(x ="Year", y = "Cumulative Extinction Probability", colour = "Female \nPopuation \nSize") 
p <- p + theme_me +
  theme(strip.text = element_text(size = 11, family = "sans", face = "bold"),
        axis.text.y = element_text(angle = 0))

#png("Figures/MVP - CDF.png", width = 8, height = 4.5, unit = "in", res = 300)
p
#dev.off()


CDF250 <- do.call(rbind, lapply(1:6, function(x){
  PVA.T.pdf[[x]][PVA.T.pdf[[x]]$year == 250,]
}))

p <- ggplot(CDF250 , aes(N, CDF))
p <- p + geom_path(aes(colour = DU))
p <- p + scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 500, 1000, 5000, 10000, 150000))
p <- p + scale_y_continuous(breaks = seq(0,1,.1))
p <- p + labs(x ="Female pop size", y = "Cumulative Extinction Probability", colour = "DU") 
p <- p + theme_me

#png("Figures/MVP - CDF.png", width = 8, height = 4.5, unit = "in", res = 300)
p
#dev.off()



#-------------------------------------------------------------------------------
# MVP Estimate
# Estiamted MVP with 5% and 1% change of extinction for all MVP and Catastrophe
# Scenarios

MVP <- function(P, a, b) 10^(-((log(1/P - 1) +a)/b))

Ps <- c(0.05, 0.01)

MVP_est.A <- lapply(1:6, function(x) {
  mvp <- sapply(Ps, function(p) {
     MVP(p, a = coef(mod.PVA[[x]])[[1]], b = coef(mod.PVA[[x]])[[2]])
  })
  names(mvp) = paste0("P",Ps)
  mvp
})
names(MVP_est.A) <- DU

mod_data <- as.data.frame(list(
  MVP = sapply(MVP_est.A, function(x) x[[2]]),
  gen_time = unlist(generation_time)
))

mod <- lm(MVP ~ gen_time, data = mod_data)

p <- ggplot(mod_data , aes(gen_time, MVP))
p <- p + geom_point()
p <- p + geom_abline(intercept = coef(mod)[1], slope = coef(mod)[2])
p <- p + labs(x ="Generation Time (years)", y = "MVP (adult females)") 
p <- p + theme_me

#png("Figures/MVP - generation time.png", width = 4, height = 3, unit = "in", res = 300)
p
#dev.off()

#-------------------------------------------------------------------------------
# MVP pop size for younger stage

# Calc pop size for given adult pp size

# YOY MVP
MVP_est.Y <- lapply(1:6, function(x) {
  sapply(MVP_est.A[[x]], function(pop) {
    stage_pop(pmx = pmx.d[[x]], stage = "YOY", pop = pop, 
              pop.stage = "adult", p.rep = p.rep)
  })
})
names(MVP_est.Y) <- DU

# Juvenile MVP
MVP_est.J <- lapply(1:6, function(x) {
  sapply(MVP_est.A[[x]], function(pop) {
    stage_pop(pmx = pmx.d[[x]], stage = "Juvenile", pop = pop, 
              pop.stage = "adult", p.rep = p.rep)
  })
})
names(MVP_est.J) <- DU

# age 1+ MVP
MVP_est.1plus <- Map("+", MVP_est.J, MVP_est.A)

# Total - all stage
MVP_est.T <- Map("+", MVP_est.Y, MVP_est.1plus)

# > 1000 mm
# Many est. of lake sturgeon pop size are
MVP_1000_est <- function(Pop, age.min, age.max, data) {
  
  t = age.min:age.max
  L = with(data, VB.growth(t, Linf, K, L0 = H0))
  
  m0 <- data$M * L[t == 26]
  
  s <- NULL
  for(i in t){
    if(i < 1){
      s <- c(s, data$syoy.1d)
    } else if(i < data$c.MO ){
      s <- c(s, (L[i+1] / L[i+2] * exp(-data$K)) ^ (m0 / data$K / data$Linf))
    } else if(i >= data$c.MO ){
      s <- c(s, exp(-data$M))
    }
  }
  
  Pop <- Pop * cumprod(s)
  sum(Pop[L >= 1000])
}

MVP_est.1000 <- lapply(1:6, function(x) {
  sapply(MVP_est.Y[[x]], MVP_1000_est,  0, 100, lh_data[[x]])
})

# Male pop size
MVP_male_est <- function(Pop, age.min, age.max, data){
  
  t = age.min:age.max
  L = with(data, VB.growth(t, Linf, K, L0 = H0))
  
  # Survival
  m0 <- data$M * L[t == 26]
  
  st <- NULL
  for(i in t){
    if(i < 1){
      st <- c(st, data$syoy.1d)
    } else if(i < data$c.MO ){
      st <- c(st, (L[i+1] / L[i+2] * exp(-data$K)) ^ (m0 / data$K / data$Linf))
    } else if(i >= data$c.MO ){
      st <- c(st, exp(-data$M))
    }
  }
  
  # Sex Ratio
  age = c(0, 25, 40, 80)   # Age
  PM = c(0.5, 1/3, 1/7, 0) # Proportion Male
  
  R <- NULL
  for(i in t){
    if(i < age[2]){
      s <- (PM[2] - PM[1])/(age[2] - age[1])
      R <- c(R, (i-age[1])*s + PM[1])
    } else if(i < age[3]){
      s <- (PM[3] - PM[2])/(age[3] - age[2])
      R <- c(R, (i-age[2])*s + PM[2])
    } else if(i < age[4]){
      s <- (PM[4] - PM[3])/(age[4] - age[3])
      R <- c(R, (i-age[3])*s + PM[3])
    } else if(i >= age[4]){
      R <- c(R, 0)
    }
  }
  R <- R/0.5
  
  Pop.F <- Pop * cumprod(st)
  Pop.M <- Pop.F * R
  
  list(YOY = Pop,
       Juvenile = sum(Pop.M[t > 0 & t < 20]),
       L1000 = sum(Pop.M[L >= 1000]),
       Adult = sum(Pop.M[t >= 20])
       
  )
}

MVP_est.Male <- lapply(1:6, function(x) {
  sapply(MVP_est.Y[[x]], MVP_male_est,  0, 100, lh_data[[x]])
})

save(MVP_est.A, MVP_est.J, MVP_est.Y, MVP_est.1plus, MVP_est.T, MVP_est.1000, MVP_est.Male, 
     file = "Results 03a - MVP estimates.Rdata")

#-------------------------------------------------------------------------------
# Output results

# MVP RESULTS
MVP_est.Y <- lapply(1:6, function(x) {
  data <- as.data.frame(list(P.05 = MVP_est.Y[[x]][[1]],
                             P.01 = MVP_est.Y[[x]][[2]],
                             Stage = "YOY",
                             Sex = "F",
                             DU = DU[[x]])) 
  data
})
MVP_est.J <- lapply(1:6, function(x) {
  data <- as.data.frame(list(P.05 = MVP_est.J[[x]][[1]],
                             P.01 = MVP_est.J[[x]][[2]],
                             Stage = "Juvenile",
                             Sex = "F",
                             DU = DU[[x]])) 
  data
})
MVP_est.A <- lapply(1:6, function(x) {
  data <- as.data.frame(list(P.05 = MVP_est.A[[x]][[1]],
                             P.01 = MVP_est.A[[x]][[2]],
                             Stage = "Adult",
                             Sex = "F",
                             DU = DU[[x]])) 
  data
})
MVP_est.1000 <- lapply(1:6, function(x) {
  data <- as.data.frame(list(P.05 = MVP_est.1000[[x]][[1]],
                             P.01 = MVP_est.1000[[x]][[2]],
                             Stage = "1000mm",
                             Sex = "F",
                             DU = DU[[x]])) 
  data
})

MVP_est.Male <- lapply(1:6, function(x) {
  data <- as.data.frame(list(P.05 = unlist(MVP_est.Male[[x]][,1]),
                             P.01 = unlist(MVP_est.Male[[x]][,2]),
                             Stage = names(MVP_est.Male[[x]][,1]),
                             Sex = "M",
                             DU = DU[[x]])) 
  data
})

MVPres <- rbind(do.call(rbind,MVP_est.Y),
                do.call(rbind,MVP_est.J),
                do.call(rbind,MVP_est.A),
                do.call(rbind,MVP_est.1000),
                do.call(rbind,MVP_est.Male))

write.table(MVPres, "Results 04 - MVP estimates.csv", sep = ",", 
            col.names = TRUE, row.names = FALSE)

#-------------------------------------------------------------------------------
# MAPV
#-------------------------------------------------------------------------------
# Minimum area for population viability - in hectares
# Estimated simply as: MAPVi = MVPi * APIi
# This produces over (or under) estimate of MVP value because
# fish growth and death is non-linear throughout the year. 
# With estimated of daily growth and survival the stage-specific date of maximum
# area usage is estiamted (d_area) 
#-------------------------------------------------------------------------------
# MAPV based on Lake Sturgeon densities
# Lake sturgeon specific density allometries
# compile lake sturgeon density at mean weight data, assume slopes from Randall
# et al. (1995) and fit two intercepts representin median and maximum densities
# use to estimate mean and minimum area requirements per individual

density <- read.csv("data/Density.csv", header = TRUE)
density <- density[density$Weight > 0,]

# Density Allometry Parameters - Randall et al. (1995)
# D = 30199.52*W^-1.01 (Lake)
# D = 78432.82*W^-0.94  (river)
parR <- c(78432.82, -0.94) # River density parameters
parL <- c(30199.52, -1.01) # Lake density parameters

# Density per unit weight
D1_L <- with(density[density$Habitat == "Lake",], Density/(Weight*1000) ^ parL[2])  # Lake
D1_R <- with(density[density$Habitat == "River",], Density/(Weight*1000) ^ parR[2]) # River

# Intercept values
mid_D_I.L <- median(D1_L) # mean Lake density intercept
mid_D_I.R <- median(D1_R) # mean Riverdensity intercept
max_D_I.L <- max(D1_L)    # max Lake density intercept
max_D_I.R <- max(D1_R)    # max River density intercept

# Density allometry parameters - Lake sturgeon specific
parL.mid <- c(mid_D_I.L, -1.01) # median density - Lake
parR.mid <- c(mid_D_I.R, -0.94) # median density - River
parL.max <- c(max_D_I.L, -1.01) # maximum density - Lake
parR.max <- c(max_D_I.R, -0.94) # maximum density - River

# clean up
rm(D1_L, D1_R, mid_D_I.L, mid_D_I.R, max_D_I.L, max_D_I.R)

# Predict relationships
W <- seq(2, 30, 0.1)
pred_data.L <- as.data.frame(list(
  Weight = rep(W, 3),
  Density = c(parL.mid[[1]] * (W * 1000) ^ parL.mid[[2]],
              parL.max[[1]] * (W * 1000) ^ parL[2],
              parL[1] * (W * 1000) ^ parL[2]), 
  Habitat = "Lake",
  Relationship = c(rep("Median", length(W)),
                   rep("Maximum", length(W)),
                   rep("Randall", length(W)))
))
pred_data.R <- as.data.frame(list(
  Weight = rep(W, 3),
  Density = c(parR.mid[[1]] * (W * 1000) ^ parR.mid[[2]],
              parR.max[[1]] * (W * 1000) ^ parR[2],
              parR[1] * (W * 1000) ^ parR[2]), 
  Habitat = "River",
  Relationship = c(rep("Median", length(W)),
                   rep("Maximum", length(W)),
                   rep("Randall", length(W)))
))
pred_data <- rbind(pred_data.L, pred_data.R)
rm(pred_data.L, pred_data.R)
pred_data$Relationship <- factor(pred_data$Relationship, levels = c("Randall", "Maximum", "Median"))

# PLot 
p <- ggplot(density, aes(Weight, Density))
p <- p + geom_path(data = pred_data, aes(Weight, Density, colour = Relationship), size = 1.25)
p <- p + geom_point(size = 1)
p <- p + scale_x_log10(breaks = c(2, 6, 10, 20))
p <- p + scale_y_log10(minor_breaks = c(0.5, 5))
p <- p + facet_wrap(~ Habitat, ncol = 2)
p <- p + labs(x = "Weight (kg)", y = "Density (fish/ha)")               
p <- p + theme_me 

#png("Figures/Density Relationship.png", width = 7, height = 3, unit = "in", res = 300)
p
#dev.off()

#clean up
rm(density, W, pred_data, p)

#-------------------------------------------------------------------------------
# MAPV Estiamtes

# Load MVP estimates
load("Results 03a - MVP estimates.Rdata")

# Daily area requirements of a Cohort

# Randall
MAPV_day.L <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                           M.e = 3, data = lh_data[[x]], parD = parL))
MAPV_day.R <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                             M.e = 3, data = lh_data[[x]], parD = parR))

# Median
MAPV_day_mid.L <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                                 M.e = 3, data = lh_data[[x]], parD = parL.mid))
MAPV_day_mid.R <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                                 M.e = 3, data = lh_data[[x]], parD = parR.mid))

# Minimum
MAPV_day_min.L <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                                 M.e = 3, data = lh_data[[x]], parD = parL.max))
MAPV_day_min.R <- lapply(1:6, function(x) lapply(MVP_est.Y[[x]][c(1,2)], MAPVd_cohort_w, age.min = 0, age.max = 100, 
                                                 M.e = 3, data = lh_data[[x]], parD = parR.max))

# Estimate MAPV value for various Lake Sturgeon population sizez
MAPV.est <- function(Pop.f, pmx, MAPVd, data, D = "Randall", M.e = 3) {
  
  YOY.F <- stage_pop(pmx = pmx, stage = "YOY", pop = Pop.f, 
                     pop.stage = "adult", p.rep = p.rep)
  
  pop.F <- (YOY.F * exp(-M.e)) * MAPVd$Survival
  pop.M <- pop.F * MAPVd$R 
  
  dMAPV <- MAPVd$API * (pop.F + pop.M)
  
  if(D == "Randall"){
    # Total area usage - all Ages
    MAPV <- aggregate(dMAPV, by = list(MAPVd$jDay), sum)
    MAPV <- max(MAPV[,2])
  } else if(D == "Density"){
    # > 1000 mm
    MAPV <- aggregate(dMAPV[MAPVd$L >= 1000], by = list(MAPVd$jDay[MAPVd$L >= 1000]), sum) # Area
    MAPV<- max(MAPV[,2])
  }
  
  MAPV
}

Pop <- seq(100, 10000, 100)
p.rep <- c(rep(0, 3), rep(1, 2))

# Randall
MAPV_pop.L <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
   MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day.L[[x]][[1]], data = lh_data[[x]], D = "Randall")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "Lake",
    DU = DU[[x]],
    Density = "Randall"
  ))
})
MAPV_pop.R <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
    MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day.R[[x]][[1]], data = lh_data[[x]], D = "Randall")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "River",
    DU = DU[[x]],
    Density = "Randall"
  ))
})

# Median
MAPV_pop_mid.L <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
    MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day_mid.L[[x]][[1]], data = lh_data[[x]], D = "Density")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "Lake",
    DU = DU[[x]],
    Density = "Median"
  ))
})
MAPV_pop_mid.R <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
    MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day_mid.R[[x]][[1]], data = lh_data[[x]], D = "Density")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "River",
    DU = DU[[x]],
    Density = "Median"
  ))
})

# Minimum
MAPV_pop_min.L <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
    MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day_min.L[[x]][[1]], data = lh_data[[x]], D = "Density")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "Lake",
    DU = DU[[x]],
    Density = "Maximum"
  ))
})
MAPV_pop_min.R <- lapply(1:6, function(x) {
  MAPV = sapply(Pop, function(p){
    MAPV.est(Pop.f = p, pmx = pmx.d[[x]], MAPVd = MAPV_day_min.R[[x]][[1]], data = lh_data[[x]], D = "Density")
  })
  as.data.frame(list(
    Pop = Pop,
    MAPV = MAPV,
    Habitat = "River",
    DU = DU[[x]],
    Density = "Maximum"
  ))
})

# Day of max area for cohort
t_area.L <- sapply(MAPV_day.L, function(x) x[[2]]$Age[which.max(x[[2]]$MAPV.T)])
t_area.R <- sapply(MAPV_day.R, function(x) x[[2]]$Age[which.max(x[[2]]$MAPV.T)])

# Stage-specific daily area requirements
# aggregate MAPV_day by day of year within each stage
MAPV_jDay.L <- lapply(1:6, function(x) lapply(MAPV_day.L[[x]], MAPV_jDay_func, data = lh_data[[x]]))
MAPV_jDay.R <- lapply(1:6, function(x) lapply(MAPV_day.R[[x]], MAPV_jDay_func, data = lh_data[[x]]))

MAPV_jDay_mid.L <- lapply(1:6, function(x) lapply(MAPV_day_mid.L[[x]], MAPV_jDay_func, data = lh_data[[x]]))
MAPV_jDay_mid.R <- lapply(1:6, function(x) lapply(MAPV_day_mid.R[[x]], MAPV_jDay_func, data = lh_data[[x]]))

MAPV_jDay_min.L <- lapply(1:6, function(x) lapply(MAPV_day_min.L[[x]], MAPV_jDay_func, data = lh_data[[x]]))
MAPV_jDay_min.R <- lapply(1:6, function(x) lapply(MAPV_day_min.R[[x]], MAPV_jDay_func, data = lh_data[[x]]))

# Stage-specific maximum area requirements
# Find that maximum stage-specific MAPV over the year
MAPV_stage.L <- lapply(MAPV_jDay.L, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))
MAPV_stage.R <- lapply(MAPV_jDay.R, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))

MAPV_stage_mid.L <- lapply(MAPV_jDay_mid.L, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))
MAPV_stage_mid.R <- lapply(MAPV_jDay_mid.R, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))

MAPV_stage_min.L <- lapply(MAPV_jDay_min.L, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))
MAPV_stage_min.R <- lapply(MAPV_jDay_min.R, function(MAPVa) lapply(MAPVa, function(x){
  MAPV <- aggregate(x$MAPV, by = list(x$Stage), max)
  colnames(MAPV) = c("Stage", "MAPV")
  MAPV
}))


# calc stage-specific day at max area used - whcih day does MAPV_Stage occur
d_area.L <- lapply(MAPV_jDay.L, function(MAPVd) lapply(MAPVd, function(x) {
  d_area <- NULL
  for(i in unique(x$Stage)){
    max <- which.max(x$MAPV[x$Stage == i])
    d_area <- c(d_area, x$jDay[x$Stage == i][max])
  }
  as.data.frame(list(Stage = unique(x$Stage), d_area = d_area))
}))
d_area.R <- lapply(MAPV_jDay.R, function(MAPVd) lapply(MAPVd, function(x) {
  d_area <- NULL
  for(i in unique(x$Stage)){
    max <- which.max(x$MAPV[x$Stage == i])
    d_area <- c(d_area, x$jDay[x$Stage == i][max])
  }
  as.data.frame(list(Stage = unique(x$Stage), d_area = d_area))
}))

#-------------------------------------------------------------------------------
# PLOTS 

# Plot - Cohort area used - over life time
MAPV_day.L <- lapply(MAPV_day.L, function(x) lapply(1:2, function(i) {x[[i]]$Habitat = "Lake"; x[[i]]}))
MAPV_day.R <- lapply(MAPV_day.R, function(x) lapply(1:2, function(i) {x[[i]]$Habitat = "River"; x[[i]]}))
MAPV_day.L <- lapply(1:6, function(x) lapply(1:2, function(i) {MAPV_day.L[[x]][[i]]$DU = DU[[x]]; MAPV_day.L[[x]][[i]]}))
MAPV_day.R <- lapply(1:6, function(x) lapply(1:2, function(i) {MAPV_day.R[[x]][[i]]$DU = DU[[x]]; MAPV_day.R[[x]][[i]]}))

p <- ggplot(rbind(do.call(rbind, lapply(MAPV_day.L, function(x) x[[2]])),
                  do.call(rbind, lapply(MAPV_day.R, function(x) x[[2]]))),
            aes(Age, MAPV.T))
p <- p + geom_path(aes(colour = as.factor(DU)), size = 1)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + labs(x = "Age", y = "Area Occupied", colour = "DU")               # Axis names
p <- p + facet_wrap(~(Habitat), ncol = 2)
p <- p + theme_me #+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

#png("Figures/MAPV - cohort.png", width = 7, height = 3, unit = "in", res = 300)
p
#dev.off()

# Plot - Stage-specific area usage - over year
MAPV_jDay.L <- lapply(MAPV_jDay.L, function(x) lapply(1:2, function(i) {x[[i]]$Habitat = "Lake"; x[[i]]}))
MAPV_jDay.R <- lapply(MAPV_jDay.R, function(x) lapply(1:2, function(i) {x[[i]]$Habitat = "River"; x[[i]]}))
MAPV_jDay.L <- lapply(1:6, function(x) lapply(1:2, function(i) {MAPV_jDay.L[[x]][[i]]$DU = DU[[x]]; MAPV_jDay.L[[x]][[i]]}))
MAPV_jDay.R <- lapply(1:6, function(x) lapply(1:2, function(i) {MAPV_jDay.R[[x]][[i]]$DU = DU[[x]]; MAPV_jDay.R[[x]][[i]]}))

p <- ggplot(rbind(do.call(rbind, lapply(2, function(x) MAPV_jDay.L[[x]][[2]])),
                  do.call(rbind, lapply(2, function(x) MAPV_jDay.R[[x]][[2]]))),
            aes(jDay, MAPV))
p <- p + geom_path(aes(colour = Stage), size = 1)
p <- p + labs(x = "Day", y = "Area Occupied") 
p <- p + facet_wrap(DU~(Habitat), ncol = 2)
#p <- p + scale_y_log10()
p <- p + theme_me #+
  #theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

#png("Figures/MAPV - daily.png", width = 7, height = 3, unit = "in", res = 300)
p
#dev.off()

# Plot - MAPV at diff pop size
MAPV_pop.L <- do.call(rbind, MAPV_pop.L)
MAPV_pop.R <- do.call(rbind, MAPV_pop.R)
MAPV_pop_mid.L <- do.call(rbind, MAPV_pop_mid.L)
MAPV_pop_mid.R <- do.call(rbind, MAPV_pop_mid.R)
MAPV_pop_min.L <- do.call(rbind, MAPV_pop_min.L)
MAPV_pop_min.R <- do.call(rbind, MAPV_pop_min.R)

p <- ggplot(rbind(MAPV_pop.L, MAPV_pop.R, MAPV_pop_mid.L, MAPV_pop_mid.R, MAPV_pop_min.L, MAPV_pop_min.R),
            aes(Pop, MAPV))
p <- p + geom_path(aes(colour = as.factor(DU)), size = 1)
p <- p + scale_x_continuous(trans = "log10",
                            minor_breaks = c(seq(20, 90, 10), seq(200, 900, 100), seq(2000, 9000, 1000)))
p <- p + scale_y_continuous(trans = "log10",
                            minor_breaks = c(seq(20, 90, 10), seq(200, 900, 100), seq(2000, 9000, 1000),
                                             seq(20000, 90000, 10000)))
p <- p + labs(x = "Female Population Size", y = "Area Occupied", colour = "DU")               # Axis names
p <- p + facet_wrap(Density~Habitat, ncol = 2, scales = "free_y")
p <- p + theme_me 
#p <- p + guides(colour = guide_legend(nrow = 1))

#png("Figures/MAPV - MAPV vs pop.png", width = 6.5, height = 7.5, unit = "in", res = 300)
p
#dev.off()


#-------------------------------------------------------------------------------
# Output results

# MAPV RESULTS
# Females and Males spatial requirements in hectares
MAPV_stage <- lapply(1:6, function(x) {
  L <- as.data.frame(list(Stage = MAPV_stage.L[[x]][[1]]$Stage,
                          P.05 = MAPV_stage.L[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage.L[[x]][[2]]$MAPV,
                          Habitat = "Lake",
                          DU = DU[[x]],
                          Estimate = "Randall")) 
  R <- as.data.frame(list(Stage = MAPV_stage.R[[x]][[1]]$Stage,
                          P.05 = MAPV_stage.R[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage.R[[x]][[2]]$MAPV,
                          Habitat = "River",
                          DU = DU[[x]],
                          Estimate = "Randall"))  
  rbind(L, R)
})
MAPV_stage <- do.call(rbind, MAPV_stage)

MAPV_stage_mid <- lapply(1:6, function(x) {
  L <- as.data.frame(list(Stage = MAPV_stage_mid.L[[x]][[1]]$Stage,
                          P.05 = MAPV_stage_mid.L[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage_mid.L[[x]][[2]]$MAPV,
                          Habitat = "Lake",
                          DU = DU[[x]],
                          Estimate = "Median")) 
  R <- as.data.frame(list(Stage = MAPV_stage_mid.R[[x]][[1]]$Stage,
                          P.05 = MAPV_stage_mid.R[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage_mid.R[[x]][[2]]$MAPV,
                          Habitat = "River",
                          DU = DU[[x]],
                          Estimate = "Median")) 
  rbind(L, R)
})
MAPV_stage_mid <- do.call(rbind, MAPV_stage_mid)

MAPV_stage_min <- lapply(1:6, function(x) {
  L <- as.data.frame(list(Stage = MAPV_stage_min.L[[x]][[1]]$Stage,
                          P.05 = MAPV_stage_min.L[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage_min.L[[x]][[2]]$MAPV,
                          Habitat = "Lake",
                          DU = DU[[x]],
                          Estimate = "Minimum")) 
  R <- as.data.frame(list(Stage = MAPV_stage_min.R[[x]][[1]]$Stage,
                          P.05 = MAPV_stage_min.R[[x]][[1]]$MAPV,
                          P.01 = MAPV_stage_min.R[[x]][[2]]$MAPV,
                          Habitat = "River",
                          DU = DU[[x]],
                          Estimate = "Minimum")) 
  rbind(L, R)
})
MAPV_stage_min <- do.call(rbind, MAPV_stage_min)

write.table(rbind(MAPV_stage, MAPV_stage_mid, MAPV_stage_min), 
            "Results 05 - MAPV estimates.csv", sep = ",", 
            col.names = TRUE, row.names = FALSE)

################################################################################
#-------------------------------------------------------------------------------
# MAPV ALTERNATIVE
#-------------------------------------------------------------------------------

# D = 30199.52*W^-1.01 (Lake)
# D = 78432.82*W^-0.94  (river)

parR <- c(78432.82, -0.94) # River density parameters
parL <- c(30199.52, -1.01) # Lake density parameters
parL.mid <- c(mid_D_I.L, -1.01)
parR.mid <- c(mid_D_I.L, -0.94)

parL.max <- c(max_D_I.L, -1.01)
parR.max <- c(max_D_I.L, -0.94)

MAPVd_w.L <- lapply(1:6, function(x){
            lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
                   data = lh_data[[x]], parD = parL)
})
MAPV_w.L <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w.L[[x]][[y]]$MAPV, max)
}))

MAPVd_w.R <- lapply(1:6, function(x){
  lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
         data = lh_data[[x]], parD = parR)
})
MAPV_w.R <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w.R[[x]][[y]]$MAPV, max)
}))

MAPVd_w_min.L <- lapply(1:6, function(x){
  lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
         data = lh_data[[x]], parD = parL.max)
})
MAPV_w_min.L <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w_min.L[[x]][[y]]$MAPV, max)
}))

MAPVd_w_min.R <- lapply(1:6, function(x){
  lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
         data = lh_data[[x]], parD = parR.max)
})
MAPV_w_min.R <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w_min.R[[x]][[y]]$MAPV, max)
}))

MAPVd_w_mid.L <- lapply(1:6, function(x){
  lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
         data = lh_data[[x]], parD = parL.mid)
})
MAPV_w_mid.L <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w_mid.L[[x]][[y]]$MAPV, max)
}))

MAPVd_w_mid.R <- lapply(1:6, function(x){
  lapply(MVP_est.Y[[x]], MAPVd_stage_w, age.min = 0, age.max = 100, 
         data = lh_data[[x]], parD = parR.mid)
})
MAPV_w_mid.R <- lapply(1:6, function(x) lapply(1:2, function(y){ 
  sapply(MAPVd_w_mid.R[[x]][[y]]$MAPV, max)
}))

# Plot - Stage-specific area usage - over year



Stage = c("YOY", "Juv", "Adt", "1Plus", "All")
plot_data <- {
  
  D <- lapply(1:6, function(x) {
    S <- lapply(1:length(Stage), function(y){
      as.data.frame(list(
        Day = 1:365,
        MAPV = MAPVd_w.L[[x]][[2]]$MAPV[[y]],
        Stage = Stage[[y]],
        DU = DU[[x]],
        Habitat = "Lake",
        P = 0.15
      ))
    })
    do.call(rbind, S)
  })
  L <-  do.call(rbind, D)
  
  D <- lapply(1:6, function(x) {
    S <- lapply(1:length(Stage), function(y){
      as.data.frame(list(
        Day = 1:365,
        MAPV = MAPVd_w.R[[x]][[2]]$MAPV[[y]],
        Stage = Stage[[y]],
        DU = DU[[x]],
        Habitat = "River",
        P = 0.15
      ))
    })
    do.call(rbind, S)
  })
  R <-  do.call(rbind, D)
  
  rbind(L, R)
  
}


p <- ggplot(plot_data[plot_data$DU == "DU3",], aes(Day, MAPV))
p <- p + geom_path(aes(colour = Stage), size = 1)
p <- p + labs(x = "Day", y = "Area Occupied") 
p <- p + facet_wrap(DU~(Habitat), ncol = 2)
#p <- p + scale_y_log10()
p <- p + theme_me #+
#theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

#png("Figures/MAPV - daily.png", width = 7, height = 3, unit = "in", res = 300)
p
#dev.off()

#-------------------------------------------------------------------------------

# MAPV Results - Minns allometry
minnsMAPVres <- as.data.frame(list(
  MAPV = as.vector(c(unlist(MAPV.Y), unlist(MAPV.J), unlist(MAPV.A), 
                     unlist(MAPV.1plus), unlist(MAPV.T))),
  Stage = c(rep("YOY", 12), rep("Juvenile", 12), rep("Adult", 12),
            rep("1+", 12), rep("All", 12)),
  Cor = rep(c(rep("ind", 2), rep("mid", 2), rep("str", 2)), 2),
  Cat = c(rep(0.1, 6), rep(0.15, 6)),
  P = rep(c(0.05, 0.01), 6)
))

write.table(minnsMAPVres, "Results 05b - MAPV estimates - Minns Allometry.csv", sep = ",", 
            col.names = TRUE, row.names = FALSE)

# MAPV Results - Mean redside dace density
meanMAPVres <- as.data.frame(list(
  MAPV = c(unlist(MAPV.Y_mean), unlist(MAPV.J_mean), unlist(MAPV.A_mean), 
           unlist(MAPV.1plus_mean), unlist(MAPV.T_mean)),
  Stage = c(rep("YOY", 12), rep("Juvenile", 12), rep("Adult", 12),
            rep("1+", 12), rep("All", 12)),
  Cor = rep(c(rep("ind", 2), rep("mid", 2), rep("str", 2)), 2),
  Cat = c(rep(0.1, 6), rep(0.15, 6)),
  P = rep(c(0.05, 0.01), 6)
))

write.table(meanMAPVres, "Results 06b - MAPV estimates - Mean Density.csv", sep = ",", 
            col.names = TRUE, row.names = FALSE)

# MAPV Results - Mi redside dace density
minMAPVres <- as.data.frame(list(
  MAPV = c(unlist(MAPV.Y_min), unlist(MAPV.J_min), unlist(MAPV.A_min), 
           unlist(MAPV.1plus_min), unlist(MAPV.T_min)),
  Stage = c(rep("YOY", 12), rep("Juvenile", 12), rep("Adult", 12),
            rep("1+", 12), rep("All", 12)),
  Cor = rep(c(rep("ind", 2), rep("mid", 2), rep("str", 2)), 2),
  Cat = c(rep(0.1, 6), rep(0.15, 6)),
  P = rep(c(0.05, 0.01), 6)
))

write.table(minMAPVres, "Results 07b - MAPV estimates - Min Density.csv", sep = ",", 
            col.names = TRUE, row.names = FALSE)

#===============================================================================
# RECOVERY PROJECTION
#===============================================================================
# predict time to recovery (MVP)  for various population sizes at different 
# population growth rates
# Begin by solving for syoy value to produce different growth rates (1.01:1.15)
# Run simulations similar to MVP to predict population size over time
# find the year when recovery occurs (Na >= MVP) for each simulation
# Find the year when there is a 95% probabilyt recovery will have occured
#-------------------------------------------------------------------------------
# YOY survival rate for each lambda

# Lambda values
lambdas <- seq(1.01, 1.15, 0.02)

# Survival rate stochastic vital rate - exlcuding YOY
params_s <- lapply(1:6, function(x) { # Survival
  lapply(1:(lh_data[[x]]$Tmax), function(i) beta_stretch_val(M_means[[x]][[i]], sd = (M_means[[x]][[i]] * lh_data[[x]]$cv.M)))
})

# OPTIMIZATION
if(FALSE){
  # Function to optimize
  M_optim_f = function (s0, target.lambda, numReps, years, pmxs, data){
    
    # mean projection matrix
    pmx_mean <- pmx_create(c(data, s0 = s0))
    
    # initial population size
    N0 <- initPop(pmx_mean, 100, "adult", p.rep = c(0, 0, 0, 1, 1))
    
    # Mean annual survival rates by age
    M0 <- -log(s0)
    
    # Vital rates
    params.s0 <- beta_stretch_val(M0, sd = (M0 * data$cv.M))
    
    # Generate random YOY survival rate
    s0.rand <- exp(-with(params.s0, rbeta(years * numReps, shape1 = a, shape2 = b) * (max - min) + min))
    s0.rand <- matrix(s0.rand, nrow = numReps, ncol = years)
    
    # insert s.yoy into projection matrices
    pmx_list <- lapply(1:numReps, function (n){ lapply(1:years, function(y) {
      pmxs[[n]][[y]][2,1] = s0.rand[n, y]; pmxs[[n]][[y]]
    })})
    
    # est. geom_mean of lambda of pmx.list
    lambdas <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$longterm) 
    #lambdas <- sapply(pmx_list, lambda_mean)
    mean_lambda <- mean(lambdas)
    
    # est distance from trarget lambda
    (mean_lambda - target.lambda)^2
  }
  
  no_cores <- detectCores() - 1
  
  cl <- makeCluster(no_cores) 
  clusterExport(cl, ls())
  clusterEvalQ(cl, library(popbio))
  
  M_optim <- parLapply(cl, 1:6, function(x){
    #cat("Location", x, "Of 6", "\n")
    
    years = 100
    numReps = 1000
    
    syoys <- lapply(lambdas, function(l){
      # Generate random variables
      pmx_rand <- pmx_build_runs(mx = A.LS,     # Projection matrix expression
                                 numReps = numReps,         # number of simulation replicates
                                 years = years,             # number of years to simulate
                                 lambda = l,                # Pop growth rate
                                 cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                                 names.Li = names_Li[[x]],  # names for growth increments
                                 dist.Li = dist_Li[[x]],    # distributions names for random length increments
                                 params.Li = params_Li[[x]],# distributions parameters fore random length increments
                                 cor.s = cor_s[[x]],        # Correlation struction for vital rates
                                 names.s = names_s[[x]],    # names for vital rates
                                 dist.s = dist_s[[x]],      # distributions names for random  vital rates
                                 params.s = params_s[[x]],  # distributions parameters fore random vital rate
                                 L.mean = L_mean[[x]],      # Initial (mean) length at age
                                 data = c(lh_data[[x]], l_min$equil[[x]], s0 = NA)) # Life history data 
      
      pmx_init <- pmx_rand$pmx
      
      rm(pmx_rand)
      gc()
      
      # Solve for syoy with lambda = 1
      syoy <- optimize(M_optim_f, c(0, 1e-1), tol = 1e-16, target.lambda = l,
                        numReps = numReps,         # number of simulation replicates
                        years = years,             # number of years to simulate
                        pmxs = pmx_init,           # random projection matrcies
                        data = c(lh_data[[x]], l_min$equil[[x]])) # Life history data 
      
      syoy$minimum
    })
    names(syoys) <- lambdas
    syoys
  })
  
  stopCluster(cl)
  
  rm(years, numReps, M_optim, M_optim_f)
}

# YOY Survival rates
syoys <- list(list(0.0001873, 0.0004585, 0.0009981, 0.0020712, 0.0040046, 0.0074005, 0.0128254, 0.0220947),
              list(0.0000578, 0.0001420, 0.0003214, 0.0008202, 0.0013904, 0.0026545, 0.0048307, 0.0084056),
              list(0.0000098, 0.0000285, 0.0000802, 0.0001852, 0.0004072, 0.0008550, 0.0016652, 0.0031087),
              list(0.0000247, 0.0000719, 0.0001942, 0.0004471, 0.0009506, 0.0019871, 0.0039020, 0.0073332),
              list(0.0000246, 0.0000630, 0.0001482, 0.0003292, 0.0006834, 0.0013205, 0.0024334, 0.0044260),
              list(0.0000629, 0.0001652, 0.0003625, 0.0008066, 0.0015601, 0.0029896, 0.0055043, 0.0095887)
)
syoys <- lapply(syoys, function(x){ names(x) = lambdas; x})
names(syoys) <- DU

# Stochastic vital rate parameters
params_s.list <- lapply(1:6, function(x) { # Survival
  lapply(1:length(lambdas), function(l){
  
    # Mean annual Survival 
    s_means <- surv_est(0:(lh_data[[x]]$Tmax - 1), c(lh_data[[x]], s0 = syoys[[x]][[l]]))
    
    # Mean annual Instantanesou Mortality
    M_means <- {
      Ms <- -log(s_means)
      names(Ms) <- paste0("M", 0:(length(s_means) - 1))
      Ms
    }
    
    lapply(M_means, function(i) beta_stretch_val(i, sd = (i * lh_data[[x]]$cv.M)))
  })
})

# Test optimization
if(FALSE){
  # simulation criteria
  numReps = 500
  years = 100
  x = 5 # DU Location for random variability plotting
  l = 8
  
  # Est. random matrices for each cor str
  system.time({
    pmx_rand <- pmx_build_runs(mx = A.LS,       # Projection matrix expression
                               numReps = numReps,         # number of simulation replicates
                               years = years,             # number of years to simulate
                               lambda = lambdas[[l]], # pop growth rate                                                ****
                               cor.Li = cor_Li[[x]],      # Correlation struction for growth increments
                               names.Li = names_Li[[x]],  # names for growth increments
                               dist.Li = dist_Li[[x]],    # distributions names for random length increments
                               params.Li = params_Li[[x]],# distributions parameters fore random length increments
                               cor.s = cor_s[[x]],        # Correlation struction for vital rates
                               names.s = names_s[[x]],    # names for vital rates
                               dist.s = dist_s[[x]],      # distributions names for random  vital rates
                               params.s = params_s.list[[x]][[l]],# distributions parameters fore random vital rate             ****
                               L.mean = L_mean[[x]],      # Initial (mean) length at age
                               data = c(lh_data[[x]], l_min$equil[[x]], s0 = syoys[[x]][[l]]), # Life history data  ****
                               sensitivity = FALSE)
  })
  
  # extract results  
  pmx_list <- pmx_rand$pmx  
  pmx.run <- pmx_create(c(lh_data[[x]], s0 = syoys[[x]][[l]], l_min$equil[[x]]))
  rm(pmx_rand) 
  
  # extract pop growth
  lambda_year <- (lapply(1:length(pmx_list), function(i) lambda.real(N0 = initPop(pmx.run, 100, "adult", c(0, 0, 0, 1, 1)), pmx_list[[i]])$annual))
  lambda_rep <- sapply(1:length(pmx_list), function(i) lambda.real(N0 = initPop(pmx.run, 100, "adult", c(0, 0, 0, 1, 1)), pmx_list[[i]])$longterm)
  
  mean(lambda_rep)   # mean lambda - overall
  lambdas[[l]]
}

#-------------------------------------------------------------------------------
# RECOVERY TIME SIMULATIONS
#-------------------------------------------------------------------------------

# Simulation runs
maxRuns <- 5      # number of independent runs
numReps <- 100    # number of replicated simulations
years <- 1500      # years to run simulation

# MVP parameters
p.cat <- c(0.15)  # prob of catastrophe
e.cat <- 0.5      # effect of catastrophe

# Mean adult pop sizes to loop through
adultPop <- c(50, 100, 250, 500, 1000)

# maturity
p.rep <- c(rep(0, 3), rep(1, 2))

# Run simulations
# Simulate annual population size for 'year' at each 'adultPop' for each of 'lambda'
# Each run is replicated 'numReps' times which is rum 'maxRuns' times

no_cores <- detectCores()  - 1

cl <- makeCluster(no_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(popbio))
# DU
rec_run <- parLapply(cl, 1:6, function(x) {
  
  # Adult population size
  lapply(as.list(adultPop), function(y) {
    
    # Lambdas
    lapply(1:length(lambdas), function(l){
      
      # matrix for current level of pop growth
      pmx.run <- pmx_create(c(lh_data[[x]], s0 = syoys[[x]][[l]], l_min$equil[[x]]))
      
      # Initial stable stage - skewed towards juveniles
      # based on Goulais bay data from Tom Pratt - not included on github
      init.SS <- {
        NJ1 <- 1052# length(which(GB_data$AGE < lh_data[[1]]$j1.Tmax))                                       # Number for pop in J1 stage
        NJ2 <- 167 # length(which(GB_data$AGE >= lh_data[[1]]$j1.Tmax & GB_data$AGE < lh_data[[1]]$j2.Tmax)) # Number for pop in J2 stage
        NA1 <- 11  # length(which(GB_data$AGE >= lh_data[[1]]$j2.Tmax & GB_data$AGE < lh_data[[1]]$a1.Tmax)) # Number for pop in A1 stage
        NA2 <- 0   # length(which(GB_data$AGE >= lh_data[[1]]$a1.Tmax))                                      # Number for pop in A2 stage
        
        # number of YOY need to produce above pop size
        NYOY <- {
          SS <- stable.stage(pmx.run)
          prop <- SS/sum(SS)
          sum(NJ1, NJ2, NA1, NA2)/sum(prop[2:5])
        }
        
        # create stabel stage dist. for staring pop
        c(NYOY, NJ1, NJ2, NA1, NA2)/sum(c(NYOY, NJ1, NJ2, NA1, NA2))
        
      }
      # Inital pop size
      N0 <- init.SS/sum(init.SS[3:5]) * y 
     
      rep_proj(mx = A.LS,                         # projection matrices as expression
               N0 = N0,                           # initial pop size
               maxRuns = maxRuns,                 # number of runs of simlulation (new matrices est each run) 
               years = years,                     # years to run simulation
               numReps = numReps,                 # Number of simulation per run
               lambda = lambdas[[l]],             # Pop growth rate
               cor.Li = cor_Li[[x]],              # Correlation struction for growth increments
               names.Li = names_Li[[x]],          # names for growth increments
               dist.Li = dist_Li[[x]],            # distributions names for random length increments
               params.Li = params_Li[[x]],        # distributions parameters fore random length increments
               cor.s = cor_s[[x]],                # Correlation struction for vital rates
               names.s = names_s[[x]],            # names for vital rates
               dist.s = dist_s[[x]],              # distributions names for random  vital rates
               params.s = params_s.list[[x]][[l]],# distributions parameters fore random vital rate
               L.mean = L_mean[[x]],              # Initial (mean) length at age  
               p.rep = p.rep,                     # Proportion mature by age class
               p.cat = p.cat,                     # Probability of castastrophe
               e.cat = e.cat,                     # effect on castastrophe on pop size (percent recuction)
               gen.time = generation_time[[x]],   # Generation Time
               data = c(lh_data[[x]], l_min$equil[[x]], s0 = syoys[[x]][[l]])# Life histoty data
      )
    }) # Lambdas
  }) # Addult Pop size
}) # DU

stopCluster(cl)

save(rec_run, file = "Results 06 - recovery population.Rdata")

# Clean up
rm(MVP_run, p.rep, adultPop, p.cat, e.cat, years, maxRuns, numReps)

#-------------------------------------------------------------------------------
# Extract Results
#-------------------------------------------------------------------------------

# Load MVP estimates
load("Results 03a - MVP estimates.Rdata")

# Load recovery simulation data
load("Results 06 - recovery population.Rdata")

# Simulation infor
lambdas <- seq(1.01, 1.15, 0.02) # pop growth rates simulated
maxRuns <- 5                     # number of independent runs
numReps <- 100                   # number of replicated simulations
years <- 1500                    # years to run simulation
adultPop <- c(50, 100, 250, 500, 1000) # inital pop sizes

# FUNCTIONS
# Time to recovery
rec_time <- function(x, threshold){
  ext <- which(x >= threshold)
  if(length(ext) < 1){
    ext = length(x)
  }
  min(ext)
}
# rec_time applied to a list
rec_time_list <- function(x, threshold){
  sapply(x, rec_time, threshold = threshold)
}

# Time to Recovery Data
# Find the year when recovery occurs for each simulation
RT_data <- lapply(1:6, function(x){
  A <- lapply(1:length(adultPop), function(a){
    L <- lapply(1:length(lambdas), function(l){
      
      as.data.frame(list(
        N = adultPop[[a]],         # inital pop size
        lambda = lambdas[[l]],     # Pop growth rate
        Time = rec_time_list(unlist(rec_run[[x]][[a]][[l]], recursive = F), threshold = MVP_est.A[[x]][[2]]), # Year Recovery occurs
        DU = DU[[x]],              # DU   
        cat = "0.15/Generation"    # Catastrophe rate
      ))
    })
    do.call(rbind, L)
  })
  do.call(rbind, A)
})

rm(rec_run)

# 95% probabilty of recovery
# Find the 95th percentile of recovery time fore each lambda within inital pop size
RT_95 <- lapply(1:6, function(x){
  A <- lapply(adultPop, function(a){
    L <- lapply(lambdas, function(l){
      
      data <- RT_data[[x]][RT_data[[x]]$N == a & RT_data[[x]]$lambda == l, ]
      
      as.data.frame(list(
        N = a,                           # inital pop size
        lambda = l,                      # pop growth rate
        Time = quantile(data$Time, 0.95),# 95th percentile of recovery times
        DU = DU[[x]],                    # DU   
        cat = "0.15/Generation"          # Catastrophe rate
      ))
    })
    do.call(rbind, L)
  })
  do.call(rbind, A)
})

#-------------------------------------------------------------------------------
# Model fit
# fit RT_95 data (95% change of recovery) with Time ~ Lambda
# To give a smooth line for plotting
# Fit one line for each pop size within DU
# Function form: y = a/(b + x)

# Model fits
rec_mod <- lapply(RT_95, function(x){
  lapply(adultPop, function(n){
    nls(Time ~ a /(lambda + b), data = x[x$N == n,], start = list(a = 100, b = 1))
  })
})

# Predicted data for plotting
xv <- seq(1.01, 1.15, 0.0001)
pred_data <- lapply(1:6, function(x){
  R <- lapply(1:length(adultPop), function(n){
    as.data.frame(list(
      N = adultPop[n],
      lambda = xv,
      Time = predict(rec_mod[[x]][[n]], list(lambda = xv)),
      DU = DU[[x]]
    ))
  })
  do.call(rbind, R)
})

# clean up
rm(rec_mod, xv)

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------
# 2 plots
# 1: density plots of time to recovery (RT_data) for each lambda for one pop size
#    and one DU. As an example of the distribution of recovery times
# 2: Line graph of time to 95% recovery probabilty for each pop size with each DU


# Density Plots

x = 1   # Location DU
n = 500 # Pop size

# 95% recovery data
# to demonstratte where data in next graph come from - plot as rug plot
rug_data <- as.data.frame(list(
  Time = RT_95[[x]]$Time[RT_95[[x]]$N == n],
  lambda = as.factor(lambdas)
))

# 1
# plot - Density plot of time to recovery accross lmabdas
p <- ggplot(RT_data[[x]][RT_data[[x]]$N == n,], aes(x = Time)) 
p <- p + geom_density(aes(fill = as.factor(lambda)), alpha = 0.75)
p <- p + geom_rug(data = rug_data, aes(x = Time, colour = lambda), size = 1, alpha = 0.75)
p <- p + scale_x_continuous(trans = "log10", breaks = c(10, 100, 1000), 
                            minor_breaks = c(seq(20, 90, 10), seq(200, 900, 100), 2000))
p <- p + guides(colour = FALSE)
p <- p + labs(x = "Time to Recovery (years)", y = "Density", fill = "Population\nGrowth Rate")           
p <- p + theme_me + 
  theme(axis.text.y = element_blank())

#png("Figures/Recovery time distribution.png", width = 6, height = 3, unit = "in", res = 300)
p
#dev.off()

# 2
# Plot - Lambda Vs Time to 95% likelihood of recovery
# Line for each pop size; panel for each DU
plot_data <- do.call(rbind, RT_95)
p <- ggplot(plot_data , aes(lambda, Time))
#p <- p + geom_point(aes(colour = as.factor(N)))
p <- p + geom_path(data = do.call(rbind, pred_data), aes(x = lambda, y = Time, colour = as.factor(N)))
p <- p + facet_wrap(~ as.factor(DU), ncol = 3)
p <- p + scale_y_continuous(limits = c(0, 300))
p <- p + scale_x_continuous(breaks = seq(1.01, 1.16, 0.03))
p <- p + labs(x = expression(bold(paste("Population Growth Rate ", (lambda)))), 
              y = "Time to Recovery (years)", colour = "Initial\nPopuation") 
p <- p + theme_me +
  theme(strip.text = element_text(size = 11, family = "sans", face = "bold"),
        axis.text.y = element_text(angle = 0))

png("Figures/Recovery time.png", width = 8, height = 4.5, unit = "in", res = 300)
p
dev.off()

# Clean up
rm(x, n, rug_data, pred_data, plot_data, RT_95, RT_data, naxRuns, numReps, years,
   lambdas)

################################################################################