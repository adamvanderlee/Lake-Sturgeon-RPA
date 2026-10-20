################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# PARAMETER FUNCTION
# functions used in estimation of parameters for warmouth model
#-------------------------------------------------------------------------------
require(nlstools)

#-------------------------------------------------------------------------------
# Growth functions

# Function to est growth based on VB relationship
# can use form with t0 or form with L0
# use of either form based on specification of t0 or L0
# can only specify 1 and other must be NA
VB.growth <- function(t, Linf, K, t0 = NA, L0 = NA, names = FALSE) {
     
    if(is.na(L0[1])){
        L <- Linf * (1 - exp(-K * (t - t0)))  # t0 form
    } else if (is.na(t0[1])) {
        L <- Linf - (Linf - L0) * exp(-K * t) # L0 form
    }
    if(names == TRUE) names(L) <- paste0("L", t)
    L
}

# function to est age from VB growth function
# can use form with t0 or form with L0
# use of either form based on specification of t0 or L0
# can only specify 1 and other must be NA
agebyVB <- function(L, Linf, K, t0 = NA, L0 = NA){
    
    if(is.na(L0[1])){
        t <- -1 / K * log( 1 - L / Linf) + t0    # t0 form
    } else if (is.na(t0[1])) {
        t <- log(-((L - Linf) / (Linf - L0)))/-K # L0 form
    }
    t
}

# Function to assign ages to length values
# Assuming: Ai when (Li + Li-1)/2 < Li < (Li + Li+1)/2
assign_age <- function(L, Tmax, Linf, K, t0 = NA, L0 = NA){
  
  # Calc age-specific ages
  Li <- VB.growth(0:Tmax, Linf, K, t0, L0)
  
  # Find the means of ages (Li + Li+1)/2
  Li_means <- 0
  for(i in 2:length(Li)){
    Li_means[i] <- (Li[i-1] + Li[i])/2
  }
  
  # Fucntion to assign ages to each length
  # Find the max Li_mean that is < Li
  assign_func <- function(x) max(which(Li_means < x) - 1)
  
  # Assign ages
  Age <- sapply(as.list(L), assign_func)
  
  # Output
  Age
}

# Function to calc length at age using growth increments
# Takes vector of current lenghts and adds increments while
# maintaining hatch size and Tmax. 
Li_growth <- function(Lt, Li){
  
  # Add annual growth increments to current length
  # Excluded length at Tmax and reincorporates L0
  L <- c(Lt[1], Lt[1:length(Li)] + Li)  
  names(L) <- paste0("L", 0:(length(Lt) - 1))
  # Output
  L
}

# Function to generate random inital length-at-age from random growth incrementts
L_init.fun <- function(n, L_mean, cor, names, dist, params){
  
  # Random grwoth increments
  Li_rand <- rdist_cor(n, cor, names, dist, params, method = "mvrnorm")
 
  # Sum growth increments to length at age iteratively
   Lt <- list(L_mean)  # Initial length
  for(i in 2:dim(Li_rand)[1]){
    Lt[i] <- list(Li_growth(Lt[[i-1]], Li_rand [i, ]))
  }
   
  # Output final length
  Lt[[n]]
}


# Functino to calc weight at age 
# Takes vector of lengths and convets to weights
# Agges names if requested
LW_fun <- function(L, params, names = TRUE){
  W <- params$a.LW * L ^ params$b.LW; 
  names(W) = NULL
  if(names == TRUE) names(W) <- paste0("W", 0:(length(L) - 1))
  W
}


# Bootstrap confidence intervals for VBGF
VBbootCI <- function(fit,          # VBGF nls fit object
                     Age,          # Ages
                     t0 = NA,      # range of t0 values if forced fit
                     Linf = NA,    # range of Linf values if forced fit
                     alpha = 0.05, # level of confidence
                     niter = 1000) # number of iterations
    {
    # Bootstrap est. of coef values for VB fit
    boot <- nlsBoot(fit, niter = niter) # execute bootstrap est
    ests <- data.frame(boot$coefboot)   # extract coef est
    niter <- dim(ests)[1]               # set niter to number of successful coef est.
    
    # If t0 was force assige coef value by random sample of uniform dist. from value range
    if(is.numeric(t0) == TRUE){ 
        ests$t0 <- runif(niter, min = min(t0), max = max(t0))
    }
    
    # If Linf was force assige coef value by random sample of uniform dist. from value range
    if(is.numeric(Linf) == TRUE) { 
        ests$Linf <- runif(niter, min = min(Linf), max = max(Linf))
    }
    
    # Est age for each bootstrap est and extract CI by age
    pv <- matrix(NA, niter, length(Age)) # initialize outpur matrix for length est.
    CI <- matrix(NA, 2, length(Age))     # initalize output matrix for CI est.
    
    # Loop through ages
    for(i in 1:length(Age)){
        pv[,i] <- VB.growth(t = Age[i], Linf = ests[,"Linf"], K = ests[,"K"], t0 = ests[,"t0"])
        qtile <- quantile(pv[,i], c(alpha/2, 1 - alpha/2))
        CI[1,i] <- qtile[1]
        CI[2,i] <- qtile[2]
    }
    # out put
    out <- list(CI = CI, bootA = pv)
}


#### Linear regression confidence intervals (FECUNDITY)
fCI <- function(fitF, Lboot, meanL){  #fitlm=loglinear fit F~L
    meanF <- predict(fitF, list(TL.mm = meanL))
    CIup <- Lboot$CI[2,]
    CIlo <- Lboot$CI[1,]
    fhi <- predict(fitF, list(TL.mm = CIup), interval="confidence")[,"upr"]
    flo <- predict(fitF, list(TL.mm = CIlo), interval="confidence")[,"lwr"]
    fsd <- (fhi-flo)/(2 * 1.96)
    (fhi - meanF) / 1.96
    (meanF - flo) / 1.96
    out <- list(f.mean=exp(meanF), f.sd=fsd, fhi=fhi, flo=flo)
}


#-------------------------------------------------------------------------------
# Survival functions

# Fucntion to est Beta - weighted catch curve for est of m0
W_catchCurve_m0 <- function(L, Age = NA, K, Linf, t0, H0, plot = FALSE) {
    
    # Histogram - to est frequency per length bin [corrected if von bert]
    if(is.numeric(Age) == FALSE){
        Age <- agebyVB(L, Linf, K, t0 = t0)
    }
    t = floor(min(Age, na.rm=TRUE)):ceiling(max(Age, na.rm=TRUE))
    if(floor(min(Age, na.rm=TRUE)) < 0) t = 0:ceiling(max(Age, na.rm=TRUE))
    Lt = VB.growth(t = t, Linf = Linf, K = K, t0 = t0)
    if(min(t) == 0) Lt[1] = H0
    breaks = log(Lt) + K * t
    histo <- hist(log(L) + K * Age, breaks = breaks, plot = FALSE)

    bin.start <- which(histo$counts == max(histo$counts)) # Indentify when count is max
    bin.end <- length(histo$counts)                       # identify when count ends
    
    # Y values for catch curve
    freq <- histo$counts[bin.start:bin.end]               # Counts
    use.freq <- freq[freq>0]                              # Removes 0 obs from counts
    
    # x values for catch curve
    xval <- histo$mids[bin.start:bin.end]                 # x values - mid points of bins in histo
    use.xval <- xval[freq>0]                              # remove xvals with 0 counts
    
    # Initial catch curve fit
    mod <- lm(log(use.freq)~use.xval)                     # log(N) ~ L [+K*age]
    
    # predictio from catch curve to use in weight mod
    expected.freq <- predict(mod)                         # extrect preds
    expected.freq[expected.freq < 0] <- 0                 # set any -ves to 0
    
    # Weighted catch curve fit - reduces impact of greater than expected
    # freq at larger sizes
    weighted.mod <- lm(log(use.freq)~use.xval, weights=expected.freq)
    
    # plot catch curve
    if(plot == TRUE){
        plot(use.xval, log(use.freq), pch=19, type="h", lwd=2)
        abline(mod)
        abline(weighted.mod, col="red")
    }
       
    # Extract Beta values from weighted catch curve and SE
    Beta <- coef(weighted.mod)[[2]]
    sdBeta <- summary(weighted.mod)$coef[2,2]

    # Output 
    out <- list(Beta=Beta, sdbeta=sdBeta, betamod=weighted.mod, N=sum(use.freq))
}

# Estimate survival between consecutive lengths using Lorrenzons
# formula and Koops est. of m0
# si = [Li*exp(-k*t)/Li+1]^(m0/linf*k)
s.m0.VB <- function(L, Beta, sdbeta, Linf, K){
    # Est m0
    m0 <- -K * Linf * Beta
    
    # Est Survival
    N <- length(L)                         # Number of L obs
    L.change <- L[1:N-1] / L[2:N] * exp(-K)# Change in Lengths between sequential obs
    si <- L.change^(m0/K/Linf)             # survival between lengths
    
    # Std Dev of survival obs
    sdm0 <- K * Linf * sdbeta              # SD of m0
    sd.si <- abs(si * log(L.change) * (1 / K / Linf)) * sdm0 # SD of Si
    sd.logsi <- 1/si * sd.si               # SD of log(si) ???
    
    # output
    out <- list(m0 = m0, si=si, sd.si=sd.si, sd.logsi=sd.logsi)
}

# Size dependent mortality 
# Weighted catch curve from length frequencies (w_cc_LF)
# Function takes length-frequency data and estimates the m0 values used to predict
# annual survival rates by body size
# Uses a modified catch curve: log(count) ~ log(L) + K*t # (weighted by predicted count)
# The slope - Beta = m0 / (K * Linf)
# survival Si = (Li / Li+1 * exp(-K)) ^ (m0 / (K * Linf)
# Data required as length counts with predcited ages
w_cc_LF <- function(data, K, Linf) {
  
  # find mode of length counts and use >= mode fore fit
  L_mode <- which.max(data$Count)
  sub_data <- data[L_mode:dim(data)[1], ]
  
  # remove 0 counts
  sub_data <- sub_data[sub_data$Count > 0, ]
  
  # Fit model
  mod <- lm(log(Count) ~ I(log(Length) + K * Age), sub_data)  
  
  # Predict counts
  expected.freq <- predict(mod)         # extrect preds
  expected.freq[expected.freq < 0] <- 0 # set any -ves to 0
  
  # Fit weighted model
  weighted.mod <- lm(log(Count) ~ I(log(Length) + K * Age), sub_data, weights = expected.freq) 
  
  # Extract Beta values from weighted catch curve 
  Beta <- coef(weighted.mod )[[2]]           # Slope 
  
  # Estmate m0 value
  m0 <- -K * Linf * Beta
  
  # Output
  out <- list(m0 = m0, model = weighted.mod, data = sub_data)
  
}

# Basic catch curve based on length at age data
catchCurve <- function(Age, Count = NA, plot = TRUE){
    
    # Calc age freq if not provided
    if(is.na(Count[1])){
        x <- hist(Age, plot = FALSE, 
             breaks = seq(floor(min(Age, na.rm = TRUE)), ceiling(max(Age, na.rm = TRUE)), 1))
        Count <- x$counts
        Age <- x$breaks
    }
    
    # Ident where catch declines
    start <- which(Count == max(Count))
    end <- length(Count)
    
    # y values for catch curve
    N <- Count[start:end]
    use.N <- N[N>0]
    
    # x values for catch curve
    A <- Age[start:end]
    use.A <- A[N>0]
    
    # fit model
    mod <- lm(log(use.N) ~ use.A)
    
    # prediction from catch curve to use in weight mod
    expected.freq <- predict(mod)                         # extrect preds
    expected.freq[expected.freq < 0] <- 0                 # set any -ves to 0
    
    # Weighted catch curve fit - reduces impact of greater than expected
    # freq at larger sizes
    weighted.mod <- lm(log(use.N)~use.A, weights=expected.freq)
    
    # Plot
    if(plot == TRUE){
        plot(use.A, log(use.N), ylab = "log(Catch)", xlab = "Age", pch=19, type="h", lwd=2)
        abline(mod)                     # inital catch curve line
        abline(weighted.mod, col="red") # weighted catch curve
    }
    
    # Extract Beta values from weighted catch curve and SE
    Intercept <- coef(weighted.mod )[[1]]
    Beta <- coef(weighted.mod )[[2]]
    sdBeta <- summary(weighted.mod )$coef[2,2]
    
    out <- list(Intercept = Intercept, Beta = Beta, sdbeta = sdBeta, 
                mod = mod, N = sum(use.N), Ages = A)
    
    out
}



# Probabilitic age assignemnt from length
P_ageByL <- function(Length, meanLength, ages, se, reps){
    
    probMat <- matrix(NA, nrow = length(Length), ncol = length(meanLength))
    for(i in 1:length(meanLength)){
        
        probMat[Length <= meanLength[i], i] = pnorm(Length[Length <= meanLength[i]],
                                                    meanLength[i],VBse.L0)
        probMat[Length > meanLength[i], i] = pnorm(Length[Length > meanLength[i]],
                                                   meanLength[i],VBse.L0, lower.tail = FALSE)
        
    }
    
    probMat <- probMat / apply(probMat,1,sum)
    
    age <- matrix(NA, nrow = dim(probMat)[1], ncol = reps)
    for(j in 1:reps) {
        for(i in 1:dim(probMat)[1]){
            age[i,j] <- sample(ages, 1, prob = probMat[i, ])
        } 
    }
    
    age
}

# Est mortality using weight catch curve 
# Output: mean M and sd of M
P_catchCurve <- function(ageMat){
    
    Ms <- NULL
    for(i in 1:dim(ageMat)[2]){
        Ms <- c(Ms,catchCurve(ageMat[,i])$Beta, plot = FALSE)
    }
    
    as.data.frame(list(mean = median(abs(Ms)), sd = sd(Ms)))
    
}


#-------------------------------------------------------------------------------
# Estimate population growth rate and predict confidence/prediction interval
# Based on data and fit from Randall et al (1995) and adjustemnt to relationship
# from Randall and Minns (2000)

r_max_pred <- function(w,         # prediction weight
                       h,         # prediction habitat (1 = lake; 2 = river)
                       interval){ # confidence or prediction interval
    
    # Data from Randall et al. (1995)
    data <- as.data.frame(list(
        Production = c(30.3, 3.8, 2.0, 228.5, 15.5, 4.7, 21.5, 30.8, 5.7, 85.7, 8.6,
                       270.0, 596.0, 198.0, 140.0, 100.0,  74.0,  52.0, 30.0, 174.5,
                       280.3, 243.5, 226.4, 375.4, 379.1, 217.6, 500.6, 297.2, 110.3,
                       128.8, 161.4, 83.4, 807.2, 228.0, 329.0, 168.7, 57.2, 147.3,
                       54.8, 76.5, 30.7, 35.3, 261.5, 98.3, 39.3, 33.3, 47.5, 26.1,
                       31.2, 53.0, 456.2, 80.3, 29.4),
        biomass = c(73.0, 7.3, 9.2, 372.2, 32.2, 11.2, 28.9, 11.6, 24.4, 277.8, 29.5,
                    161.2, 198.0,  95.0,  75.0,  59.0,  51.0, 129.0, 62.0, 124.4,
                    274.8, 228.4, 149.9, 376.0, 278.7, 127.5, 217.5, 187.2, 111.2,
                    241.0, 307.5, 50.2, 310.5, 142.5, 86.6, 45.6, 10.8, 40.9, 49.8,
                    42.5, 38.4, 37.4, 173.0, 71.0, 30.5, 21.3, 38.5, 21.5, 27.1, 36.0,
                    233.4, 80.2, 15.4),
        weight = c(18.0246914, 3.5540409, 25.8426966, 6.6786291, 24.2105263, 37.8378378,
                   17.8947368, 5.7058534, 219.8198198, 13.6712598, 64.2701525, 1.24,
                   0.2544987, 1.4843750, 3.0, 2.4583333, 5.6666667, 18.4285714, 15.5,
                   2.4777915, 2.3869707, 1.4815487, 0.7197976, 3.7586469, 1.3582997,
                   2.0675564, 1.6603053, 0.6602313, 2.8333376, 13.9872316, 6.1428743,
                   2.1183222, 9.3047648, 4.2462529, 4.4194948, 0.4892704, 1.6410880,
                   0.8880686, 7.9935795,   0.9173916, 5.9222702, 10.0700054, 0.8397202,
                   1.3179144, 1.9223497, 2.1759117, 3.3804548, 2.6022755, 3.7215051,
                   2.7067669, 5.1344098, 12.0819524, 0.6273679),
        habitat = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1)
    ))
    
    # fit regression model
    mod = lm(log10(Production/biomass) ~ log10(weight) + habitat, data = data)
    
    # x prediction values
    pred.data = data.frame(weight = w, habitat = h)
    
    # make prediction with desired interval estiamte
    rmax = (predict(mod, pred.data, interval = interval))
    
    # retransform and multiply by 2 (Randall and Minns 2000)
    10^rmax*2
}


#-------------------------------------------------------------------------------
# OTHER

# function to calc geometric mean between sequenital values in a vector
geom_mean_seq <- function(vec) {
    n = length(vec)
    sqrt(vec[2:n]*vec[1:(n-1)])
}

# geometric mean of a vector
geom_mean <- function(vec, na.rm = FALSE) {
    exp(mean(log(vec), na.rm = na.rm))
}


