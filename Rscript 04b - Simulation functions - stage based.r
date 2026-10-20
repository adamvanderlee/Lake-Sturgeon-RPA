################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# Lake Sturgeon
# simulation functions
#===============================================================================
# Distribution functions

# Function to estimate shape parameters of the beta distribution
# from the defined mean and variance of data
beta_val <- function(mean, var){
  
  a <- mean * ((mean * (1 - mean)) / var - 1)
  b <- (1 - mean) * ((mean * (1 - mean)) / var - 1)
  
  list(a = a, b = b)
}

# Function to estimate shape parameters of the streched beta distribution
# from the defined mean and variance of data
# Streched beta is a beta dist'n scaled to outside (0,1) limits
# shape parameters give regulare beat dist'n which must be scaled with
# (max - min) + min when generating values with rbeta
beta_stretch_val <- function(mean, var = NA, sd = NA){
  
  if(is.na(var[1])) var = sd^2    # calc variance if not given
  if(is.na(sd[1])) sd = sqrt(var) # calc sd if not given

  # estimate mins and max of distribution
  min = qnorm(1e-6, mean, sd = sd)
  max = qnorm(1 - 1e-6, mean, sd = sd)
  
  # Scale means and variance to between 0 and 1
  mean.b <- (mean - min) / (max - min)
  var.b <- var * (1 / (max - min))^2
  
  # estimate beta dist'n shape parameters
  a <- mean.b * ((mean.b * (1 - mean.b)) / var.b - 1)
  b <- (1 - mean.b) * ((mean.b * (1 - mean.b)) / var.b - 1)
  
  # output
  list(a = a, b = b, min = min, max = max)
}

# Function to estimate the parameters of a log-normal distribution
# given the log(mean) and log(sd) or log(var) used in rlnorm function
# [required for accurate converstion between lnorm and beta dist'n]
lnorm_var <- function(meanlog, sdlog = NA, varlog = NA){
  if(is.na(varlog)) varlog <- sdlog^2
  
  mean <- exp(meanlog + 1/2 * varlog)
  var <- exp(2 * meanlog + varlog) * (exp(varlog) - 1)
  
  list(mean = mean, var = var)
}

# Function to generate an ARI correlation structure
AR1 <- function(n, rho, names = NA){
  C <- diag(n)
  C <- rho ^ abs(row(C) - col(C))
  
  if(is.character(names)) colnames(C) <- names
  
  C
}

# Function to build correlation matrix for intra-annual variability among vital rates
# Assumes n-1 survival terms and 1 fecundity term. Survival rates rate correlated with
# each other following an AR1 structure and a values of rho.s. We assume YOY survival
# is correlated with the fecudnity term with a value of rho.f
cor_mat <- function(n, rho.s, rho.f, names = NA) {
  
  C <- diag(n)                                  # Create a idenity matrix of dimenions n*n
  C[1:(n-1), 1:(n-1)] <- AR1((n-1), rho.s)      # Assiing an AR1 structure among survival terms
  C[1, n] <- C[n, 1] <- rho.f                   # Assign correlation between fecundity and YOY survival
  
  if(is.character(names)) colnames(C) <- names  # Assign names to matrix if provided
   
  C                                             # OUtput
}

# Function to generate a set of random standard normal values (mean = 0, sd = 1)
# that have a specified correlation between 
# Based on code from: 
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable/15035#15035
rnorm_cor <- function(n, cor, names = NA) {
  
  # number of variables
  vr_num <- sqrt(length(cor))
  
  # Choleski Decomposition of correlation matrix
  CC <- chol(cor)
  
  X <- replicate(vr_num, rnorm(n, mean = 0, sd = 1))
  
  # induce correlation (does not change X1)
  df <- X %*% CC
  
  if(is.character(names)) colnames(df) <- names
  
  df
  
}

# Function that estimates correlated random values with various distriubtional shapes
# Applies normal cummulative distribution function to each variable and then the 
# inverse CDF of desired distribution
rdist_cor <- function(n, cor, names = NA, dist, params, method = "Choleski") {
  
  # number of variables
  vr_num <- sqrt(length(cor))
  
  # Generate correlated random standard values
  if(method == "Choleski"){
    X <- rnorm_cor(n, cor, names)
  } else if( method == "mvrnorm") {
    require(MASS)
    X <- mvrnorm(n, mu = rep(0, vr_num), Sigma = cor)
  } else {
    warning("method not recognized")
  }
  
  # Convert to probabilities
  pX <- pnorm(X, mean = 0, sd = 1)
  
  # USe probability vectors to prodice random variabiles with alt dist'n
  qX <- matrix(NA, nrow = n, ncol = vr_num)
  
  for(i in 1:dim(pX)[2]){
    qX[,i] <- rdist_change(pX[,i], dist[i], params[[i]])
  }
  
  if(is.character(names)) colnames(qX) <- names
  
  qX
  
}

# Function that takes probability fucntion from a random standard normal vector
# and converts to random values of alternative dist'n
rdist_change <- function(x, dist, params) {
  
  if(dist == "norm"){
    qX <- qnorm(x, mean = params$mean, sd = params$sd)
  } else if(dist == "lnorm") {
    qX <- qlnorm(x, meanlog = params$meanlog, sd = params$sdlog)
  } else if(dist == "beta") {
    qX <- qbeta(x, shape1 = params$a, shape2 = params$b)
  } else if(dist == "beta_stretch"){
    qX <- qbeta(x, shape1 = params$a, shape2 = params$b) * (params$max - params$min) + params$min
  } else if(dist == "uniform") {
    qX <- qunif(x, min = params$min, max = params$max)
  }
  
  qX
}


# Function that creates a correlation matrix for simulations with within and between
# year correlations among variables. Must start with two n_var*n_var dimension
# Assumes constant within year correlation over time and an AR1 structure among years
# correaltion matrices
# See Morris and Doak (2002) pg 295
btwYear_Cor <- function(n_var, # Number of variables
                        years, # Number of years to simluation
                        C,     # Within year correlation matrix among varialbes
                        B      # Between year correlation matrix among variables
){
  
  # initalize correlation matrix with dimention (n_vars*years)^2
  X <- matrix(0, ncol = n_vars*years, nrow = n_vars*years)
  
  # Assign sub-matrices to blocks of X
  for(i in 1:years){   # Loop through years
    for(j in 1:years){ # loop through years
      
      if(i == j){ 
        mx <- C    # When i = j (diagonals of X) used within year correaltion structure
      } else { 
        mx <- B    # Else use btween year correlation structure
        }
      
      # Lower left matrix strucutre
      # raise to the power of i - j
      if(i > j) {
        expo <- i - j
        mx <- mx^expo
      }
      
      # Upper right matrix structure
      # Transpose and raise to the power of j-i
      if(i < j) {
        expo <- j - i
        mx <- t(mx)^expo
      }
      
      # assign sub-matrices
      for(ii in 1:n_vars) {
        for(jj in 1:n_vars) {
          
          X[ii + n_vars*(i - 1), jj + n_vars * (j - 1)] <- mx[ii, jj]
          
        } # close jj
      }   # close ii
    }     # close j
  }       # close i
  
  # Output
  X
  
} # close function

#===============================================================================

# Build projection matrices vital rates of a single replicate as a list of
# length years
vr_gen <- function(Lis,      # Matrix of random length increments
                   Ms,       # MAtrix of random vital rates
                   L.init,   # Initial length at age
                   data)  # life history data (needs b.F)
{
  
  # Calculate length at age from growth increments (assuming inital mean size at age)
  Lt <- list(L.init)  # Initial length
  for(i in 2:dim(Lis)[1]){
    Lt[i] <- list(Li_growth(Lt[[i-1]], Lis[i, ]))
  }
  
  # Convert lengths to weight
  Wt <- lapply(Lt, function(x) LW_fun(x, data, names = TRUE))
  
  # Estimate age-specific fecundity
  #fs <- lapply(Wt, function(x) Ex_means(prefix = "f", suffix = c(10:34, "a"), data = c(data, x)))
  fs <- lapply(Wt, function(x) with(data, fec_est(x[-length(x)], a = a.F, b.F, names = TRUE)))

  # Extract annual survival rates - organize into lise
  ss <- lapply(unlist(apply(Ms, 1, list), recursive = FALSE), function(x) exp(-x))
  
  # output
  list(Lt = Lt, fs = fs, ss = ss)
  
}


# Function to estimate elasticity values for vit
# Build projection matrices of a numRep replicates as a list of length years
# Produced nmat random variables for growth length increments and vital rates
# requires random variation info for Li and vr includeing:
#    Correlation matrix, variable names, distributions, and parameter values
pmx_build_runs <- function(mx,       # Projection matrix expression
                           numReps,  # number of simulation replicates
                           years,    # number of years to simulate
                           lambda,   # pop growth rate
                           cor.Li,   # Correlation struction for growth increments
                           names.Li, # names for growth increments
                           dist.Li,  # distributions names for random length increments
                           params.Li,# distributions parameters fore random length increments
                           cor.s,    # Correlation struction for vital rates
                           names.s,  # names for vital rates
                           dist.s,   # distributions names for random  vital rates
                           params.s, # distributions parameters fore random vital rate
                           L.mean,   # Initial (mean) length at age
                           data,     # Life history data 
                           sensitivity = FALSE # Boolean: 1 - for sensitivity (dont genrate matricies); 0 - compile matrices
){
  
  # Sim list of initial age-specific length values (numReps long)
  L.init <- replicate(numReps, L_init.fun(n = 100, L_mean = L.mean, cor = cor.Li, names = names.Li, 
                       dist = dist.Li, params = params.Li), simplify = FALSE)

  # Calc number of parameter combos needed
  nmats = numReps * years
  
  # random growth increments
  Li_rand <- rdist_cor(n = nmats, cor = cor.Li, names = names.Li, 
                       dist = dist.Li, params = params.Li, method = "mvrnorm")
  
  # Random vital rates
  M_rand <- rdist_cor(n = nmats, cor = cor.s, names = names.s, 
                      dist = dist.s,  params = params.s, method = "mvrnorm")
  
  matNum <- 1:nmats # create vector of parameter combos IDs
  sample <- sample(matNum, nmats, replace = FALSE)
  sample <- lapply(1:numReps, function(x) matrix(sample, ncol = years, byrow = TRUE)[x,])
  
  # Create matrices for sampled parameters 
  vr_list <- lapply(1:numReps, function(x) {
    vr_gen(Lis = Li_rand[sample[[x]],], Ms = M_rand[sample[[x]],], L.init = L.init[[x]], data)
  })
  
  # clean up
  rm(Li_rand, M_rand)
  
  # extract vital rates
  Lts <- lapply(1:numReps, function(x) vr_list[[x]]$Lt)
  ss <- lapply(1:numReps, function(x) vr_list[[x]]$ss)
  fs <- lapply(1:numReps, function(x) vr_list[[x]]$fs)
  
  rm(vr_list)
  
  # Calc maturity ogive
  MO <-  mat_ogive(0:(data$Tmax-1), b = data$b.MO, c = data$c.MO, names = TRUE)

  # Calc age-specific fertility
  Fs <- lapply(1:numReps, function(x) {
    lapply(1:years, function(y){
      Fertility(data, fs[[x]][[y]], MO)
    })
  })

  # Aggregate survival to stage level
  s_stage <- lapply(1:numReps, function(x) lapply(1:years, function(y){
    with(data, list(
      syoy = ss[[x]][[y]][[1]],
      sj1 = exp(-mean(-log(ss[[x]][[y]][2:j1.Tmax]))),
      sj2 = exp(-mean(-log(ss[[x]][[y]][(j1.Tmax+1):j2.Tmax]))),
      sa1 = exp(-mean(-log(ss[[x]][[y]][(j2.Tmax+1):a1.Tmax]))),
      sa2 = exp(-mean(-log(ss[[x]][[y]][(a1.Tmax+1):(Tmax)])))
    ))
  }))
  
  # Aggregate Fertility to stage level - weighted mean using stage dist's
  F_stage <- lapply(1:numReps, function(x) lapply(1:years, function(y){
    
    # Age-specific mortality rate
    #M <- -log(ss[[x]][[y]])
    M <- -log(surv_est(0:(data$Tmax - 1), data))
    
    # stage stage dist'n within stage - use to weight fecundity est.
    wj <- exp(-cumsum(M[(data$j1.Tmax + 1):data$j2.Tmax]))
    w1 <- exp(-cumsum(M[(data$j2.Tmax + 1):data$a1.Tmax]))
    w2 <- exp(-cumsum(M[(data$a1.Tmax + 1):data$Tmax]))
    
    # aggregate Fertility to stage estimat
    list(
      Fj = sum(Fs[[x]][[y]][(data$j1.Tmax + 1):data$c.MO] * wj) / sum(wj), # Pre adult - first half of maturity ogive
      F1 = sum(Fs[[x]][[y]][(data$c.MO + 1):data$a1.Tmax] * w1) / sum(w1), # young adult - Tmat to half adult stage
      F2 = sum(Fs[[x]][[y]][(data$a1.Tmax+1):(data$Tmax)] * w2) / sum(w2)  # old adult - ent young adult to Tmax
      )
  }))
  
  # Calc stage-specific transition probabilities
  G_stage <- lapply(1:numReps, function(x) lapply(1:years, function(y){
    c(Ex_means(prefix = "Gj", suffix = 1:2, data = c(data, s_stage[[x]][[y]], lambda = lambda)),
      Ex_means(prefix = "Ga", suffix = 1:2, data = c(data, s_stage[[x]][[y]], lambda = lambda)))
  }))
    
  
  if(any(with(data, c(M0, Mj, Ma, Mf)) < 1)){
    if(data$M0 < 1) {
      
      s_stage <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          s_stage[[n]][[y]]$syoy = s_stage[[n]][[y]]$syoy * data$M0 
          s_stage[[n]][[y]]
          })
      })
      
    } else if(data$Mj < 1) { 
      s_stage <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          s_stage[[n]][[y]]$sj1 = s_stage[[n]][[y]]$sj1 * data$Mj 
          s_stage[[n]][[y]]$sj2 = s_stage[[n]][[y]]$sj2 * data$Mj 
          s_stage[[n]][[y]]
        })
      })
    } else if(data$Ma < 1) { 
      s_stage <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          s_stage[[n]][[y]]$sa1 = s_stage[[n]][[y]]$sa1 * data$Ma 
          s_stage[[n]][[y]]$sa2 = s_stage[[n]][[y]]$sa2 * data$Ma 
          s_stage[[n]][[y]]
        })
      })
    } else if(data$Mf < 1) {  # Fecundity reduction
      F_stage <- lapply(1:numReps, function (n){
        lapply(1:years, function(y) {
          F_stage[[n]][[y]]$Fj = F_stage[[n]][[y]]$Fj * data$Mf 
          F_stage[[n]][[y]]$F1 = F_stage[[n]][[y]]$F1 * data$Mf 
          F_stage[[n]][[y]]$F2 = F_stage[[n]][[y]]$F2 * data$Mf 
          F_stage[[n]][[y]]
        })
      })
    }
  }
  
  # Evaluate projection matrices for each year within replicate
  if(sensitivity){
    pmx <- NA
  } else {
    pmx <- lapply(1:numReps, function(x){
      lapply(1:years, function(y) pmx_eval(mx, 
                                           vars = c(data, s_stage[[x]][[y]], F_stage[[x]][[y]], G_stage[[x]][[y]])))
    })
  }

  # clean up
  rm(MO, sample, fs, ss)
  
  # output
  list(pmx = pmx, Lts = Lts, Fs = F_stage, ss = s_stage, Gs = G_stage)
}

#===============================================================================
# Elasticity functions

# function to calc elasticity values for vital rates
elasticity <- function (x, mx){
  require(popbio)
  sensitivity <- vitalsens(elements = mx, x) # calc vital rate sensitivity
  elast <- sensitivity$elasticity            # extract elasticities
  names(elast) <- rownames(sensitivity)      # assign names
  
  elast                                          # output
}

# Function to apply elasticity over a list and organize as a matrix 
# with each vital rate in a column
elasticity_list <- function (x, mx) { t(sapply(x, elasticity, mx)) }

# Function to randomly generate list of projection matrices and est. mean elasticity
# for vital rates of replicate simluatioin
elasticity_lambdas <- function(mx = A.LS,            # Projection matrix expression
                               numReps = 500,        # number of simulation replicates
                               years = 100,          # number of years to simulate
                               lambda,               # pop growth rate
                               cor.Li = cor_Li,      # Correlation struction for growth increments
                               names.Li = names_Li,  # names for growth increments
                               dist.Li = dist_Li,    # distributions names for random length increments
                               params.Li = params_Li,# distributions parameters fore random length increments
                               cor.s = cor_s,        # Correlation struction for vital rates
                               names.s = names_s,    # names for vital rates
                               dist.s = dist_s,      # distributions names for random  vital rates
                               params.s = params_s,  # distributions parameters fore random vital rate
                               L.mean = L_mean,      # Initial (mean) length at age
                               data = lh_data)       # Life history data 
{

  # gegnerate random vital rate
  pmx_rand <- pmx_build_runs(mx, numReps, years, lambda,
                             cor.Li, names.Li, dist.Li, params.Li,
                             cor.s, names.s, dist.s, params.s,
                             L.mean, data, sensitivity = TRUE)
  
  # extract list of vital rate values
  ss <- pmx_rand$ss
  Fs <- pmx_rand$Fs
  Gs <- pmx_rand$Gs
  
  rm(pmx_rand)
  
  vr_list <- lapply(1:numReps, function(x) lapply(1:years, function(y) {
    c(ss[[x]][[y]], Fs[[x]][[y]], Gs[[x]][[y]])
  }))
  rm(ss, Fs, Gs)
  
  # estimate elasticity of each vital rate (output a list of matrices of year specific elasticity)
  
  elast_list_year <- lapply(1:numReps, function(x) lapply(1:years, function(y){
    elasticity(vr_list[[x]][[y]], mx = mx)
  }))
  rm(vr_list)
  
  # Output
  elast_list_year
}



# Function to take the mean value of elasticities over years within replicates
elasticity_mean <- function(e){
  
  # take the mean of elasticity over years within a replicate
  # output is a list mean elasticity values for each vital rate for each replicate
  elast_list_mean <- lapply(e, function(x) apply(x, 2, mean))
  
  # find teh confidence intervals of elasticity over years within a replicate
  # output is a list mean elasticity values for each vital rate for each replicate
  elast_list_UCI <- lapply(e, function(x) apply(x, 2, quantile, probs = c(0.975)))
  elast_list_LCI <- lapply(e, function(x) apply(x, 2, quantile, probs = c(0.025)))
  
  # output - list of length numReps
  list(means = elast_list_mean, LCI = elast_list_LCI,  UCI = elast_list_UCI)
  
}


#===============================================================================

# Function to generate initial age specific population size base on the stable
# stage distn and a desired pop size

initPop <- function(pmx,         # projection matrix
                    pop,         # desired pop size
                    pop.stage,   # stage pop repsents - character ["total", "adult","juvenile","yoy","age1+"]
                    p.rep        # Proportion adult by age class              
            ) {
    
    # needs popbio library
    require(popbio)
    
    # Total population size
    if(pop.stage == "total") {
        N <- stable.stage(pmx) * pop                      # Calc pop based on stable stage dist
    }
    
    # Adult population size (Tmat to Tmax)
    if(pop.stage == "adult") {
        SS <- stable.stage(pmx)                           # Stable stage distn
        Na <- SS * p.rep / sum(SS * p.rep) * pop          # Calc adult pop
        Nj <- SS * (1-p.rep)/ sum(SS * p.rep) * pop       # Calc juv pop
        N <- Nj + Na                                      # Sum adult and juv pop by age
        
    }
    
    # juvenile popuation size (age 1 to Tmat)
    if(pop.stage == "juvenile") {
        SS <- stable.stage(pmx)                           # Stable stage distn
        p.juv <- 1 - p.rep                                # Calc proportion juvenile by stage
        p.juv[1] <- 0                                     # Assign 0 to YOY [YOY dont count as juvs]
        Nj <- SS * p.juv / sum(SS * p.juv) * pop          # Calc juv pop
        Na <- SS * (1 - p.juv)/ sum(SS * p.juv) * pop     # Calc non-juv pop
        N <- Nj + Na                                      # Sum juv and non-juv pop by age
    }

    # YOY population size (age 0 to 1)
    if(pop.stage == "yoy") {
        SS <- stable.stage(pmx)                           # Stable stage distn
        p.yoy <- c(1, rep(0, length(SS) - 1))             # creat vector indentifying YOY
        N1 = p.yoy * pop                                  # Calc YOY pop
        N1plus = SS * (1 - p.yoy)/ sum(SS * p.yoy) * pop  # Calc age 1+ pop
        N <- N1 + N1plus                                  # Sum YOY and 1+ pop by age
    }
    
    # Pop from age 1 onward (age 1 to Tmax)
    if(pop.stage == "age1+") {
        SS <- stable.stage(pmx)                           # Stable stage distn
        p.1plus <- c(0, rep(1, length(SS) - 1))           # Create vectore identifying non YOYs
        N1plus <- SS * p.1plus / sum(SS * p.1plus) * pop  # Calc age 1+ pop
        N1 <- SS * (1 - p.1plus)/ sum(SS * p.1plus) * pop # calc YOY pop
        N <- N1 + N1plus                                  # sum YOY and 1+ pop by age
    }
      
    # Output pop size by age class
    N
}

#===============================================================================
# Function to calc. number of fish is desired stage when another stage pop is
# known

stage_pop <- function(pmx,         # projection matrix
                      stage,       # stage you want pop for (YOY, Juv, Adt)
                      pop,         # desired pop size
                      pop.stage,   # stage pop repsents - character ["total", "adult","juvenile","yoy","age1+"]
                      p.rep        # Proportion adult by age class              
) {
  
  # needs popbio library
  require(popbio)
  
  # age-specific pop
  N <- initPop(pmx, pop, pop.stage, p.rep)
  if(length(N) > length(p.rep)){
    N <- colSums(matrix(N,ncol = 5, byrow=T))
  }
  
  if(stage %in% c("YOY", "yoy")) {
    N <- ((1-p.rep) * N)[1]
  } else if(stage %in% c("Juv", "J", "Juvenile")) {
    N <- sum(((1-p.rep) * N)[-1])
  } else if(stage %in% c("Adt", "A", "Adult")) {
    N <- sum(p.rep * N)
  }
  (1-p.rep) 
  
  # Output pop size by age class
  N
}


#===============================================================================
# MVP - Function to estimate minimum viable population size for a given
# probability of persistence from parameter values for exticntion probability vs
# adults population size power curve

MVP <- function(P, a, b){
  mvp <- exp(log(P/a)/b)
  names(mvp) <- paste0("P[",P,"]")
  mvp
}

#===============================================================================
# Area per individual
# function to estimate area per indiviuals based on allometric relationships 
# of lakes and river (Randall et al. 1995, Minns 2003)
# Lakes: API = exp(-10.37)*TL^2.58 in mm
# River: API = exp(-13.28)*L^2.904 in mm

API <- function(L, habitat) {
  if(habitat %in% c("Lake", "lake", "L")) {
    api <-exp(-10.37)*L^2.58  # Lake (mm)
  }else if (habitat %in% c("River", "river", "R")){
    api <-  exp(-13.28)*L^2.904 #  River (mm)
  }
  api
}

#===============================================================================
# MAPV 
# Function to estimate the daily area requirement by a population to determin
# when maximum MAPV is. Combines growth and survival with MVP to est. area
# requirements per day of life. Specific to redisde dace in that survival
# switches from Lorrenzen to constant at age-2

# Growth estimated from length at age and API estiamte using Minns
# API = exp(-10.37) * L ^ 2.58 (Lake)
# API = exp(-13.28) * L ^ 2.904 (River)
MAPVd_cohort_l <- function(MVP, age.min, age.max, params, parAPI = c(exp(-13.28), 2.904)) {
  
  params$m0 <- params$m0.1
  
  # LENGTH
  t <- (age.min * 365):((age.max+1) * 365) / 365 # calc age in years for earch day
  t <- t[-length(t)]                             # 
  Ld <- with(params, VB.growth(t, Linf, K, t0))  # calc daily length
  names(Ld) = NULL
  
  # SURVIVAL
  sd <- with(params, (Ld[1:(length(Ld) -1)] * exp(-K * 1/365)/Ld[2:length(Ld)])^(m0 / K / Linf))
  
  # AREA
  APId <- parAPI[1] * Ld^parAPI[2] # Daily API
  MVPd <- MVP * cumprod(sd)        # Daily MVP - inital MVP (pop at age.min with total adult pop at MVP) time cumulative survival
  
  MAPV <- APId * c(MVP, MVPd)                     # Daily MAPV - daily API times daily MVP
  as.data.frame(list(jDay = rep(1:365, age.max + 1), Age = t, 
                     MAPV = MAPV/10000, API = APId)) # Output in hectares
  
}

# Growth estimated from weight at age and API estiamte using Tandall et al. (1995) density (inverse)
# D = 30199.52*W^-1.01 (Lake)
# D = 78432.82*W^-0.94  (river)
MAPVd_cohort_w <- function(MVP, age.min, age.max, M.e = 0, data, parD = c(78432.82, -0.94)) {
  
  # Egg mortality
  MVP <- MVP * exp(-M.e)
  
  # LENGTH
  t <- (age.min * 365):((age.max) * 365) / 365 # calc age in years for earch day
  t <- t[-length(t)]                             # 
  Ld <- with(data, VB.growth(t, Linf, K, L0 = H0))  # calc daily length
  Wd <- with(data, (a.LW * 1000) * Ld ^ b.LW)
  names(Ld) = NULL
  
  # SURVIVAL
  
  # YOY survival
  # Est daily relationship over first year
  # assume M inversly related to length (Lorenzen 2000); solve for m0.yoy to give
  # correct syoy value; Optimize survival relationship
  # Optimization function
  
  m0.yoy_optim = function (m0.0, target, data){
    s.yoy <- with(data, (H0 * exp(-K)/VB.growth(1, Linf, K, L0 = H0))^(m0.0/K/Linf))
    (-log(s.yoy) - target)^2
  }
  
  # Run optimization
  data$m0.yoy <-  optimize(m0.yoy_optim, c(1, 1000), tol = 1e-16, 
                      target = -log(data$syoy.1d) - M.e, data = data)$minimum
  
  # Juvenile survival
  # calc m0 for juvenile survival
  data$m0 <- with(data, VB.growth(c.MO, Linf, K, L0 = H0) * M)
  
  # Survival schedule over lifetime (daily)
  sd <- sapply(1:(length(t)-1), function(x) {
    if(t[[x]] <= 1){
      with(data, (Ld[x] * exp(-K * 1/365)/Ld[x + 1]) ^ (m0.yoy / K / Linf))
    } else if(t[[x]] < data$c.MO){
      with(data, (Ld[x] * exp(-K * 1/365)/Ld[x + 1]) ^ (m0 / K / Linf))
    } else {
      exp(-data$M/365)
    }
  })
  
  # Sex ratio
  # Proportion Male - COSEWIC 2006
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
  R <- R[-1]
  
  # AREA
  APId <- 1/(parD[1] * Wd ^ parD[2])
  MVPd.F <- MVP * cumprod(sd)        # Daily MVP - inital MVP (pop at age.min with total adult pop at MVP) time cumulative survival
  MVPd.M <- MVPd.F * R
  
  MAPV.F <- APId * c(MVP, MVPd.F)      # Daily MAPV - daily API times daily MVP in hectares
  MAPV.M <- APId * c(MVP, MVPd.M)      # Daily MAPV - daily API times daily MVP in hectares
  
  as.data.frame(list(jDay = rep(1:365, age.max), Age = t,  L = Ld, Survival = cumprod(c(1,sd)), R = c(1,R), API = APId, 
                     MAPV.F =  MAPV.F,  MAPV.M =  MAPV.M,  MAPV.T = MAPV.F + MAPV.M,
                     Pop.F = c(MVP, MVPd.F), Pop.M = c(MVP, MVPd.M), Pop.T = c(MVP * 2, MVPd.F + MVPd.M)
  ))
}

# Stage-specific daily area requirements
# aggregate MAPV_day by day of year within each stage
MAPV_jDay_func <- function(y, data){
  
  # assing ages to each value
  y$Age_year <- trunc(y$Age)
  
  # Total area usage - all Ages
  MAPV_T <- with(y, aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_T <- with(y, aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # YOY stage - only age 0
  MAPV_Y <- with(y[y$Age_year == 0,], aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_Y <- with(y[y$Age_year == 0,], aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # Age 1+ - all ages excluding YOY
  MAPV_1Plus <- with(y[y$Age_year >= 1,], aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_1Plus <- with(y[y$Age_year >= 1,], aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # Juvenile Stage 
  y.j <- y[y$Age_year >= 1 & y$Age_year < data$c.MO,]
  MAPV_J <- with(y.j, aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_J <- with(y.j, aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # Adult stage 
  y.a <- y[y$Age_year >= data$c.MO,]
  MAPV_A <- with(y.a, aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_A <- with(y.a, aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # > 1000 mm
  y.1000 <- y[y$L >= 1000,]
  MAPV_1000 <- with(y.1000, aggregate(MAPV.T, by = list(jDay), sum)) # Area
  Pop_1000 <- with(y.1000, aggregate(Pop.T, by = list(jDay), sum))   # Pop size
  
  # Combine into data frame
  MAPV <- rbind(MAPV_Y, MAPV_J, MAPV_A, MAPV_1000, MAPV_1Plus, MAPV_T)
  colnames(MAPV) <- c("jDay", "MAPV")
  MAPV$Pop <- c(Pop_Y[, 2], Pop_J[, 2], Pop_A[, 2], Pop_1000[, 2], Pop_1Plus[, 2], Pop_T[, 2])
  MAPV$Stage = c(rep("YOY", 365), rep("Juvenile", 365), rep("Adult", 365), rep("L1000", 365),
                 rep("Age-1+", 365), rep("All Ages", 365))
  MAPV$Stage = factor(MAPV$Stage, levels = c("YOY", "Juvenile", "Adult", "L1000", "Age-1+", "All Ages"))
  
  # Output
  MAPV
  
}


# Function to estimate the space requirements of a life stage daily. Used stage-
# specific stable stage dist'n and tracks daily growth and pop size to est.
# daily weighted mean length and calc. daily API and MAPV (ha). Take the max value
# to get max stage-specific MAPV.

# Growth estimated from weight at age and API estiamte using Tandall et al. (1995) density (inverse)
# D = 30199.52*W^-1.01 (Lake)
# D = 78432.82*W^-0.94  (river)
MAPVd_stage_w <- function(MVP, age.min, age.max, data, parD = c(78432.82, -0.94)) {
  
  # LENGTH
  t <- (age.min * 365):((age.max) * 365) / 365 # calc age in years for earch day
  t <- t[-length(t)]                             # 
  Ld <- with(data, VB.growth(t, Linf, K, L0 = H0))  # calc daily length
  Wd <- with(data, (a.LW * 1000) * Ld ^ b.LW)
  names(Ld) = NULL
  
  # SURVIVAL
  
  # YOY survival
  # Est daily relationship over first year
  # assume M inversly related to length (Lorenzen 2000); solve for m0.yoy to give
  # correct syoy value; Optimize survival relationship
  # Optimization function
  
  m0.yoy_optim = function (m0.0, target, data){
    s.yoy <- with(data, (H0 * exp(-K)/VB.growth(1, Linf, K, L0 = H0))^(m0.0/K/Linf))
    (-log(s.yoy) - target)^2
  }
  
  # Run optimization
  data$m0.yoy <-  optimize(m0.yoy_optim, c(1, 1000), tol = 1e-16, 
                           target = -log(data$syoy.1), data = data)$minimum
  
  # Juvenile survival
  # calc m0 for juvenile survival
  data$m0 <- with(data, VB.growth(c.MO, Linf, K, L0 = H0) * M)
  
  # Survival schedule over lifetime (daily)
  sd <- sapply(1:(length(t)-1), function(x) {
    if(t[[x]] <= 1){
      with(data, (Ld[x] * exp(-K * 1/365)/Ld[x + 1]) ^ (m0.yoy / K / Linf))
    } else if(t[[x]] < data$c.MO){
      with(data, (Ld[x] * exp(-K * 1/365)/Ld[x + 1]) ^ (m0 / K / Linf))
    } else {
      exp(-data$M/365)
    }
  })
  
  Nd <- MVP * cumprod(sd)        # Daily MVP - inital MVP (pop at age.min with total adult pop at MVP) time cumulative survival
  
  # Matrix of daily mean age-specific weight and pop size - rows: ages (0:100) - cols: day of year
  Wd.mat <- matrix(Wd, ncol = 365, byrow = TRUE)
  Nd.mat <- matrix(c(Nd, 0) , ncol = 365, byrow = TRUE)
  
  # Mean daily stage-specific weights 
  # weighted mean of daily weights by population size within specificed stage
  
  # YOY
  W.YOY <- Wd.mat[1, ]
  
  # JUVENILE - age 1 to c.MO
  W.Juv  <- sapply(1:365, function(x){
    sum(Wd.mat[, x][2:(data$c.MO - 1)] * Nd.mat[, x][2:(data$c.MO - 1)]) / sum(Nd.mat[, x][2:(data$c.MO - 1)])
  })
  
  # Adult - age c.MO to Tmax
  W.Adt  <- sapply(1:365, function(x){
    sum(Wd.mat[, x][data$c.MO:data$Tmax] * Nd.mat[, x][data$c.MO:data$Tmax]) / sum(Nd.mat[, x][data$c.MO:data$Tmax])
  })
  
  # AGE 1+ - age 1 to Tmax
  W.1plus  <- sapply(1:365, function(x){
    sum(Wd.mat[, x][-1] * Nd.mat[, x][-1]) / sum(Nd.mat[, x][-1])
  })
  
  # All ages
  W.All <- sapply(1:365, function(x){
    sum(Wd.mat[, x] * Nd.mat[, x]) / sum(Nd.mat[, x])
  })
  
  # Population size - stage specific
  # YOY
  N.YOY <- Nd.mat[1 ,]                                  # YOY
  N.Juv <- apply(Nd.mat[2:(data$c.MO - 1), ], 2, sum)   # JUVENILE - age 1 to c.MO
  N.Adt <- apply(Nd.mat[data$c.MO:data$Tmax, ], 2, sum) # Adult - age c.MO to Tmax
  N.1plus <- apply(Nd.mat[-1, ], 2, sum)                # AGE 1+ - age 1 to Tmax
  N.All <- apply(Nd.mat, 2, sum)                        # All ages
  
  # Density
  D.YOY <- (parD[1] * W.YOY ^ parD[2])     # YOY
  D.Juv <- (parD[1] * W.Juv ^ parD[2])     # JUVENILE - age 1 to c.MO
  D.Adt <- (parD[1] * W.Adt ^ parD[2])     # Adult - age c.MO to Tmax
  D.1plus <- (parD[1] * W.1plus ^ parD[2]) # AGE 1+ - age 1 to Tmax
  D.All <- (parD[1] * W.All ^ parD[2])     # All ages
  
  # MAPV - Area useage - Population/density
  MAPV.YOY <- N.YOY / D.YOY       # YOY
  MAPV.Juv <- N.Juv / D.Juv       # JUVENILE - age 1 to c.MO
  MAPV.Adt <- N.Adt / D.Adt       # Adult - age c.MO to Tmax
  MAPV.1plus <- N.1plus / D.1plus # AGE 1+ - age 1 to Tmax
  MAPV.All <- N.All / D.All       # All ages
  
  # Output
  list(
    N = list(
      YOY = N.YOY,
      Juv = N.Juv,
      Adt = N.Adt,
      onePlus = N.1plus,
      All = N.All
    ),
    Weight = list(
      YOY = W.YOY,
      Juv = W.Juv,
      Adt = W.Adt,
      onePlus = W.1plus,
      All = W.All
    ),
    Density = list(
      YOY = D.YOY,
      Juv = D.Juv,
      Adt = D.Adt,
      onePlus = D.1plus,
      All = D.All
    ),
    MAPV = list(
      YOY = MAPV.YOY,
      Juv = MAPV.Juv,
      Adt = MAPV.Adt,
      onePlus = MAPV.1plus,
      All = MAPV.All
    ))
}

#===============================================================================
# Transient Allowable Harm

# Function to simulate the impact of harm on population growth rate
# outputs a comparison of lambda values between an unharmed and harm populations
T_harm_impact <- function(harm,                 # proportional change to survival
                          mx = A.w,             # Projection matrix expression
                          numReps = 1000,       # number of simulation replicates
                          years = 10,           # number of years to simulate
                          lambda,               # pop growth rate
                          cor.Li = cor_Li,      # Correlation struction for growth increments
                          names.Li = names_Li,  # names for growth increments
                          dist.Li = dist_Li,    # distributions names for random length increments
                          params.Li = params_Li,# distributions parameters fore random length increments
                          cor.s = cor_s,        # Correlation struction for vital rates
                          names.s = names_s,    # names for vital rates
                          dist.s = dist_s,      # distributions names for random  vital rates
                          params.s = params_s,  # distributions parameters fore random vital rate
                          L.mean = L_mean,      # Initial (mean) length at age
                          data = data)          # Life history data 
{
  
  yearNum = years + 10
  # Create random matrices for unharmed simulation
  pmx_rand <- pmx_build_runs(mx, numReps, yearNum, lambda, 
                             cor.Li, names.Li,dist.Li, params.Li,
                             cor.s, names.s, dist.s, params.s,
                             L.mean, data, sensitivity = FALSE)
  ss_list_init <- pmx_rand$ss
  Fs_list_init <- pmx_rand$Fs
  Gs_list_init <- pmx_rand$Gs

  # clean up
  rm(pmx_rand)
  
  # Apply harm - est. mean pop growth rate
  lambda_harm <- lapply(harm, function(h){
    
    L <- sapply(1:numReps, function(n){
      
      # Extract vital rates for rep
      ss <- ss_list_init[[n]]
      Fs <- Fs_list_init[[n]]
      Gs <- Gs_list_init[[n]]
      
      # reduce survival in first matrix  by harm
      ss[[11]] <- unlist(ss[[11]]) * (1 - h)
      
      # evaluate harm matrces
      pmx_list_harm <- lapply(1:yearNum, function(y){
        pmx_eval(mx, c(ss[[y]],  Fs[[y]], Gs[[y]]))
      })
      
      # Inital population stabel structure
      #pop.harm <- list(stable.stage(pmx_list_harm[[1]]))
      
      # Annual pop size
      # loops through years winin reps
      #for(j in 2:yearNum) {
      #  pop.harm[j] <- list(pmx_list_harm[[j-1]] %*% pop.harm[[j-1]])      # harm pop
      #}
      
      # Calc pop growth rate for each replicate over years
      #pop.harm <- sapply(pop.harm, sum)
      #(pop.harm[yearNum] / pop.harm[11]) ^ (1/(yearNum - 11))
      #geom_mean(pop.harm[12:yearNum] / pop.harm[11:19])
      geom_mean(lambda.real(stable.stage(pmx_list_harm[[1]]), pmx_list_harm)$annual[11:yearNum])
      
    })
  })
}


#===============================================================================
# ALLOWABLE HARM


# Function to estimate allowable ahrm from vital rate elasticities
allowable_harm <- function(Elast, Lambda){
  
  Ah <- 1/Elast * ((1 - Lambda) / Lambda)
  Ah
}

# Function to simulate the impact of harm on population growth rate
# outputs a comparison of lambda values between an unharmed and harm populations
A_harm_impact <- function(s_means,              # mean unharmed survival rates
                          harm,                 # proportional change to survival
                          mx = A.w,             # Projection matrix expression
                          numReps = 500,        # number of simulation replicates
                          years = 100,          # number of years to simulate
                          N0 = N0,              # stage population vector
                          lambda,               # pop grwoth rate
                          cor.Li = cor_Li_ind,  # Correlation struction for growth increments
                          names.Li = names_Li,  # names for growth increments
                          dist.Li = dist_Li,    # distributions names for random length increments
                          params.Li = params_Li,# distributions parameters fore random length increments
                          cor.s = cor_s,        # Correlation struction for vital rates
                          names.s = names_s,    # names for vital rates
                          dist.s = dist_s,      # distributions names for random  vital rates
                          L.mean = L_mean,      # Initial (mean) length at age
                          data = lh_data)       # Life history data 
{
  
  # HARM 
  s_means <- s_means + s_means * harm
  
  M_means <- {
    Ms <- -log(s_means)
    names(Ms) <- paste0("M", 0:(length(s_means) - 1))
    Ms
  }
  
  # Est beta shape paramters
  beta.M <- beta_stretch_val(M_means, sd = (M_means * data$cv.M))
  
  # Vital rates
  params.s <- lapply(1:data$Tmax, function(i) beta_stretch_val(M_means[[i]], sd = (M_means[[i]] * data$cv.M)))

  # Create random matrices
  pmx_rand <- pmx_build_runs(mx, numReps, years, lambda,
                             cor.Li ,names.Li, dist.Li, params.Li,
                             cor.s, names.s, dist.s, params.s,
                             L.mean, data, sensitivity = FALSE)
  pmx_list <- pmx_rand$pmx
  
  # clean up
  rm(pmx_rand, s_means, M_means, beta.M, params.s)
  
  # Annual lambda alues
  lambdas_year_new <- sapply(unlist(pmx_list, recursive = F), lambda)
  
  # est. geom_mean of lambda of pmx.list
  lambdas_mean_new <- sapply(pmx_list, lambda_mean)
  lambda_new <- mean(lambdas_mean_new)
  
  # est. geom_mean of lambda of pmx.list
  lambdas_mean_new <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$annual) 
  lambdas_mean_new <- sapply(1:length(pmx_list), function(i) lambda.real(N0, pmx_list[[i]])$longterm) 
  lambda_new <- mean(lambdas_mean_new)
  
  # OUTPUT
  list(
    lambda_year = lambdas_year_new, # Annual lambda values
    lambda_mean = lambdas_mean_new, # Mean lambda values over 100 years simulatinos
    lambda = lambda_new             # Mean lambda value over replications
    )
}


#===============================================================================
# Function to projection a population 
# Includes stochasticity and catastrophe
# For use in MVP and recovery time analysis 
# function included in rep_proj function

# Function to project population size including density dependence effects
Projection <- function(pmx_list,    # List of random projection matrices of length years
                       N0,          # initial pop size as stage-structure vector
                       years,       # years to run simulation
                       p.rep,       # Proportion mature by age class
                       p.cat,       # Probability of castastrophe
                       e.cat,       # effect on castastrophe on pop size (percent recuction)
                       gen.time    # Generation Time
){
  
  # needs popbio library
  require(popbio)

  # Extinction
  Ext = 0
  
  # set starting pop
  N = N0
  
  # Castastrophes
  Catastrophe = sample(c(1,0),years,replace = TRUE,
                       prob = c(p.cat / gen.time, 1 - p.cat  /gen.time))
  
  # initialize output vectors
  Nvec <- rep(NA,years + 1); Nvec[1] <- sum(N * p.rep)
  lambdas <- rep(NA,years)
  
  # loop through years
  for(t in 1:years){
    
    if(Catastrophe[t] == 1){
      N <- N * e.cat
      Na <-  sum(N * p.rep)
      Nvec[t+1] <- Na         # Number of mature fish in pop
      lambdas[t] <- e.cat # pop growth rate
      next
    } 
    # Projection matrix
    A <- pmx_list[[t]]  #  
    
    # project the population 1 year ahead.
    N <- as.vector(A %*% N )                   
    
    if(Catastrophe[t] == 1){
      N <- N * e.cat
    }
    
    # Adult pop size
    Na <-  sum(N * p.rep)
    
    Nvec[t+1] <- Na         # Number of mature fish in pop
    lambdas[t] <- lambda(A) # pop growth rate
 
  }
  # output
  list(adult.pop = as.data.frame(list(year = 0:years,
                                      Cat = c(0,Catastrophe),
                                      N = Nvec)),
       lambdas = lambdas)
}

#===============================================================================
# Function that replicates Projection 
# replicates simulations numReps times and model runs maxRuns times
# Begins by generating random matricesthen executes Projection function
# outputs a list (maxRuns in length) of matrices (years by numReps) of adult
# population sizes
# also indicates when extinction took place for 

rep_proj <- function(mx,           # Projection matrix expression
                     N0,           # initial pop size
                     maxRuns,      # number of runs of simlulation (new matrices est each run)  
                     years,        # number of years to simulate
                     numReps,      # number of simulation replicates
                     lambda,       # pop growth rate
                     cor.Li,       # Correlation struction for growth increments
                     names.Li,     # names for growth increments
                     dist.Li,      # distributions names for random length increments
                     params.Li,    # distributions parameters fore random length increments
                     cor.s,        # Correlation struction for vital rates
                     names.s,      # names for vital rates
                     dist.s,       # distributions names for random  vital rates
                     params.s,     # distributions parameters fore random vital rate
                     L.mean,       # Initial (mean) length at age
                     p.rep,        # Proportion mature by age class
                     p.cat,        # Probability of castastrophe
                     e.cat,        # effect on castastrophe on pop size (percent recuction)
                     gen.time,     # Generation Time
                     data          # Life histoty data
) {
  

  # Loop though simuation
  adtPop <- lapply(1:maxRuns, function(i){
    
    # generate matrices with random variability
    pmx_rand <- pmx_build_runs(mx, numReps, years, lambda, 
                               cor.Li, names.Li, dist.Li, params.Li, 
                               cor.s, names.s, dist.s, params.s,
                               L.mean, data)
    pmx_rand <- pmx_rand$pmx
    
    # clean up
    gc()
    
    # loop through replicates
    Na <- lapply(1:numReps, function(j) {
    
      proj_pop <- Projection(pmx_list = pmx_rand[[j]], # List of vandom variables to build matrix
                             N0 = N0,                  # initial pop size
                             years = years,            # years to run simulation
                             p.rep = p.rep,            # Proportion mature by age class
                             p.cat = p.cat,            # Probability of castastrophe
                             e.cat = e.cat,            # effect on castastrophe on pop size (percent recuction)
                             gen.time = gen.time)      # Generation Time
      
      # extract adult pop
      Na <- proj_pop$adult.pop$N
      
      # clean up
      rm(proj_pop)
      gc()
      
      # output
      Na
      
      })
    # clean up
    rm(pmx_rand)
    gc()
    
    # output
    Na
    })

  # output
  adtPop
}

#===============================================================================
# Function used in MVP analysis
# used output from rep_proj
# identifies if population is below extinction threshols threshold for each run

# Counts number of extirpations (NA < threshold)
is_ext <- function(x, threshold, years){
  ifelse(sum(x[1:years] <= threshold) > 1, 1, 0)
}

# Applies is_ext to a list of MVP outputs
is_ext_list <- function(x, threshold, years) {
  sapply(x,  is_ext, threshold = threshold, years = years)
}

# Takes MVP sim output and organizes into table; solves MVP power function and
# Est. MVP values fore 5% and 1% persistence
MVP_solve <- function(MVP_res,   # MVP Simluation output
                      init_pop,  # Initial adult pop size in simluation
                      threshold, # Extinction threshold 
                      Cat,       # Catastrophe rates
                      DU,        # Designatable unit
                      numReps,   # Number of simluation replicates 
                      maxRuns,   # Number of independent mode runs
                      start,     # list of inital values fore nls fit
                      fit_range) # survival probs over which to fit relationship
  {    
  
  # MVP data
  data <-  as.data.frame(list(
    Cat = Cat,
    Threshold = threshold,
    DU = DU,
    N = sort(rep(init_pop, maxRuns)),
    P = unlist(lapply(MVP_res, function(x) sapply(x, is_ext_list, threshold = threshold))) / numReps
  ))
  
  # MVP model
  MVP_mod <- nls(P ~ a * N^b, data = data[data$P <= fit_range[2] & data$P >= fit_range[1],], 
                 start = start)
  
  # MVP value
  MVP_val <- MVP(c(0.05, 0.01), coef(MVP_mod)[[1]], coef(MVP_mod)[[2]])
  
  # Output
  list(data = data,
       model = MVP_mod,
       MVP = MVP_val)
  
}

#===============================================================================
# Function used in recovery time analysis
# used output from rep_proj
# counts the number of times the population is above a threshold each year for
# each run
prop_rec <- function(x, threshold){
     length(which(x > threshold)) 
}

prop_rec_list <- function(x, threshold) {
      apply(x, 1, prop_rec, threshold = threshold)
}

#===============================================================================
